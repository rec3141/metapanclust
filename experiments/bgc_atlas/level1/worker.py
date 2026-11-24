"""Level-1 test on BGC gene cluster families (BGC Atlas R2.0) within one GTDB species.

usage: worker.py PANGBANK_ID H5 ATLAS_TSV OUTDIR [--max-genomes 400] [--seed 1]

Genomes: the Atlas genomes of the species (GTDB collection, MIMAG medium/high) that are also in the PanGBank
pangenome (accessions with the GB_/RS_ prefix stripped), at most --max-genomes.
Units (targets): primary families carried by 2 .. N-2 of those genomes.
Leakage guard: PanGBank gene families with a gene overlapping any BGC region of the species (same contig name,
overlapping coordinates) are removed from the gene-content distance, and each target family is removed from the
BGC-content distance used to predict it.

Predictors (10-fold CV over genomes, same folds for all; bits = -log2 P(observed state)):
  prev            Beta(1,1) prevalence of the training genomes
  ker_gc, tree_gc gene-content relatedness (Jaccard on accessory PanGBank families, BGC-overlapping families
                  removed): kernel (bandwidth by training LOO) and gain/loss tree on UPGMA
  ker_bgc, tree_bgc  BGC-content relatedness (Jaccard on the other families)
Detection (tip model of the trees, and how predictions are turned into P(observed)):
  *_c    P(observed | present) = CheckM2 completeness
  *_frag  completeness x P(a region of the family's median length lies intact on one of the genome's contigs)
          (metapanclust.detect.p_intact with the genome's PanGBank contig lengths), relative to a 1 kb gene. Strict:
          wrong for the Atlas, where contig-edge fragments are still called and assigned; kept as a contrast
  *_seen  1 - P(region entirely lost) (metapanclust.detect.p_lost)
Writes OUTDIR/<id>.targets.tsv (per target: carriers, tier, completeness of member BGCs, bits per predictor) and
OUTDIR/<id>.genomes.tsv (per genome: source type, completeness, N50, bits per predictor).
"""
import argparse, os, sys, time
import numpy as np, pandas as pd

sys.path.insert(0, '/data/dev/metapanclust')
from metapanclust.tree import upgma, tree_posterior, jaccard_dist
from metapanclust.detect import p_intact, p_lost

ap = argparse.ArgumentParser()
ap.add_argument('pid'); ap.add_argument('h5'); ap.add_argument('atlas'); ap.add_argument('out')
ap.add_argument('--max-genomes', type=int, default=400); ap.add_argument('--seed', type=int, default=1)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
T0 = time.time()
EPS = 0.002
QS = np.array([0.03, 0.1, 0.3, 1.0, 3.0])


def log(m): print(f'[{time.time() - T0:6.0f}s] {m}', flush=True)


def bits(p, y):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return -(y * np.log2(p) + (1 - y) * np.log2(1 - p))


# ---------------------------------------------------------------- data
import tables
H = tables.open_file(a.h5)
A = H.root.annotations
pb_genomes = np.array([x.decode() for x in A.genomes.col('name')])
at = pd.read_csv(a.atlas, sep='\t')
gq = at.drop_duplicates('genome').set_index('genome')
pool = [g for g in pb_genomes if g in gq.index]
rng = np.random.default_rng(a.seed)
if len(pool) > a.max_genomes:
    pool = list(rng.choice(pool, a.max_genomes, replace=False))
pool = sorted(pool)
N = len(pool)
gidx = {g: i for i, g in enumerate(pool)}
log(f'{a.pid}: {len(pb_genomes)} PanGBank genomes, {gq.shape[0]} Atlas genomes, {N} in both')
if N < 20:
    H.close(); log('fewer than 20 shared genomes; skipped'); sys.exit(0)
ctg_id = A.contigs.col('ID'); ctg_gen = [x.decode() for x in A.contigs.col('genome')]
ctg_name = [x.decode() for x in A.contigs.col('name')]
ctg_len = A.contigs.col('length')
ctg = pd.DataFrame({'cid': ctg_id, 'genome': ctg_gen, 'name': ctg_name, 'len': ctg_len})
ctg = ctg[ctg.genome.isin(gidx)]
genes = pd.DataFrame({'gene': A.genes.col('ID'), 'cid': A.genes.col('contig'), 'gd': A.genes.col('genedata_id')})
genes = genes[genes.cid.isin(set(ctg.cid))]
gd = A.genedata
gdid = gd.col('genedata_id'); st = np.zeros(gdid.max() + 1, np.int64); sp = np.zeros_like(st)
st[gdid] = gd.col('start'); sp[gdid] = gd.col('stop')
genes['start'] = st[genes.gd.values]; genes['stop'] = sp[genes.gd.values]
gf = pd.DataFrame({'gene': H.root.geneFamilies.col('gene'), 'fam': H.root.geneFamilies.col('geneFam')})
genes = genes.merge(gf, on='gene').merge(ctg[['cid', 'genome', 'name']], on='cid')
FI = H.root.geneFamiliesInfo
fams = np.array([x.decode() for x in FI.col('name')]); part = np.array([x.decode()[:1] for x in FI.col('partition')])
H.close()
fpos = {f.encode(): i for i, f in enumerate(fams)}
genes['f'] = genes.fam.map(fpos)
genes['g'] = genes.genome.map(gidx)

# BGC regions of the pool genomes; mask PanGBank families overlapping them
reg = at[at.genome.isin(gidx) & at.bgc_id.notna()].copy()
reg['g'] = reg.genome.map(gidx)
names_pb = set(ctg.name)
hit_names = reg.contig.isin(names_pb).mean() if len(reg) else 0
mask_f = set()
for (g, c), r in reg.groupby(['genome', 'contig']):
    gg = genes[(genes.genome == g) & (genes.name == c)]
    if not len(gg): continue
    s_, e_ = gg.start.values, gg.stop.values
    for a0, b0 in zip(r.start_pos.values, r.end_pos.values):
        mask_f.update(gg.f.values[(s_ <= b0) & (e_ >= a0)].tolist())
log(f'{len(reg)} BGC regions; contig names matched for {hit_names:.1%}; {len(mask_f)} PanGBank families masked')
R = np.zeros((len(fams), N), bool); R[genes.f.values, genes.g.values] = True
acc = (part != 'P') & (R.sum(1) >= 2) & (R.sum(1) <= N - 1)
acc[list(mask_f)] = False
D_gc = jaccard_dist(R[acc].T)

# family presence (primary family) and median region length per family
reg_f = reg.dropna(subset=['gcf_id'])
fam_ids = sorted(reg_f.gcf_id.unique())
fi = {f: i for i, f in enumerate(fam_ids)}
Y_all = np.zeros((N, len(fam_ids)), bool)
Y_all[reg_f.g.values, reg_f.gcf_id.map(fi).values] = True
k = Y_all.sum(0)
tg = np.flatnonzero((k >= 2) & (k <= N - 2))
log(f'{len(fam_ids)} families in the pool, {len(tg)} targets; gene-content features {acc.sum()}')
if not len(tg):
    sys.exit(0)
span = reg_f.groupby('gcf_id').region_len.median().reindex(fam_ids).values
comp = np.clip(np.array([gq.checkm_completeness.get(g, 100) for g in pool]) / 100, 0.3, 0.99)
src = np.array([gq.source_type.get(g, '') for g in pool])
L_by_g = ctg.groupby('genome').len.apply(lambda v: np.sort(v.values)).to_dict()
n50 = np.array([np.sort(L_by_g[g])[::-1][np.searchsorted(np.cumsum(np.sort(L_by_g[g])[::-1]), L_by_g[g].sum() / 2)]
                for g in pool])
# fragmentation-aware detection per (genome, target): completeness x intact(span)/intact(1 kb)
# two detection models per (genome, target):
#   intact: completeness x P(region intact on one contig) relative to a 1 kb gene (strict: a fragment = absent)
#   seen:   1 - P(region entirely lost) (a fragment still counts, as antiSMASH calls and the Atlas assigns
#           contig-edge regions)
DET = np.zeros((N, len(tg))); SEEN = np.zeros((N, len(tg)))
for i, g in enumerate(pool):
    L = L_by_g[g].astype(float)
    f0 = p_intact(1000, L, 1.0)
    DET[i] = np.clip(comp[i] * p_intact(np.maximum(span[tg], 1000), L, 1.0) / max(f0, 1e-9), 0.01, comp[i])
    SEEN[i] = np.clip(1 - p_lost(np.maximum(span[tg], 1000), L, comp[i]), 0.01, 0.999)
Y = Y_all[:, tg].astype(float)


# ---------------------------------------------------------------- predictors
def kernel(D, Y, tr, te, p0):
    med = np.median(D[np.triu_indices(len(D), 1)]); med = med if med > 0 else max(D.max(), 1e-6)
    best = None
    for q in QS:
        W = np.exp(-D / (q * med))
        Wtt = W[np.ix_(tr, tr)].copy(); np.fill_diagonal(Wtt, 0)
        loo = (Wtt @ Y[tr] + p0) / (Wtt.sum(1, keepdims=True) + 1)
        b = bits(loo, Y[tr]).mean(0)
        Wte = W[np.ix_(te, tr)]; pte = (Wte @ Y[tr] + p0) / (Wte.sum(1, keepdims=True) + 1)
        if best is None: best = [b, pte]
        else:
            m = b < best[0]; best[0] = np.where(m, b, best[0]); best[1][:, m] = pte[:, m]
    return best[1]


def tree_det(tree, Y, tr, det):
    """tree posterior with a per-(genome,target) detection probability: run per target column groups by
    using the tip likelihoods directly (tree_posterior takes one completeness per tip, so call per target)"""
    out = np.zeros(Y.shape)
    for t in range(Y.shape[1]):
        out[:, t] = tree_posterior(tree, Y[:, [t]], tr, det[:, t])[:, 0]
    return out


fold = rng.permutation(np.arange(N) % 10)
tree_gc = upgma(D_gc)
P = {m: np.zeros((N, len(tg))) for m in ['prev', 'ker_gc', 'tree_gc_c', 'tree_gc_frag', 'tree_gc_seen', 'ker_bgc',
                                          'tree_bgc_c', 'tree_bgc_frag', 'tree_bgc_seen']}
COMP = np.repeat(comp[:, None], len(tg), 1)
# BGC-content distance without the target: computed per target (few families per species)
Db_t = [jaccard_dist(np.delete(Y_all, tg[t], axis=1)) for t in range(len(tg))]
tb_t = [upgma(D) for D in Db_t]
for fo in range(10):
    te = fold == fo; tr = ~te
    p0 = (Y[tr].sum(0) + 1) / (tr.sum() + 2)
    P['prev'][te] = p0
    P['ker_gc'][te] = kernel(D_gc, Y, tr, te, p0)
    P['tree_gc_c'][te] = tree_det(tree_gc, Y, tr, COMP)[te]
    P['tree_gc_frag'][te] = tree_det(tree_gc, Y, tr, DET)[te]
    P['tree_gc_seen'][te] = tree_det(tree_gc, Y, tr, SEEN)[te]
    for t in range(len(tg)):
        Db, tb = Db_t[t], tb_t[t]
        yt = Y[:, [t]]
        P['ker_bgc'][te, t] = kernel(Db, yt, tr, te, p0[[t]])[:, 0]
        P['tree_bgc_c'][te, t] = tree_posterior(tb, yt, tr, COMP[:, t])[te, 0]
        P['tree_bgc_frag'][te, t] = tree_posterior(tb, yt, tr, DET[:, t])[te, 0]
        P['tree_bgc_seen'][te, t] = tree_posterior(tb, yt, tr, SEEN[:, t])[te, 0]
    log(f'fold {fo + 1}/10')
# turn P(present) into P(observed): prevalence and kernels predict observed state directly; trees predict presence
Q = {}
for m, p in P.items():
    if m.startswith('tree'):
        d = DET if m.endswith('frag') else (SEEN if m.endswith('seen') else COMP)
        Q[m] = p * d + (1 - p) * EPS
    else:
        Q[m] = p
B = {m: bits(q, Y) for m, q in Q.items()}
tier = pd.Series(['F' if '-F-' in str(x) else 'S' for x in fam_ids])
cplt = reg_f.assign(c=reg_f.is_complete.astype(str).isin(['t', 'True', 'true'])).groupby('gcf_id').c.mean().reindex(fam_ids).values
pd.DataFrame(dict(id=a.pid, gcf=[fam_ids[i] for i in tg], carriers=Y.sum(0).astype(int), span=span[tg],
                  frac_complete=np.round(cplt[tg], 3), **{f'bits_{m}': B[m].mean(0).round(5) for m in B})
             ).to_csv(f'{a.out}/{a.pid}.targets.tsv', sep='\t', index=False)
pd.DataFrame(dict(id=a.pid, genome=pool, source=src, completeness=comp.round(3), n50=n50,
                  **{f'bits_{m}': B[m].mean(1).round(5) for m in B})
             ).to_csv(f'{a.out}/{a.pid}.genomes.tsv', sep='\t', index=False)
log(' '.join(f'{m}={B[m].mean():.4f}' for m in B))
