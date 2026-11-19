"""Identity likelihood-ratio merge test (docs/theory/merge_test.md eq. 8) vs surprise-per-Mb excess.

Rebuilds the binning_sim_genus merges with the same RNG seed and draw order as
experiments/binning_sim_genus/03_simulate.py, so the bins are the same ones; each merge is scored
on the *added* genes with:
  excess     surprise-per-Mb rule (recomputed, strict cut 0.80)
  lr         mean per-gene log LR, eq. (8): H0 = this species' held-out strains (g_s, mu_s fitted on
             held-out genomes other than the base and the partner), H1 = congener at 95% ANI
             (mu_A = max(kappa*0.05, 1.5 mu_s), g_A = g_s - 1.6*0.05), alpha = 2, kappa = 0.65
  lrc        the same with divergence clipped at XC = 0.10 (mu's refitted on clipped values)
  xclip/xmed mean of clipped divergence / median divergence of added genes with a family hit
  xbar       mean aa divergence (1 - identity) of added genes that hit a focal family (UMP part)
  for95      fraction of added genes without a family at the 95%-identity hard cut
  for80      the same at 80% (the current rule's foreign count)
Identity is to the best-hit family representative (the existing mmseqs tables), not the nearest member.
Also writes per-genome along-contig autocorrelation / variance inflation of divergence.
usage: 01_simulate.py BSG_WORKDIR FOCAL OUTDIR
"""
import math, random, sys
from collections import Counter, defaultdict
import numpy as np, pandas as pd

FRAG = 50_000
FRACS = [0.05, 0.10, 0.20, 0.40]
REPS = 3
NPIECE = 10
W, A, OUT = sys.argv[1:4]
IDENT = 0.8
ALPHA, KAPPA = 2.0, 0.65
XC = 0.10   # divergence clip for the robust variants (paralog / distant best hits)
rng = random.Random(11)

G = pd.read_csv(f'{W}/genomes.tsv', sep='\t', header=None, usecols=range(5),
                names=['acc', 'species', 'genus', 'role', 'tag'], dtype=str).fillna('')
foc = G[G.tag == A]
species_A = foc.species.iloc[0]; genus_A = foc.genus.iloc[0]
train = foc[foc.role == 'train'].acc.tolist(); n = len(train)

pan = f'{W}/pan/{A}'
fam_genome_copies = defaultdict(Counter)
for line in open(f'{pan}/gene_families.tsv'):
    fam, gene = line.split('\t')[:2]
    fam_genome_copies[fam][gene.rsplit('_CDS_', 1)[0]] += 1
fams = list(fam_genome_copies); fam_idx = {f: i for i, f in enumerate(fams)}
k = np.array([len(fam_genome_copies[f]) for f in fams]); p_arr = (k + 1) / (n + 2)
kmulti = {f: sum(1 for c in fam_genome_copies[f].values() if c > 1) for f in fams}
persistent = set(l.strip() for l in open(f'{pan}/partitions/persistent.txt'))
markers = [f for f in persistent if kmulti[f] / n < 0.05]
marker_idx = np.array([fam_idx[f] for f in markers])
bits_foreign = -math.log2(1 / (n + 2))

hits = pd.read_csv(f'{W}/hits/{A}.m8', sep='\t', header=None, names=['q', 'fam', 'pid', 'qcov', 'tcov', 'bits'])
hits = hits.sort_values('bits', ascending=False).drop_duplicates('q')
hits['pid'] = hits.pid / (100 if hits.pid.max() > 1.0 else 1)
hits['strict'] = (hits.pid >= IDENT) & (hits.qcov >= 0.8) & (hits.tcov >= 0.8)
hits['s95'] = (hits.pid >= 0.95) & (hits.qcov >= 0.8) & (hits.tcov >= 0.8)
rel = dict(zip(hits.q, hits.fam)); strict = dict(zip(hits.q[hits.strict], hits.fam[hits.strict]))
s95 = set(hits.q[hits.s95]); divg = dict(zip(hits.q, 1 - hits.pid))

use = G[(G.role != 'train') & ((G.genus == genus_A) | (G.role == 'far'))]
genomes = {}
autoc = []
for r in use.itertuples():
    genes = []
    for line in open(f'{W}/prot/{r.acc}.faa'):
        if line[0] == '>':
            q = line[1:].strip()
            g, contig, start = q.rsplit('|', 2)
            # gene: (strict fam, relaxed fam, divergence or nan, in95)
            genes.append((contig, int(start), fam_idx.get(strict.get(q), -1), fam_idx.get(rel.get(q), -1),
                          divg.get(q, np.nan), q in s95))
    L, name = {}, None
    for line in open(f'{W}/genomes/{r.acc[3:]}.fna'):
        if line[0] == '>':
            name = line[1:].split()[0]; L[name] = 0
        else:
            L[name] += len(line) - 1
    pieces = {f'{c}:{i // FRAG}': min(FRAG, ln - i) for c, ln in L.items() for i in range(0, ln, FRAG)}
    gp = defaultdict(list)
    for c, st, fs, fr, x, i95 in genes:
        gp[f'{c}:{st // FRAG}'].append((fs, fr, x, i95))
    if r.species == species_A:
        src = 'self'
    elif r.role == 'far':
        src = 'far'
    else:
        src = 'genus:' + r.species.split(' ', 1)[1]
    genomes[r.acc] = dict(src=src, role=r.role, gp=gp, pieces=pieces, length=sum(L.values()))
    # ---- clustering of divergence along contigs ----
    if src == 'far':
        continue
    gs = sorted(genes, key=lambda t: (t[0], t[1]))
    xs = defaultdict(list); fo = defaultdict(list)
    for c, st, fs, fr, x, i95 in gs:
        fo[c].append(0.0 if fr >= 0 else 1.0)
        if fr >= 0:
            xs[c].append(x)
    allx = np.concatenate([np.array(v) for v in xs.values()]); mx, vx = allx.mean(), allx.var()
    allf = np.concatenate([np.array(v) for v in fo.values()]); mf, vf = allf.mean(), allf.var()
    row = dict(A=A, genome=r.acc, src=src, n_shared=len(allx), mean_x=mx, sd_x=allx.std(), foreign=mf)
    for lag in (1, 2, 5, 10, 20):
        num = den = 0.0
        for v in xs.values():
            v = np.array(v)
            if len(v) > lag:
                num += ((v[:-lag] - mx) * (v[lag:] - mx)).sum(); den += len(v) - lag
        row[f'acf_x_{lag}'] = num / den / vx if den else np.nan
    # variance inflation of block means: phi = Var(block mean) / (Var(x)/block size)
    for B in (10, 25, 50):
        bm = []; bf = []
        for c in xs:
            v = np.array(xs[c])
            for i in range(0, len(v) - B + 1, B):
                bm.append(v[i:i + B].mean())
            v = np.array(fo[c])
            for i in range(0, len(v) - B + 1, B):
                bf.append(v[i:i + B].mean())
        row[f'phi_x_{B}'] = np.var(bm) / (vx / B) if len(bm) > 5 else np.nan
        row[f'phi_for_{B}'] = np.var(bf) / (vf / B) if len(bf) > 5 and vf > 0 else np.nan
    autoc.append(row)

held = sorted(g for g, d in genomes.items() if d['src'] == 'self' and d['role'] == 'held') \
    if sum(1 for d in genomes.values() if d['src'] == 'self' and d['role'] == 'held') >= 6 else \
    [g for g, d in genomes.items() if d['src'] == 'self' and d['role'] == 'held']


def h0_params(exclude):
    """g_s, mu_s from held-out strains other than `exclude`"""
    tot = shared = 0; sx = sc = 0.0
    for g in held:
        if g in exclude:
            continue
        for lst in genomes[g]['gp'].values():
            for fs, fr, x, i95 in lst:
                tot += 1
                if fr >= 0:
                    shared += 1; sx += x; sc += min(x, XC)
    return shared / tot, sx / shared, sc / shared


H0 = {}


def lr_coefs(x, y):
    key = frozenset((x, y))
    if key not in H0:
        gs, mus, muc = h0_params(key)
        muA = max(KAPPA * 0.05, 1.5 * mus); gA = max(gs - 1.6 * 0.05, 0.05)
        mucA = max(KAPPA * 0.05, 1.5 * muc)
        H0[key] = (math.log((1 - gA) / (1 - gs)), math.log(gA / gs) + ALPHA * math.log(mus / muA),
                   ALPHA * (1 / mus - 1 / muA), gs, mus,
                   math.log(gA / gs) + ALPHA * math.log(muc / mucA), ALPHA * (1 / muc - 1 / mucA), muc)
    return H0[key]


def score(parts):
    cnt = np.zeros(len(fams), int); cnt_rel = np.zeros(len(fams), int)
    ng = nfor = 0; ln = 0
    for g, keep in parts:
        d = genomes[g]
        for pc in keep:
            ln += d['pieces'][pc]
            for fs, fr, x, i95 in d['gp'].get(pc, ()):
                ng += 1
                if fs >= 0:
                    cnt[fs] += 1
                else:
                    nfor += 1
                if fr >= 0:
                    cnt_rel[fr] += 1
    pres = cnt > 0
    mc = cnt[marker_idx]; comp = float((mc > 0).mean())
    S = float(-np.log2(p_arr[pres]).sum() + bits_foreign * nfor) + float(-np.log2(1 - p_arr[~pres] * comp).sum())
    return dict(mb=ln / 1e6, S_tot=S, S_perMb=S / (ln / 1e6))


def idstats(g, keep, coefs):
    c0, c1, c2 = coefs[:3]; c1c, c2c = coefs[5:7]
    m0 = m = 0; T = Tc = 0.0; n95 = n80 = 0; tot = 0; xs = []
    for pc in keep:
        for fs, fr, x, i95 in genomes[g]['gp'].get(pc, ()):
            tot += 1
            if fr >= 0:
                m += 1; T += x; Tc += min(x, XC); xs.append(x)
            else:
                m0 += 1
            n95 += not i95; n80 += fs < 0
    lam = m0 * c0 + m * c1 + c2 * T
    lamc = m0 * c0 + m * c1c + c2c * Tc
    return dict(add_genes=tot, lr=lam / max(tot, 1), lr_sum=lam, lrc=lamc / max(tot, 1), xbar=T / max(m, 1),
                xclip=Tc / max(m, 1), xmed=float(np.median(xs)) if xs else np.nan, for95=n95 / max(tot, 1),
                for80=n80 / max(tot, 1), m_noham=m0)


rows = []


def emit(kind, base_kind, frac, rep, x, y, base, add):
    b = score([(x, base)]); mm = score([(x, base), (y, add)])
    add_mb = mm['mb'] - b['mb']
    co = lr_coefs(x, y)
    r = dict(A=A, kind=kind, base=base_kind, frac=frac, rep=rep, genome=x, other=y, add_mb=add_mb,
             excess=(mm['S_tot'] - b['S_tot']) / add_mb - b['S_perMb'], g_s=co[3], mu_s=co[4], mu_s_clip=co[7], **idstats(y, add, co))
    rows.append(r)


half = {g: ('calib' if i < len(held) // 2 else 'test') if len(held) >= 6 else 'all' for i, g in enumerate(held)}
by_src = defaultdict(list)
for g, d in genomes.items():
    if d['src'] != 'self':
        by_src[d['src']].append(g)


def draw(g, target):
    pcs = list(genomes[g]['pieces']); rng.shuffle(pcs)
    out, tot = [], 0
    for x in pcs:
        if tot >= target:
            break
        out.append(x); tot += genomes[g]['pieces'][x]
    return out


def remove(g, frac):
    pcs = list(genomes[g]['pieces']); rng.shuffle(pcs)
    gone, tot = set(), 0
    for x in pcs:
        if tot >= frac * genomes[g]['length']:
            break
        gone.add(x); tot += genomes[g]['pieces'][x]
    return [x for x in genomes[g]['pieces'] if x not in gone]


for x in held:
    full = list(genomes[x]['pieces'])
    for rep in range(REPS):
        b70 = remove(x, 0.3)
        sources = [('strain', [h for h in held if h != x and half[h] == half[x]])] + sorted(by_src.items())
        for kind, pool in sources:
            y = rng.choice(pool)
            for f in FRACS:
                for bk, base in (('full', full), ('70', b70)):
                    emit(kind, bk, f, rep, x, y, base, draw(y, f * genomes[x]['length']))
    b70 = remove(x, 0.3); gone = [pc for pc in genomes[x]['pieces'] if pc not in set(b70)]
    pools = [('piece_own', [(x, pc) for pc in gone]),
             ('piece_strain', [(h, pc) for h in held if h != x for pc in genomes[h]['pieces']])]
    pools += [(f'piece_{s}', [(g, pc) for g in gs for pc in genomes[g]['pieces']]) for s, gs in sorted(by_src.items())]
    for kind, pool in pools:
        pool = [(g, pc) for g, pc in pool if genomes[g]['pieces'][pc] >= 20_000 and genomes[g]['gp'].get(pc)]
        for rep, (g, pc) in enumerate(rng.sample(pool, min(NPIECE, len(pool)))):
            emit(kind, '70', round(genomes[g]['pieces'][pc] / 1000, 1), rep, x, g, b70, [pc])

pd.DataFrame(rows).to_csv(f'{OUT}/sim_{A}.tsv', sep='\t', index=False)
pd.DataFrame(autoc).to_csv(f'{OUT}/autocorr_{A}.tsv', sep='\t', index=False)
print(A, len(rows), 'bins')
