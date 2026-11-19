"""Simulated bins scored against a species pangenome built without them.

For species A (pangenome rebuilt without 6 held-out genomes) and an unrelated species B:
  clean      held-out A genome alone
  incomplete A with 10-50% of its length removed (random pseudo-contigs)
  contam_B   A plus pseudo-contigs of a held-out B genome at 5/10/20/40% of A's length
  strain     A plus pseudo-contigs of another held-out A genome at the same fractions
  strain_inc A at 70% plus pseudo-contigs of another held-out A genome (complementary merge)
Pseudo-contigs: real contigs cut into <= FRAG bp windows, genes assigned by start position, so
complete genomes can be subsampled like MAGs.

Scores per bin (families from `ppanggolin projection`, strict = 80% id / 80% cov, relaxed = 30%/50%):
  comp            fraction of single-copy-marker persistent families present (CheckM-like)
  dup / dup_rel   extra copies over those markers / markers (strict / relaxed homology)
  foreign         fraction of genes in no pangenome family
  S_*             surprise in bits: present -log2 p, absent -log2(1 - p c) (c = comp, or 1 for
                  S_abs1), foreign gene -log2 1/(n+2), multicopy -log2 q (q = share of carriers
                  that are multicopy); p = (k+1)/(n+2) over the n training genomes.
usage: 02_simulate.py WORKDIR SPECIES_A SPECIES_B OUT_TSV
  SPECIES_B = the other species, or near_<A> for the related Chesterfield bins in held_near_<A>.tsv
"""
import csv, glob, gzip, math, os, random, sys
from collections import Counter, defaultdict
import numpy as np, pandas as pd

FRAG = 50_000
FRACS_ADD = [0.05, 0.10, 0.20, 0.40]
FRACS_DEL = [0.10, 0.20, 0.30, 0.40, 0.50]
REPS = 5
W, A, B, OUT = sys.argv[1:5]
rng = random.Random(1)


def read_lengths(fna):
    L, name = {}, None
    op = gzip.open if fna.endswith('.gz') else open
    with op(fna, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].split()[0]; L[name] = 0
            else:
                L[name] += len(line.strip())
    return L


# ---- pangenome of A (training genomes only) ----
pan = f'{W}/pan_{A}'
train = [l.split('\t')[0] for l in open(f'{W}/train_{A}.tsv')]
n = len(train)
fam_genome_copies = defaultdict(Counter)       # family -> genome -> copies
for line in open(f'{pan}/gene_families.tsv'):
    fam, gene = line.split('\t')[:2]
    g = gene.rsplit('_CDS_', 1)[0]
    fam_genome_copies[fam][g] += 1
fams = list(fam_genome_copies)
k = {f: len(fam_genome_copies[f]) for f in fams}
p = {f: (k[f] + 1) / (n + 2) for f in fams}
kmulti = {f: sum(1 for c in fam_genome_copies[f].values() if c > 1) for f in fams}
q = {f: (kmulti[f] + 1) / (k[f] + 2) for f in fams}
persistent = set(l.strip() for l in open(f'{pan}/partitions/persistent.txt'))
# single-copy markers: persistent and multicopy in < 5% of training genomes (PPanGGOLiN dup_margin)
markers = {f for f in persistent if kmulti[f] / n < 0.05}
p_new = 1 / (n + 2)
# modules
mod_fams = defaultdict(set)
for r in csv.DictReader(open(f'{pan}/functional_modules.tsv'), delimiter='\t'):
    mod_fams[r['module_id']].add(r['family_id'])
mod_k = Counter()
for r in csv.DictReader(open(f'{pan}/modules_in_genomes.tsv'), delimiter='\t'):
    if float(r['completion']) >= 0.5:
        mod_k[r['module_id']] += 1
mod_p = {m: (mod_k[m] + 1) / (n + 2) for m in mod_fams}
fam_arr = np.array(fams)
p_arr = np.array([p[f] for f in fams])
fam_idx = {f: i for i, f in enumerate(fams)}


# ---- projected genomes: genes with (pseudo-contig, family strict, family relaxed) ----
def load_proj(projdir, genome):
    fs = glob.glob(f'{projdir}/{genome}/*.tsv') + glob.glob(f'{projdir}/{genome}*.tsv')
    fs = [f for f in fs if 'summary' not in os.path.basename(f) and 'statistic' not in os.path.basename(f)]
    tab = [f for f in fs if 'table' in f or os.path.basename(f) == f'{genome}.tsv']
    df = pd.read_csv(tab[0] if tab else fs[0], sep='\t', dtype=str)
    return df


held = {}
for s in (A, B):
    for line in open(f'{W}/held_{s}.tsv'):
        g, fna = line.rstrip('\n').split('\t')[:2]
        held[g] = (s, fna)
genomes = {}
for g, (s, fna) in held.items():
    L = read_lengths(fna)
    tag = 'near' if s.startswith('near') else ''
    st = load_proj(f'{W}/proj{tag}_{A}', g)
    rl = load_proj(f'{W}/proj{tag}rel_{A}', g)
    rl = rl.set_index('gene')['family']
    st['family_rel'] = st['gene'].map(rl)
    st['start'] = st['start'].astype(int)
    st['piece'] = st['contig'] + ':' + (st['start'] // FRAG).astype(str)
    pieces = {}
    for c, ln in L.items():
        for i in range(0, ln, FRAG):
            pieces[f'{c}:{i // FRAG}'] = min(FRAG, ln - i)
    genomes[g] = dict(species=s, genes=st[['piece', 'family', 'family_rel']], pieces=pieces,
                      length=sum(L.values()))


def score(gene_df, length):
    fam_cnt = Counter(gene_df['family'])
    rel_cnt = Counter(gene_df['family_rel'].dropna())
    ng = len(gene_df)
    known = {f: c for f, c in fam_cnt.items() if f in fam_idx}
    nforeign = sum(c for f, c in fam_cnt.items() if f not in fam_idx)
    comp = sum(1 for f in markers if f in known) / len(markers)
    dup = sum(known[f] - 1 for f in markers if f in known) / len(markers)
    dup_rel = sum(rel_cnt[f] - 1 for f in markers if rel_cnt.get(f, 0) > 1) / len(markers)
    pres = np.zeros(len(fams), bool)
    for f in known:
        pres[fam_idx[f]] = True
    S_pres = float(-np.log2(p_arr[pres]).sum())
    S_abs = float(-np.log2(1 - p_arr[~pres] * comp).sum())
    S_abs1 = float(-np.log2(np.clip(1 - p_arr[~pres], 1e-12, 1)).sum())
    S_for = -math.log2(p_new) * nforeign
    S_dup = sum(-math.log2(q[f]) for f, c in known.items() if c > 1)
    S_mod = 0.0
    for m, mf in mod_fams.items():
        present = sum(1 for f in mf if f in known) / len(mf) >= 0.5
        S_mod += -math.log2(mod_p[m]) if present else -math.log2(1 - mod_p[m] * comp)
    mb = length / 1e6
    return dict(mb=mb, genes=ng, comp=comp, dup=dup, dup_rel=dup_rel, foreign=nforeign / ng,
                S_pres=S_pres, S_abs=S_abs, S_abs1=S_abs1, S_for=S_for, S_dup=S_dup, S_mod=S_mod,
                S_tot=S_pres + S_abs + S_for, S_totdup=S_pres + S_abs + S_for + S_dup)


def subset(g, keep):
    d = genomes[g]
    return d['genes'][d['genes']['piece'].isin(keep)], sum(d['pieces'][x] for x in keep)


def draw(g, target_len):
    """random pseudo-contigs of g until their length reaches target_len"""
    pcs = list(genomes[g]['pieces']); rng.shuffle(pcs)
    out, tot = [], 0
    for x in pcs:
        if tot >= target_len:
            break
        out.append(x); tot += genomes[g]['pieces'][x]
    return out


def remove(g, frac):
    pcs = list(genomes[g]['pieces']); rng.shuffle(pcs)
    target = frac * genomes[g]['length']; gone, tot = set(), 0
    for x in pcs:
        if tot >= target:
            break
        gone.add(x); tot += genomes[g]['pieces'][x]
    return [x for x in genomes[g]['pieces'] if x not in gone]


rows = []
As = [g for g in genomes if genomes[g]['species'] == A]
Bs = [g for g in genomes if genomes[g]['species'] == B]


def emit(kind, frac, rep, x, other, parts):
    """parts: list of (genome, pieces); first is the base bin, the rest is the added piece"""
    dfs, ln = [], 0
    for g, keep in parts:
        d, l = subset(g, keep); dfs.append(d); ln += l
    merged = score(pd.concat(dfs), ln)
    r = dict(A=A, B=B, kind=kind, frac=frac, rep=rep, genome=x, other=other, **merged)
    if len(parts) > 1:                         # base alone and added piece alone, for merge deltas
        base = score(dfs[0], subset(*parts[0])[1])
        add_mb = ln / 1e6 - base['mb']
        for key in ('S_tot', 'S_totdup', 'S_mod', 'S_abs1'):
            r[f'base_{key}_perMb'] = base[key] / base['mb']
            r[f'd_{key}_perMbAdded'] = (merged[key] - base[key]) / add_mb
        r['base_mb'] = base['mb']; r['base_comp'] = base['comp']; r['base_dup'] = base['dup']
    rows.append(r)


for x in As:
    allp = list(genomes[x]['pieces'])
    emit('clean', 0, 0, x, '', [(x, allp)])
    for rep in range(REPS):
        for f in FRACS_DEL:
            emit('incomplete', f, rep, x, '', [(x, remove(x, f))])
        y = rng.choice(Bs)
        for f in FRACS_ADD:
            emit('contam_B', f, rep, x, y, [(x, allp), (y, draw(y, f * genomes[x]['length']))])
        y = rng.choice([a for a in As if a != x])
        for f in FRACS_ADD:
            emit('strain', f, rep, x, y, [(x, allp), (y, draw(y, f * genomes[x]['length']))])
        base70 = remove(x, 0.3)
        for f in FRACS_ADD:
            emit('strain_inc', f, rep, x, y, [(x, base70), (y, draw(y, f * genomes[x]['length']))])
        yb = rng.choice(Bs)
        for f in FRACS_ADD:
            emit('contam_B_inc', f, rep, x, yb, [(x, base70), (yb, draw(yb, f * genomes[x]['length']))])

# single-contig moves: one pseudo-contig offered to the bin (x at 70%), from x's own removed part,
# from another A genome, or from B; frac column holds the contig length in kb
NPIECE = 15
for x in As:
    keep = remove(x, 0.3); gone = [pc for pc in genomes[x]['pieces'] if pc not in set(keep)]
    for src_kind, pool in (('piece_own', [(x, pc) for pc in gone]),
                           ('piece_strain', [(y, pc) for y in As if y != x for pc in genomes[y]['pieces']]),
                           ('piece_B', [(y, pc) for y in Bs for pc in genomes[y]['pieces']])):
        for rep, (y, pc) in enumerate(rng.sample(pool, min(NPIECE, len(pool)))):
            if subset(y, [pc])[0].empty:
                continue
            emit(src_kind, round(genomes[y]['pieces'][pc] / 1000, 1), rep, x, y, [(x, keep), (y, [pc])])

df = pd.DataFrame(rows)
df.to_csv(OUT, sep='\t', index=False)
print(f'{A} vs {B}: training genomes {n}, families {len(fams)}, persistent {len(persistent)}, '
      f'markers {len(markers)}, modules {len(mod_fams)}; {len(df)} simulated bins -> {OUT}')
