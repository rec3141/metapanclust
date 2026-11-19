"""Simulated bins for one focal species A, scored against A's pangenome (built without the held-out genomes).

Bins (base = held-out A genome x, whole or reduced to 70% of its length):
  clean        x alone
  strain       base + pseudo-contigs of another held-out A genome at 5/10/20/40% of x's length
  genus:<sp>   base + pseudo-contigs of a genome of congeneric species <sp> at the same fractions
  far          base + pseudo-contigs of a Flavobacterium psychrophilum genome (other phylum)
  piece_*      one 50-kb pseudo-contig offered to x at 70%: own removed part, strain, genus:<sp>, far
Pseudo-contigs: contigs cut into <= 50-kb windows, genes kept with their window.

Gene -> family: best mmseqs hit of the gene against A's family representatives;
strict = >= 80% identity and >= 80% coverage of both (PPanGGOLiN default), relaxed = >= 30% / 50%.
A gene with no strict hit is foreign. Scores as in experiments/binning_sim/02_simulate.py:
comp, dup (strict), dup_rel (relaxed), foreign fraction, S_tot surprise in bits
(present -log2 p, absent -log2(1 - p comp), foreign gene -log2 1/(n+2); p = (k+1)/(n+2)).
Merge excess = (S(bin+piece) - S(bin)) / Mb added - S(bin)/Mb; reject if > 0.
excess_content: the same with S_content = present + foreign terms only (no absence term).
usage: 03_simulate.py WORKDIR FOCAL_TAG OUT_TSV [STRICT_IDENTITY, default 0.8]
  env GENOMES: genome table (default WORKDIR/genomes.tsv); in a plan table (grex/) every row,
  contaminants included, carries its focal tag
"""
import math, os, random, sys
from collections import Counter, defaultdict
import numpy as np, pandas as pd

FRAG = 50_000
FRACS = [0.05, 0.10, 0.20, 0.40]
REPS = 3
NPIECE = 10
W, A, OUT = sys.argv[1:4]
IDENT = float(sys.argv[4]) if len(sys.argv) > 4 else 0.8   # strict family cut (identity)
GENOMES = os.environ.get('GENOMES', f'{W}/genomes.tsv')
rng = random.Random(11)

G = pd.read_csv(GENOMES, sep='\t', header=None, usecols=range(5), names=['acc', 'species', 'genus', 'role', 'tag'],
                dtype=str).fillna('')
foc = G[G.tag == A]
species_A = foc.species.iloc[0]
genus_A = foc.genus.iloc[0]
train = foc[foc.role == 'train'].acc.tolist()
n = len(train)

# ---- pangenome ----
pan = f'{W}/pan/{A}'
fam_genome_copies = defaultdict(Counter)
for line in open(f'{pan}/gene_families.tsv'):
    fam, gene = line.split('\t')[:2]
    fam_genome_copies[fam][gene.rsplit('_CDS_', 1)[0]] += 1
fams = list(fam_genome_copies)
fam_idx = {f: i for i, f in enumerate(fams)}
k = np.array([len(fam_genome_copies[f]) for f in fams])
p_arr = (k + 1) / (n + 2)
kmulti = {f: sum(1 for c in fam_genome_copies[f].values() if c > 1) for f in fams}
persistent = set(l.strip() for l in open(f'{pan}/partitions/persistent.txt'))
markers = [f for f in persistent if kmulti[f] / n < 0.05]
marker_idx = np.array([fam_idx[f] for f in markers])
bits_foreign = -math.log2(1 / (n + 2))

# ---- genes of non-training genomes, with strict / relaxed family ----
hits = pd.read_csv(f'{W}/hits/{A}.m8', sep='\t', header=None,
                   names=['q', 'fam', 'pid', 'qcov', 'tcov', 'bits'])
hits = hits.sort_values('bits', ascending=False).drop_duplicates('q')
pid = hits.pid / (100 if hits.pid.max() > 1.0 else 1)
hits['strict'] = (pid >= IDENT) & (hits.qcov >= 0.8) & (hits.tcov >= 0.8)
rel = dict(zip(hits.q, hits.fam))
strict = dict(zip(hits.q[hits.strict], hits.fam[hits.strict]))

if foc.role.isin(['contam', 'far']).any():   # plan mode: every row carries its focal tag
    use = foc[foc.role != 'train']
else:                                        # local mode: contaminants are untagged rows of the genus
    use = G[(G.role != 'train') & ((G.genus == genus_A) | (G.role == 'far'))]
genomes = {}
for r in use.itertuples():
    genes = []
    for line in open(f'{W}/prot/{r.acc}.faa'):
        if line[0] == '>':
            q = line[1:].strip()
            g, contig, start = q.rsplit('|', 2)
            genes.append((f'{contig}:{int(start) // FRAG}', fam_idx.get(strict.get(q), -1),
                          fam_idx.get(rel.get(q), -1)))
    L, name = {}, None
    for line in open(f'{W}/genomes/{r.acc[3:]}.fna'):
        if line[0] == '>':
            name = line[1:].split()[0]; L[name] = 0
        else:
            L[name] += len(line) - 1
    pieces = {f'{c}:{i // FRAG}': min(FRAG, ln - i) for c, ln in L.items() for i in range(0, ln, FRAG)}
    gp = defaultdict(list)
    for pc, fs, fr in genes:
        gp[pc].append((fs, fr))
    if r.species == species_A:
        src = 'self'
    elif r.role == 'far':
        src = 'far'
    else:
        src = 'genus:' + r.species.split(' ', 1)[1]
    genomes[r.acc] = dict(src=src, role=r.role, gp=gp, pieces=pieces, length=sum(L.values()),
                          strict_frac=np.mean([fs >= 0 for _, fs, _ in genes]))


def score(parts):
    """parts: list of (genome, piece list) -> score dict"""
    cnt = np.zeros(len(fams), int); cnt_rel = np.zeros(len(fams), int)
    ng = nfor = 0; ln = 0
    for g, keep in parts:
        d = genomes[g]
        for pc in keep:
            ln += d['pieces'][pc]
            for fs, fr in d['gp'].get(pc, ()):
                ng += 1
                if fs >= 0:
                    cnt[fs] += 1
                else:
                    nfor += 1
                if fr >= 0:
                    cnt_rel[fr] += 1
    pres = cnt > 0
    mc = cnt[marker_idx]
    comp = float((mc > 0).mean())
    dup = float(np.clip(mc - 1, 0, None).sum() / len(markers))
    dup_rel = float(np.clip(cnt_rel[marker_idx] - 1, 0, None).sum() / len(markers))
    S_content = float(-np.log2(p_arr[pres]).sum() + bits_foreign * nfor)
    S = S_content + float(-np.log2(1 - p_arr[~pres] * comp).sum())
    return dict(mb=ln / 1e6, genes=ng, comp=comp, dup=dup, dup_rel=dup_rel,
                foreign=nfor / max(ng, 1), S_tot=S, S_perMb=S / (ln / 1e6), S_content=S_content)


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


rows = []


def emit(kind, base_kind, frac, rep, x, y, base, add):
    b = score([(x, base)])
    m = score([(x, base), (y, add)])
    add_mb = m['mb'] - b['mb']
    r = dict(A=A, kind=kind, base=base_kind, frac=frac, rep=rep, genome=x, other=y, half=half[x], ident=IDENT, **m,
             base_perMb=b['S_perMb'], base_comp=b['comp'], add_mb=add_mb,
             d_S_perMbAdded=(m['S_tot'] - b['S_tot']) / add_mb,
             d_foreign=m['foreign'] - b['foreign'], d_dup=m['dup'] - b['dup'], d_dup_rel=m['dup_rel'] - b['dup_rel'])
    r['excess'] = r['d_S_perMbAdded'] - r['base_perMb']
    # variant without the absence term: content surprise (present + foreign) per Mb
    r['excess_content'] = (m['S_content'] - b['S_content']) / add_mb - b['S_content'] / b['mb']
    rows.append(r)


held = [g for g, d in genomes.items() if d['src'] == 'self' and d['role'] == 'held']
if len(held) >= 6:
    held = sorted(held)
# with >= 6 held-out genomes, split them into a calibration and a test half; strain partners are
# drawn within the half, so a threshold set on calibration merges can be tested on the other half
half = {g: ('calib' if i < len(held) // 2 else 'test') if len(held) >= 6 else 'all' for i, g in enumerate(held)}
by_src = defaultdict(list)
for g, d in genomes.items():
    if d['src'] != 'self':
        by_src[d['src']].append(g)
for x in held:
    full = list(genomes[x]['pieces'])
    rows.append(dict(A=A, kind='clean', base='full', frac=0, rep=0, genome=x, other='', half=half[x], ident=IDENT, **score([(x, full)])))
    for rep in range(REPS):
        b70 = remove(x, 0.3)
        sources = [('strain', [h for h in held if h != x and half[h] == half[x]])] + sorted(by_src.items())
        for kind, pool in sources:
            y = rng.choice(pool)
            for f in FRACS:
                for bk, base in (('full', full), ('70', b70)):
                    emit(kind, bk, f, rep, x, y, base, draw(y, f * genomes[x]['length']))
    # single-contig moves onto the 70% base
    b70 = remove(x, 0.3); gone = [pc for pc in genomes[x]['pieces'] if pc not in set(b70)]
    pools = [('piece_own', [(x, pc) for pc in gone]),
             ('piece_strain', [(h, pc) for h in held if h != x for pc in genomes[h]['pieces']])]
    pools += [(f'piece_{s}', [(g, pc) for g in gs for pc in genomes[g]['pieces']]) for s, gs in sorted(by_src.items())]
    for kind, pool in pools:
        pool = [(g, pc) for g, pc in pool if genomes[g]['pieces'][pc] >= 20_000 and genomes[g]['gp'].get(pc)]
        for rep, (g, pc) in enumerate(rng.sample(pool, min(NPIECE, len(pool)))):
            emit(kind, '70', round(genomes[g]['pieces'][pc] / 1000, 1), rep, x, g, b70, [pc])

df = pd.DataFrame(rows)
df.to_csv(OUT, sep='\t', index=False)
sf = pd.DataFrame([(g, d['src'], d['strict_frac']) for g, d in genomes.items()], columns=['genome', 'src', 'strict_frac'])
sf.to_csv(OUT.replace('.tsv', '_strictfrac.tsv'), sep='\t', index=False)
print(f'{A}: n={n} families={len(fams)} markers={len(markers)} held={len(held)} bins={len(df)}')
