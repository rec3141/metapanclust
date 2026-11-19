"""Group conserved loci into cluster families by content similarity, with Leiden instead of connected components.

Two loci A, B are compared by matched genes: a gene of A is matched if B has a gene of the same MMseqs2
family or a Spacedust hit partner (Spacedust's search is more sensitive than 30%-identity families, so
distant genera still match). similarity = (mA + mB) / (nA + nB) (BiG-SCAPE-like). An edge needs
similarity >= SMIN and >= MIN_SHARED matched genes on both sides, so one shared hub gene (ABC ATPase,
response regulator) cannot link two loci. Leiden (CPM, weights = similarity) cuts the graph into groups
whose internal mean similarity exceeds RES; no single-linkage chaining. Each family is then split again
at RES_SUB into subfamilies (e.g. peptide ABC transporter paralogues, lineage variants of nuo).
usage: 03_families.py GENES_TSV LOCI_TSV BIN_TAXONOMY CLUSTERS_TSV OUT_PREFIX [smin=0.4] [res=0.1] [min_shared=2] [seed=1] [best=0.9] [res_sub=0.4]
-> OUT_PREFIX_catalogue.tsv (one row per cluster family), OUT_PREFIX_members.tsv (one row per locus)"""
import sys
from collections import Counter, defaultdict
import numpy as np, pandas as pd, scipy.sparse as sp
import igraph as ig, leidenalg as la

genes_tsv, loci_tsv, tax_tsv, sd_tsv, out = sys.argv[1:6]
SMIN = float(sys.argv[6]) if len(sys.argv) > 6 else 0.4
RES = float(sys.argv[7]) if len(sys.argv) > 7 else 0.1
MIN_SHARED = int(sys.argv[8]) if len(sys.argv) > 8 else 2
SEED = int(sys.argv[9]) if len(sys.argv) > 9 else 1
BEST = float(sys.argv[10]) if len(sys.argv) > 10 else 0.9    # 0 = keep all edges
RES_SUB = float(sys.argv[11]) if len(sys.argv) > 11 else 0.4

G = pd.read_csv(genes_tsv, sep='\t', keep_default_na=False, dtype={'idx': int})
L = pd.read_csv(loci_tsv, sep='\t', keep_default_na=False)
genus = {}
for l in open(tax_tsv):
    b, t = l.rstrip('\n').split('\t')[:2]; r = t.split(';'); genus[b] = r[5] if len(r) > 5 and r[5] != 'g__' else b
gi = G.set_index(['contig', 'idx'])
fam_of = gi.family.to_dict(); prod_of = gi['product'].to_dict(); strand_of = gi.strand.to_dict()

loc_genes, loc_fams, loc_order = [], [], []
gene_locus = {}
for k, r in enumerate(L.itertuples()):
    gs = [(r.contig, i) for i in range(r.first, r.last + 1)]
    for g in gs: gene_locus[g] = k
    fs = [fam_of[g] for g in gs]; st = [strand_of[g] for g in gs]
    if st.count('-') > st.count('+'):                  # orient so most genes read left to right
        fs = fs[::-1]; st = ['+' if s == '-' else '-' for s in st[::-1]]; gs = gs[::-1]
    loc_genes.append(gs)
    loc_order.append([f if s == '+' else f + '(-)' for f, s in zip(fs, st)])
    loc_fams.append(set(fs))
n_loc = len(L)

# Spacedust hit partners between loci: (a, b) -> (matched genes of a, matched genes of b), a < b
hits = defaultdict(lambda: (set(), set()))
for l in open(sd_tsv):
    if l[0] != '>': continue
    q, t = l[1:].split('\t', 2)[:2]
    q = q.rsplit('_', 3); t = t.rsplit('_', 3)
    gq, gt = (q[0], int(q[1])), (t[0], int(t[1]))
    a, b = gene_locus.get(gq), gene_locus.get(gt)
    if a is None or b is None or a == b: continue
    if a > b: a, b, gq, gt = b, a, gt, gq
    h = hits[(a, b)]; h[0].add(gq); h[1].add(gt)

# candidate pairs from shared gene families (>= MIN_SHARED)
fid = {f: k for k, f in enumerate(sorted({f for fs in loc_fams for f in fs}))}
rows = np.repeat(np.arange(n_loc), [len(f) for f in loc_fams])
cols = np.array([fid[f] for fs in loc_fams for f in fs])
X = sp.csr_matrix((np.ones(len(cols), dtype=np.float32), (rows, cols)), shape=(n_loc, len(fid)))
I = sp.triu(X @ X.T, k=1).tocoo()
cand = set(zip(I.row[I.data >= MIN_SHARED].tolist(), I.col[I.data >= MIN_SHARED].tolist())) | set(hits)

src, dst, w = [], [], []
for a, b in cand:
    ha, hb = hits.get((a, b), (set(), set()))
    ma = {g for g in loc_genes[a] if fam_of[g] in loc_fams[b]} | ha
    mb = {g for g in loc_genes[b] if fam_of[g] in loc_fams[a]} | hb
    if min(len(ma), len(mb)) < MIN_SHARED: continue
    s = (len(ma) + len(mb)) / (len(loc_genes[a]) + len(loc_genes[b]))
    if s >= SMIN: src.append(a); dst.append(b); w.append(s)
n_all = len(w)
if BEST > 0:
    # best-match filter: keep a cross-bin edge only if it is (within BEST of) the best match of one of its
    # loci in the other bin; drop same-bin edges. Paralogous systems that share only a hub gene (ABC
    # transporters) then link to their ortholog, not to every paralog.
    lbin = L.bin.to_numpy(); top = defaultdict(float)
    for a, b, x in zip(src, dst, w):
        if lbin[a] != lbin[b]:
            top[(a, lbin[b])] = max(top[(a, lbin[b])], x); top[(b, lbin[a])] = max(top[(b, lbin[a])], x)
    kept = [(a, b, x) for a, b, x in zip(src, dst, w) if lbin[a] != lbin[b] and
            (x >= BEST * top[(a, lbin[b])] or x >= BEST * top[(b, lbin[a])])]
    src, dst, w = (list(t) for t in zip(*kept)) if kept else ([], [], [])
print(f'{n_loc} loci; {len(hits)} Spacedust-linked and {len(cand)} candidate pairs; {n_all} edges '
      f'(sim>={SMIN}, matched>={MIN_SHARED}), {len(w)} after best-match filter', file=sys.stderr)

g = ig.Graph(n=n_loc, edges=list(zip(src, dst)), edge_attrs={'weight': w})
cc = sorted((len(x) for x in g.connected_components()), reverse=True)
print(f'for comparison, connected components of the same graph: largest {cc[:3]} loci', file=sys.stderr)
part = la.find_partition(g, la.CPMVertexPartition, weights='weight', resolution_parameter=RES, seed=SEED, n_iterations=-1)
memb = np.array(part.membership)
sizes = Counter(memb.tolist())
ren = {c: k for k, c in enumerate(sorted(sizes, key=lambda c: (-sizes[c], c)))}
L['family'] = [f'CF{ren[c]:05d}' for c in memb]
sub_of = np.zeros(n_loc, dtype=int)
for c in sizes:
    vs = np.flatnonzero(memb == c)
    if len(vs) < 3: continue
    sg = g.induced_subgraph(vs.tolist())
    sp_ = la.find_partition(sg, la.CPMVertexPartition, weights='weight', resolution_parameter=RES_SUB, seed=SEED, n_iterations=-1)
    sm = np.array(sp_.membership); ss = Counter(sm.tolist())
    rr = {x: k for k, x in enumerate(sorted(ss, key=lambda x: (-ss[x], x)))}
    sub_of[vs] = [rr[x] for x in sm]
L['subfamily'] = [f'{f}.{k}' for f, k in zip(L.family, sub_of)]
L['genus'] = L.bin.map(lambda b: genus.get(b, b))
L['gene_families'] = [','.join(o.replace('(-)', '') for o in ordr) for ordr in loc_order]
L['order'] = [' '.join(o) for o in loc_order]
L['products'] = [' | '.join(prod_of[g] for g in gs) for gs in loc_genes]
wdeg = defaultdict(float)
for a, b, x in zip(src, dst, w):
    if memb[a] == memb[b]: wdeg[a] += x; wdeg[b] += x

cat = []
for fam, sub in L.groupby('family', sort=True):
    idx = sub.index.to_numpy(); n = len(idx)
    cnt = Counter(f for i in idx for f in loc_fams[i])
    core = [f for f, c in cnt.most_common() if c >= 0.5 * n] if n > 1 else sorted(loc_fams[idx[0]])
    med = max(idx, key=lambda i: (wdeg[i], -i))         # medoid: largest summed within-family similarity
    prods = []
    for f in core:                                     # most common Bakta product per core gene family
        pc = Counter(prod_of[g] for i in idx for g in loc_genes[i] if fam_of[g] == f)
        prods.append(pc.most_common(1)[0][0])
    mr = L.loc[med]
    cat.append(dict(family=fam, n_loci=n, n_bins=sub.bin.nunique(), n_genera=sub.genus.nunique(),
                    median_genes=float(sub.n_genes.median()), n_subfamilies=sub.subfamily.nunique(),
                    loci_per_bin=round(n / sub.bin.nunique(), 2), n_core=len(core),
                    medoid_locus=mr['locus'], typical_order=mr['order'], medoid_products=mr['products'],
                    core_families=','.join(core), core_products=' | '.join(prods), member_loci=','.join(sub.locus)))
C = pd.DataFrame(cat).sort_values(['n_loci', 'family'], ascending=[False, True])
C.to_csv(f'{out}_catalogue.tsv', sep='\t', index=False)
L.to_csv(f'{out}_members.tsv', sep='\t', index=False)
multi = C[C.n_loci > 1]
print(f'{len(C)} families ({len(multi)} with >=2 loci holding {multi.n_loci.sum()} loci; {int((C.n_genera >= 2).sum())} in >=2 genera); '
      f'largest {C.n_loci.max()} loci / {C.n_bins.max()} bins', file=sys.stderr)
