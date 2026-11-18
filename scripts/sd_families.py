"""Pairwise Spacedust clusters -> per-bin loci (overlapping gene-index ranges merged) -> cluster families
(connected components of loci linked by a pairwise cluster). First-pass catalogue, no splitting of hubs.
usage: sd_families.py CLUSTERS_TSV OUT_PREFIX"""
import sys
from collections import defaultdict, Counter
sd, out = sys.argv[1:3]
def parse(p):
    c, i, a, b = p.rsplit('_', 3); return c, int(i)
clusters = []                                   # list of (qbin, [(contig, idx)], tbin, [(contig, idx)])
cur = None
for l in open(sd):
    if l.startswith('#'):
        f = l[1:].rstrip('\n').split('\t'); cur = [f[1].replace('.gff3', ''), [], f[2].replace('.gff3', ''), []]; clusters.append(cur)
    elif l.startswith('>'):
        q, t = l[1:].split('\t')[:2]; cur[1].append(parse(q)); cur[3].append(parse(t))
# loci: per (bin, contig) merge index intervals of cluster sides
iv = defaultdict(list)                           # (bin, contig) -> [(lo, hi, side_id)]
sides = []
for k, (qb, qg, tb, tg) in enumerate(clusters):
    for b, g in ((qb, qg), (tb, tg)):
        byc = defaultdict(list)
        for c, i in g: byc[c].append(i)
        for c, idx in byc.items():
            sid = len(sides); sides.append(k); iv[(b, c)].append((min(idx), max(idx), sid))
side_locus = {}; loci = []
for (b, c), lst in iv.items():
    lst.sort(); lo = hi = None; cur_id = None
    for a, z, sid in lst:
        if lo is None or a > hi + 1:
            cur_id = len(loci); loci.append((b, c, a, z)); lo, hi = a, z
        else:
            hi = max(hi, z); loci[cur_id] = (b, c, lo, hi)
        side_locus[sid] = cur_id
# families: union loci linked by a cluster
parent = list(range(len(loci)))
def find(x):
    while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
    return x
by_cluster = defaultdict(set)
for sid, k in enumerate(sides): by_cluster[k].add(side_locus[sid])
for ls in by_cluster.values():
    ls = list(ls)
    for x in ls[1:]: parent[find(x)] = find(ls[0])
fam = defaultdict(list)
for i in range(len(loci)): fam[find(i)].append(i)
with open(f'{out}_loci.tsv', 'w') as o:
    o.write('locus\tfamily\tbin\tcontig\tfirst_gene\tlast_gene\n')
    for f, ls in fam.items():
        for i in ls:
            b, c, a, z = loci[i]; o.write(f'L{i}\tF{f}\t{b}\t{c}\t{a}\t{z}\n')
nb = Counter(len({loci[i][0] for i in ls}) for ls in fam.values())
glen = sorted(loci[i][3] - loci[i][2] + 1 for i in range(len(loci)))
print(f'{len(clusters)} pairwise clusters -> {len(loci)} loci -> {len(fam)} families')
print(f'locus length in genes: median {glen[len(glen) // 2]}, 90th pct {glen[int(len(glen) * 0.9)]}, max {glen[-1]}')
bins_per = sorted((len({loci[i][0] for i in ls}) for ls in fam.values()), reverse=True)
print('largest families span bins:', bins_per[:10], '| families in 2 bins:', nb[2], ', 3-9:', sum(v for k, v in nb.items() if 3 <= k <= 9), ', 10+:', sum(v for k, v in nb.items() if k >= 10))
big = max(fam.values(), key=len); print(f'largest family: {len(big)} loci across {len({loci[i][0] for i in big})} bins (hub merging to split later)')
