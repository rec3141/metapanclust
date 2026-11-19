"""(3b) Are assembly breaks placed at random relative to conserved clusters? Fraction of CDS that sit in
any Spacedust cluster instance, and fraction annotated 'hypothetical protein', by distance (in genes)
from the nearest contig end, Chesterfield bins. Contigs with >= 40 CDS only, so the middle is far from
both ends."""
import glob, os
import numpy as np, pandas as pd
import common as C

member = set()
for c, mem in C.sd_instances(f'{C.D}/chesterfield/sd_bins_mmseqs/clusters.tsv'):
    for s, e in mem:
        member.add((c, s, e))
rows = []
for g in glob.glob(f'{C.BINS}/*.gff3'):
    by = {}
    for l in open(g):
        f = l.split('\t')
        if len(f) < 9 or f[2] != 'CDS':
            continue
        by.setdefault(f[0], []).append((int(f[3]), int(f[4]), 'hypothetical protein' in f[8]))
    for c, x in by.items():
        x.sort(); n = len(x)
        if n < 40:
            continue
        for i, (s, e, h) in enumerate(x):
            rows.append((min(i, n - 1 - i), (c, s, e) in member, h))
r = pd.DataFrame(rows, columns=['d', 'in_cluster', 'hypothetical'])
r['dist'] = np.minimum(r.d, 20)
t = r.groupby('dist').agg(cds=('d', 'size'), in_cluster=('in_cluster', 'mean'), hypothetical=('hypothetical', 'mean')).round(3)
print(t.to_string())
