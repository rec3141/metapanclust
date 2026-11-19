"""Share of each focal-bin contig covered (>= 98% identity) by contigs of a bin from another co-assembly.
usage: marine_check.py FOCAL_CONTIGS_TSV MMSEQS_M8 OUT_TSV   (m8 from mmseqs easy-search --search-type 3)"""
import sys
import numpy as np, pandas as pd
o = pd.read_csv(sys.argv[1], sep='\t', index_col=0)
h = pd.read_csv(sys.argv[2], sep='\t', header=None, names='q t pid aln mm gap qs qe ts te ev bits'.split())
h = h[h.pid >= (0.98 if h.pid.max() <= 1 else 98)]
cov = {}
for q, d in h.groupby('q'):
    m = np.zeros(int(o.len[q]), bool)
    for a, b in zip(d.qs, d.qe): m[min(a, b) - 1:max(a, b)] = True
    cov[q] = m.mean()
o['other_cov98'] = pd.Series(cov).reindex(o.index).fillna(0)
o['grp'] = np.where(o.set == 'module', 'module:' + o['class'].str.split('+').str[0], o.set)
print(o.groupby('grp').other_cov98.describe()[['count', '25%', '50%', '75%']].round(2))
o[['set', 'class', 'other_cov98']].to_csv(sys.argv[3], sep='\t')
