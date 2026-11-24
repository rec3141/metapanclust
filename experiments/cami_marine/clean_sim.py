"""Cleaning by novelty on CAMI bins: remove contigs whose genes are mostly novel to the reference and report the
change in true contamination and completeness.  usage: clean_sim.py DIR ANALYSIS_BINS.tsv OUTSUB
(ANALYSIS_BINS.tsv from analyze.py; contigs without a called gene stay in the bin and are not counted.)"""
import sys
import numpy as np, pandas as pd

d, ab, sub = sys.argv[1:4]
X = pd.read_csv(ab, sep='\t')
gt = pd.read_csv(f'{d}/gsa_pooled_mapping_short.binning', sep='\t', comment='@', header=None,
                 names=['contig', 'genome', 'taxid', 'len'])
glen = gt.groupby('genome').len.sum(); g = gt.set_index('contig')
rules = {'novel>0.5': lambda c: c.novel_frac > 0.5,
         'novel>0.5 & unanchored': lambda c: (c.novel_frac > 0.5) & (c.persistent == 0),
         'novel>0.8': lambda c: c.novel_frac > 0.8}
rows = []
for b, maj, kind in zip(X.bin, X.genome, X.contam_kind):
    c = pd.read_csv(f'{d}/{sub}/{b}.contigs.tsv', sep='\t')
    c['len'] = c.contig.map(g.len); c['gen'] = c.contig.map(g.genome)
    bp0 = c.len.sum(); nat0 = c.len[c.gen == maj].sum()
    for rn, f in rules.items():
        k = c[~f(c)]
        bp = k.len.sum(); nat = k.len[k.gen == maj].sum()
        rows.append(dict(rule=rn, kind=kind, contam0=1 - nat0 / bp0, contam1=1 - nat / max(bp, 1),
                         dcompl=(nat - nat0) / glen[maj]))
R = pd.DataFrame(rows)
R['dcontam'] = R.contam1 - R.contam0
print(R.groupby(['rule', 'kind']).agg(n=('kind', 'size'), contam_before=('contam0', 'median'),
                                      contam_after=('contam1', 'median'), mean_contam_drop=('dcontam', 'mean'),
                                      mean_compl_loss=('dcompl', 'mean')).round(4).to_string())
