"""Calibrated alternative to the zero-excess merge rule: reject a merge whose added contigs cost more
bits per Mb than the 95th percentile of same-species (strain) merges into the same focal species.
The threshold is set on the held-out strain merges themselves, so this is an optimistic bound on
what a per-species calibration (e.g. from training genomes) could achieve.
usage: 07_threshold.py RESULTS_DIR   (reads sim_*.tsv and pairs.tsv; prints a markdown table)
"""
import glob, sys
import numpy as np, pandas as pd

R = sys.argv[1]
d = pd.concat([pd.read_csv(f, sep='\t') for f in glob.glob(f'{R}/sim_*.tsv') if 'strictfrac' not in f])
d = d[(d.kind != 'clean') & ~d.kind.str.startswith('piece_')]
P = pd.read_csv(f'{R}/pairs.tsv', sep='\t')
rows = []
for (A, base), g in d.groupby(['A', 'base']):
    s = g[g.kind == 'strain']
    T = np.quantile(s.d_S_perMbAdded, 0.95)
    for kind, c in g[g.kind != 'strain'].groupby('kind'):
        rows.append(dict(A=A, base=base, kind=kind, T=T, rej_cal=(c.d_S_perMbAdded > T).mean(),
                         rej_zero=(c.excess > 0).mean()))
X = pd.DataFrame(rows)
X['contaminant'] = X.kind.str.replace('genus:', '', regex=False)
P['A'] = P.focal.str.split(' ').str[1].str.lower()
P['contaminant'] = P.contaminant.str.split(' ').str[1]
P.loc[P.contaminant == 'psychrophilum', 'contaminant'] = 'far'
X = X.merge(P[['A', 'contaminant', 'ani', 'strict_frac']], on=['A', 'contaminant'], how='left')
X['ani_bin'] = pd.cut(X.ani.fillna(0), [0, 80, 85, 90, 93, 95, 101],
                      labels=['<80 / no skani', '80-85', '85-90', '90-93', '93-95', '>=95'], right=False)
t = X.pivot_table(index='ani_bin', columns='base', values=['rej_zero', 'rej_cal'], aggfunc='mean', observed=True).round(3)
print('| ANI bin | ' + ' | '.join(f'{a} ({b} base)' for a, b in t.columns) + ' |')
print('|' + '---|' * (len(t.columns) + 1))
for i, r in t.iterrows():
    print(f'| {i} | ' + ' | '.join(str(v) for v in r.values) + ' |')
X.round(3).to_csv(f'{R}/calibrated_threshold.tsv', sep='\t', index=False)
