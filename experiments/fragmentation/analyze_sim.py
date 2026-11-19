"""(1)+(2) Summarise the fragmentation simulation and check the detection model against it."""
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
import common as C
from detection import p_intact, p_lost

d = pd.read_csv(f'{C.OUT}/sim_clusters.tsv.gz', sep='\t')
d['partial'] = 1 - d.intact - d.lost
d['n50b'] = pd.cut(d.N50, [0, 10e3, 20e3, 50e3, 100e3, 250e3, 1e7], labels=['<10k', '10-20k', '20-50k', '50-100k', '100-250k', '>250k'])
d['spanb'] = pd.cut(d.span, [0, 2e3, 5e3, 10e3, 20e3, 1e7], labels=['<2kb', '2-5kb', '5-10kb', '10-20kb', '>20kb'])
d['kb'] = pd.cut(d.k, [1, 2, 3, 5, 9, 1000], labels=['2', '3', '4-5', '6-9', '10+'])
pd.set_option('display.width', 200)

print('## clusters and bins')
print(d.groupby('set').agg(instances=('span', lambda x: len(x) // d.bin.nunique()), median_span=('span', 'median'), median_k=('k', 'median')))
print(d.drop_duplicates('bin').n50b.value_counts().sort_index().to_string())

print('\n## fragmentation only (completeness ignored): fraction intact, by N50 x span (both sets)')
print(d.pivot_table(index='n50b', columns='spanb', values='frag_intact', aggfunc='mean', observed=True).round(3))
print('\n## fragmentation only, by N50 x genes')
print(d.pivot_table(index='n50b', columns='kb', values='frag_intact', aggfunc='mean', observed=True).round(3))
print('\n## with completeness: intact / partial / lost / detected(>=2 genes on a contig), by N50')
print(d.groupby(['set', 'n50b'], observed=True)[['completeness', 'intact', 'partial', 'lost', 'detected', 'end_touch']].mean().round(3))
print('\n## operons vs sd6 by span (fragmentation only)')
print(d.pivot_table(index='spanb', columns='set', values='frag_intact', aggfunc='mean', observed=True).round(3))

# model check: predicted per bin x cluster
pred_i, pred_f, pred_l = np.zeros(len(d)), np.zeros(len(d)), np.zeros(len(d))
for b, ix in d.groupby('bin').indices.items():
    L = list(C.fasta_lengths(f'{C.BINS}/{b}.fna').values())
    c = d.completeness.iloc[ix[0]]
    s = d.span.values[ix]
    pred_f[ix] = p_intact(s, L, 1.0); pred_i[ix] = p_intact(s, L, c); pred_l[ix] = p_lost(s, L, c)
d['pred_frag'], d['pred_intact'], d['pred_lost'] = pred_f, pred_i, pred_l
g = d.groupby(['bin', 'spanb'], observed=True)[['frag_intact', 'pred_frag', 'intact', 'pred_intact', 'lost', 'pred_lost']].mean()
print('\n## model vs simulation (bin x span-class means, n=%d cells)' % len(g))
for a, p in (('frag_intact', 'pred_frag'), ('intact', 'pred_intact'), ('lost', 'pred_lost')):
    e = g[p] - g[a]
    print(f'{a:12s} mean abs err {e.abs().mean():.4f}  max {e.abs().max():.4f}  bias {e.mean():+.4f}  r={np.corrcoef(g[a], g[p])[0,1]:.4f}')
e = d.pred_intact - d.intact
print(f'per cluster x bin (50 reps each): intact MAE {e.abs().mean():.4f}')

# one-parameter rule of thumb: P_intact ~ c * exp(-s / (a * N50))
def loss(a):
    return ((d.completeness * np.exp(-d.span / (a * d.N50)) - d.intact) ** 2).mean()
a = minimize_scalar(loss, bounds=(0.05, 20), method='bounded').x
r = d.completeness * np.exp(-d.span / (a * d.N50))
print(f'\n## rule of thumb P_intact ~ c*exp(-s/({a:.2f}*N50)): MAE {np.abs(r - d.intact).mean():.4f} '
      f'(exact formula {np.abs(d.pred_intact - d.intact).mean():.4f})')
d.drop(columns=['n50b', 'spanb', 'kb']).to_csv(f'{C.OUT}/sim_with_pred.tsv.gz', sep='\t', index=False)
