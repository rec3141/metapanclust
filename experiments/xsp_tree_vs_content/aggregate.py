"""Aggregate the between-species tree-vs-content run.  usage: aggregate.py OUTDIR/res > results.txt"""
import glob, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr, wilcoxon

res = sys.argv[1]
SP = pd.concat([pd.read_csv(f, sep='\t') for f in sorted(glob.glob(f'{res}/*.species.tsv'))])
DG = pd.concat([pd.read_csv(f, sep='\t') for f in sorted(glob.glob(f'{res}/*.diag.tsv'))])
MC = pd.concat([pd.read_csv(f, sep='\t') for f in sorted(glob.glob(f'{res}/*.matcor.tsv'))])
X = SP.merge(DG.drop(columns=['species', 'genus', 'gc']), on=['clade', 'id'])
preds = [c[5:] for c in SP.columns if c.startswith('bits_')]
pd.set_option('display.width', 220); pd.set_option('display.max_columns', 40)
print(f'# between-species tree vs gene content: {SP.clade.nunique()} clades, {len(SP)} species\n')

print('## 1. % of prevalence loss removed (leave-one-species-out; mean over species within clade)\n')
cm = SP.groupby('clade')[[f'bits_{m}' for m in preds]].mean()
pct = pd.DataFrame({m: 100 * (1 - cm[f'bits_{m}'] / cm['bits_prev']) for m in preds if m != 'prev'})
pct.insert(0, 'species', SP.groupby('clade').size())
pct.insert(1, 'prev_bits', cm['bits_prev'].round(3))
print(pct.round(1).to_string())
print('\nmedian over clades:'); print(pct.drop(columns=['species', 'prev_bits']).median().round(1).to_string())

print('\n## 2. Best relatedness predictor (lowest bits among ker/tree x gtdb/aai/jac)\n')
rel = [m for m in preds if m.split('_')[0] in ('ker', 'tree')]
X['best_rel'] = X[[f'bits_{m}' for m in rel]].idxmin(1).str[5:]
print('species:'); print(X.best_rel.value_counts().to_string())
cb = pct[rel].idxmax(1)
print('\nclades:'); print(cb.value_counts().to_string())
for a_, b_ in (('tree_gtdb', 'tree_jac'), ('tree_gtdb', 'tree_aai'), ('tree_aai', 'tree_jac')):
    d = X[f'bits_{a_}'] - X[f'bits_{b_}']
    print(f'{a_} - {b_}: median {d.median():+.4f} bits/target, {a_} better for {(d < 0).mean():.0%} of species, '
          f'Wilcoxon p {wilcoxon(d).pvalue:.2g}')

print('\n## 3. Context beyond relatedness (% of prevalence loss, median over clades)\n')
for r in ('relseq', 'relall'):
    print(f'{r}: rel {pct[r].median():.1f}, ctx (lr) {pct["lr"].median():.1f}, avg {pct[f"avg_{r}"].median():.1f}; '
          f'avg - max(lr, rel) median {(pct[f"avg_{r}"] - pct[["lr", r]].max(axis=1)).median():+.1f}')
d = X.bits_lr - X.bits_relall
print(f'lr better than relall for {(d < 0).mean():.0%} of species')

print('\n## 4. Distance correlations (Spearman over species pairs; median and range over clades)\n')
MC['pair'] = MC.a + ' ~ ' + MC.b
print(MC.groupby('pair').spearman.describe()[['count', 'min', '50%', 'max']].round(3).to_string())

print('\n## 5. GC: is the GTDB tree or gene content better for GC-shifted species?\n')
X['abs_gcdev'] = X.gc_dev.abs()
X['gtdb_minus_jac'] = X.bits_tree_gtdb - X.bits_tree_jac          # > 0: gene content predicts better
X['gtdb_minus_aai'] = X.bits_tree_gtdb - X.bits_tree_aai
rho = spearmanr(X.abs_gcdev, X.gtdb_minus_jac)
print(f'pooled: Spearman(|GC dev|, bits_tree_gtdb - bits_tree_jac) = {rho.statistic:.3f} (p {rho.pvalue:.2g}, '
      f'n {len(X)}); positive = gene content gains on GC-shifted species')
pc = X.groupby('clade').apply(lambda g: spearmanr(g.abs_gcdev, g.gtdb_minus_jac).statistic, include_groups=False)
print(f'per clade: median {pc.median():.3f}, {(pc > 0).sum()}/{len(pc)} positive')
X['gc_class'] = pd.cut(X.abs_gcdev, [-1, 2, 5, 10, 100], labels=['<2', '2-5', '5-10', '>10'])
print(X.groupby('gc_class', observed=True).agg(n=('id', 'size'), tree_gtdb=('bits_tree_gtdb', 'mean'),
      tree_aai=('bits_tree_aai', 'mean'), tree_jac=('bits_tree_jac', 'mean'),
      jac_better=('gtdb_minus_jac', lambda v: (v > 0).mean()), r2t_z=('r2t_gtdb_z', 'mean')).round(3).to_string())
rho = spearmanr(X.abs_gcdev, X.r2t_gtdb_z)
print(f'|GC dev| vs GTDB root-to-tip z: Spearman {rho.statistic:.3f} (p {rho.pvalue:.2g}) (long branches from '
      f'composition?)')

print('\n## 6. Discordant species\n')
# tree outlier, content typical ("MGI-like"): far from everything on the GTDB tree, close to something by content
X['tree_out_content_in'] = ((X.nnpct_gtdb >= 0.8) | (X.r2t_gtdb_z > 1.5)) & (X.nnpct_jac <= 0.5)
X['content_out_tree_in'] = (X.nnpct_jac >= 0.8) & (X.nnpct_gtdb <= 0.5)
X['nn_differ'] = X.nn_gtdb != X.nn_jac
for f in ('tree_out_content_in', 'content_out_tree_in', 'nn_differ'):
    g = X[X[f]]; o = X[~X[f]]
    print(f'{f}: {len(g)} species; |GC dev| {g.abs_gcdev.median():.2f} vs {o.abs_gcdev.median():.2f}; '
          f'tree_gtdb {g.bits_tree_gtdb.mean():.3f} tree_aai {g.bits_tree_aai.mean():.3f} '
          f'tree_jac {g.bits_tree_jac.mean():.3f} lr {g.bits_lr.mean():.3f}; '
          f'gene content better for {(g.gtdb_minus_jac > 0).mean():.0%}')
nd = X[X.nn_differ]
print(f'\nnearest neighbour differs (gtdb vs jac): {len(nd)} species; |GC diff| to GTDB nn '
      f'{nd.nn_gcdiff_gtdb.median():.2f} vs to content nn {nd.nn_gcdiff_jac.median():.2f} '
      f'(content nn closer in GC for {(nd.nn_gcdiff_jac < nd.nn_gcdiff_gtdb).mean():.0%})')
cols = ['clade', 'species', 'gc', 'gc_dev', 'r2t_gtdb_z', 'nn_gtdb', 'nn_gcdiff_gtdb', 'nn_jac', 'nn_gcdiff_jac',
        'bits_tree_gtdb', 'bits_tree_jac', 'bits_lr']
print('\ntree outlier, content typical:')
print(X[X.tree_out_content_in].sort_values('abs_gcdev', ascending=False)[cols].head(25).to_string(index=False))
print('\nlargest |GC dev| per clade:')
print(X.sort_values('abs_gcdev', ascending=False).drop_duplicates('clade')[cols].to_string(index=False))
print('\nNitrososphaeria (Marine Group I = Nitrosopumilaceae):')
print(X[X.clade == 'Nitrososphaeria'].sort_values('gc')[cols].to_string(index=False))
