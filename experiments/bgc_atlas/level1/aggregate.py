"""Summarise the BGC-family level-1 test.  usage: aggregate.py LEVEL1_DIR > results.txt"""
import glob, sys
import numpy as np, pandas as pd
from scipy.stats import wilcoxon

d = sys.argv[1]
T = pd.concat([pd.read_csv(f, sep='\t') for f in glob.glob(f'{d}/res/*.targets.tsv')])
G = pd.concat([pd.read_csv(f, sep='\t') for f in glob.glob(f'{d}/res/*.genomes.tsv')])
S = pd.read_csv(f'{d}/species.tsv', sep='\t').set_index('id')
preds = [c[5:] for c in T.columns if c.startswith('bits_')]
pd.set_option('display.width', 220)
print(f'# BGC-family level-1 test: {T.id.nunique()} species, {len(T)} families (targets), '
      f'{G.genome.nunique()} genomes\n')
print('## 1. % of prevalence bits removed (median over species of the per-species mean)\n')
sp = T.groupby('id')[[f'bits_{m}' for m in preds]].mean()
pct = pd.DataFrame({m: 100 * (1 - sp[f'bits_{m}'] / sp['bits_prev']) for m in preds if m != 'prev'})
print(pct.median().round(1).to_string())
print('\nbest relatedness predictor per species:')
rel = ['ker_gc', 'tree_gc_c', 'ker_bgc', 'tree_bgc_c']
print(sp[[f'bits_{m}' for m in rel]].idxmin(axis=1).str[5:].value_counts().to_string())
dd = sp.bits_tree_gc_c - sp.bits_tree_bgc_c
print(f'\ngene-content tree vs BGC-content tree: gene content better in {(dd < 0).mean():.0%} of species, '
      f'median {dd.median():+.4f} bits (Wilcoxon p {wilcoxon(dd).pvalue:.2g})')
print('\n## 2. Detection models (tree on gene content): completeness only vs strict intact vs 1 - P(lost)\n')
for m in ('tree_gc_frag', 'tree_gc_seen'):
    x = sp[f'bits_{m}'] - sp.bits_tree_gc_c
    print(f'{m} - tree_gc_c: median {x.median():+.4f} bits, worse for {(x > 0).mean():.0%} of species')
print('\nby genome source (mean bits per genome):')
print(G.groupby('source')[[f'bits_{m}' for m in ('prev', 'tree_gc_c', 'tree_gc_seen', 'tree_gc_frag')]].mean().round(4)
      .assign(genomes=G.groupby('source').size()).to_string())
G['n50_class'] = pd.cut(G.n50, [0, 2e4, 1e5, 1e6, 1e9], labels=['<20 kb', '20-100 kb', '0.1-1 Mb', '>1 Mb'])
print(G.groupby('n50_class', observed=True)[[f'bits_{m}' for m in ('prev', 'tree_gc_c', 'tree_gc_seen', 'tree_gc_frag')]]
      .mean().round(4).assign(genomes=G.groupby('n50_class', observed=True).size()).to_string())
print('\n## 3. By family tier and by species\n')
T['tier'] = np.where(T.gcf.astype(str).str.contains('-F-'), 'fragmented', 'seed')
print(T.groupby('tier')[[f'bits_{m}' for m in ('prev', 'tree_gc_c', 'tree_bgc_c')]].mean().round(4)
      .assign(targets=T.groupby('tier').size()).to_string())
out = pct[['ker_gc', 'tree_gc_c', 'ker_bgc', 'tree_bgc_c']].round(1).join(S[['species', 'phylum', 'n_genomes', 'n_fams']])
print('\n' + out.sort_values('tree_gc_c').to_string())
