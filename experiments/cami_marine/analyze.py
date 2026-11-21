"""CAMI II marine: does metapanclust surprise predict true bin contamination beyond CheckM2?

usage: analyze.py DIR > results.txt
DIR holds bins.tsv (truth), sp/genome_species.tsv, out/<bin>.{summary,contigs}.tsv, checkm2_out/quality_report.tsv,
gsa_pooled_mapping_short.binning (contig -> genome).
"""
import glob, os, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_predict, GroupKFold

d = sys.argv[1]
pd.set_option('display.width', 220); pd.set_option('display.max_columns', 30)
T = pd.read_csv(f'{d}/bins.tsv', sep='\t')
sp = pd.read_csv(f'{d}/sp/genome_species.tsv', sep='\t').set_index('genome')
S = pd.concat([pd.read_csv(f, sep='\t') for f in glob.glob(f'{d}/out/*.summary.tsv')])
S['bin'] = S['bin'].str.replace('.fna', '', regex=False)
ck = pd.read_csv(f'{d}/checkm2_out/quality_report.tsv', sep='\t').rename(
    columns={'Name': 'bin', 'Completeness': 'ck_compl', 'Contamination': 'ck_contam'})
X = T.merge(S, on='bin').merge(ck[['bin', 'ck_compl', 'ck_contam']], on='bin', how='left')
X['species'] = X.genome.map(sp.species)
X['second_species'] = X.second_genome.map(sp.species)
X['contam_kind'] = np.where(X.true_contam < 0.01, 'clean',
                            np.where(X.second_species == X.species, 'same_species', 'other_species'))
X['novel_per_Mb'] = X.novel_genes / X.Mb
X['novel_excess_per_Mb'] = (X.novel_genes - X.novel_expected) / X.Mb
print(f'# CAMI II marine: {len(X)} bins scored (of {len(T)}; bins whose majority genome has a PanGBank species '
      f'pangenome), {X.binner.nunique()} binners, {X.genome.nunique()} genomes, {X.species.nunique()} species\n')
print(X.contam_kind.value_counts().to_string(), '\n')


def auc(y, s):
    ok = np.isfinite(s)
    return roc_auc_score(y[ok], s[ok]) if len(set(y[ok])) == 2 else np.nan


print('## 1. Detecting contaminated bins (AUC; bins >= 50% true completeness)\n')
Y = X[X.true_compl >= 0.5].copy()
preds = ['ck_contam', 'novel_bits', 'novel_excess_per_Mb', 'Zplus', 'rare_present', 'excess_presence']
# combined: CheckM2 + metapanclust, logistic, CV grouped by genome
rows = []
for thr in (0.02, 0.05, 0.10):
    y = (Y.true_contam >= thr).astype(int).values
    r = dict(threshold=thr, n=len(Y), positives=int(y.sum()))
    for p in preds: r[p] = auc(y, Y[p].values.astype(float))
    F = np.column_stack([Y.ck_contam.fillna(0), np.log1p(Y.novel_bits), Y.Zplus, np.log1p(Y.rare_present)])
    g = Y.genome.astype('category').cat.codes.values
    pr = cross_val_predict(LogisticRegression(C=1.0, max_iter=1000), F, y, cv=GroupKFold(5), groups=g,
                           method='predict_proba')[:, 1]
    r['combined_cv'] = roc_auc_score(y, pr)
    pr0 = cross_val_predict(LogisticRegression(C=1.0, max_iter=1000), F[:, :1], y, cv=GroupKFold(5), groups=g,
                            method='predict_proba')[:, 1]
    r['checkm2_only_cv'] = roc_auc_score(y, pr0)
    rows.append(r)
print(pd.DataFrame(rows).round(3).to_string(index=False))
print('\nby contamination kind (>= 5% contamination vs clean < 1%):')
for kind in ('same_species', 'other_species'):
    Z = Y[(Y.contam_kind == kind) & (Y.true_contam >= 0.05) | (Y.contam_kind == 'clean')]
    y = (Z.contam_kind == kind).astype(int).values
    print(f'  {kind}: n {len(Z)}, pos {y.sum()}; ' + ', '.join(f'{p} {auc(y, Z[p].values.astype(float)):.3f}'
                                                            for p in preds))
print('\nSpearman with true contamination:')
print('  ' + ', '.join(f'{p} {spearmanr(Y.true_contam, Y[p], nan_policy="omit").statistic:.3f}' for p in preds))

print('\n## 2. Completeness: estimate vs truth (mean absolute error, percentage points)\n')
for nm, col, sc in (('metapanclust', 'completeness_est', 100), ('CheckM2', 'ck_compl', 1)):
    e = (X[col] * sc - X.true_compl * 100)
    print(f'  {nm}: MAE {e.abs().mean():.2f}, bias {e.mean():+.2f}, Spearman {spearmanr(X[col], X.true_compl).statistic:.3f}')

print('\n## 3. Contig level: foreign contigs inside contaminated bins\n')
gt = pd.read_csv(f'{d}/gsa_pooled_mapping_short.binning', sep='\t', comment='@', header=None,
                 names=['contig', 'genome', 'taxid', 'len']).set_index('contig')
rows = []
for b in X[X.true_contam >= 0.05].bin:
    f = f'{d}/out/{b}.contigs.tsv'
    if not os.path.exists(f): continue
    c = pd.read_csv(f, sep='\t')
    maj = X.loc[X.bin == b, 'genome'].iloc[0]
    c['foreign'] = c.contig.map(gt.genome) != maj
    c['foreign_kind'] = np.where(~c.foreign, 'native', np.where(
        c.contig.map(gt.genome).map(sp.species) == sp.species.get(maj), 'same_species', 'other_species'))
    c['len'] = c.contig.map(gt.len); c['bin'] = b
    rows.append(c)
C = pd.concat(rows)
C['div'] = 100 - C.median_pident
C['unanchored'] = (C.persistent == 0).astype(int)
for kind in ('same_species', 'other_species'):
    Z = C[C.foreign_kind.isin(['native', kind])]
    y = (Z.foreign_kind == kind).astype(int).values
    print(f'  {kind}: {len(Z)} contigs, {y.sum()} foreign; AUC novel_frac {auc(y, Z.novel_frac.values):.3f}, '
          f'divergence {auc(y, Z["div"].fillna(Z["div"].max()).values):.3f}, presence_excess '
          f'{auc(y, Z.presence_excess.values):.3f}, unanchored {auc(y, Z.unanchored.values.astype(float)):.3f}')
    w = Z.len
    print(f'     bp share of foreign contigs with novel_frac > 0.5: {w[(y == 1) & (Z.novel_frac > .5)].sum() / max(w[y == 1].sum(), 1):.3f}; '
          f'of native: {w[(y == 0) & (Z.novel_frac > .5)].sum() / max(w[y == 0].sum(), 1):.3f}')

print('\n## 4. Per binner: median true contamination, CheckM2, novel bits\n')
print(X.groupby('binner').agg(n=('bin', 'size'), true_contam=('true_contam', 'median'), ck=('ck_contam', 'median'),
                               novel_bits=('novel_bits', 'median'), Zplus=('Zplus', 'median')).round(3).to_string())
X.to_csv(f'{d}/analysis_bins.tsv', sep='\t', index=False)
