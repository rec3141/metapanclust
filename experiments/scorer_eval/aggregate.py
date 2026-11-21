"""Summarise the metapanclust v0 held-out evaluation.  usage: aggregate.py OUT_DIR > results.txt

OUT_DIR/<genus>/shard_<i>/{bins.tsv, contigs.tsv.gz, units_<genome>_<comp>.tsv.gz}
"""
import glob, os, sys
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score

out = sys.argv[1]
pd.set_option('display.width', 220); pd.set_option('display.max_columns', 30)


def auc(pos, neg):
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0: return np.nan
    return roc_auc_score(np.r_[np.ones(len(pos)), np.zeros(len(neg))], np.r_[pos, neg])


def bt(q, y):
    q = np.clip(q, 1e-9, 1 - 1e-9)
    return -(y * np.log2(q) + (1 - y) * np.log2(1 - q))


B, C, U = [], [], []
for g in sorted(os.listdir(out)):
    for sh in glob.glob(f'{out}/{g}/shard_*'):
        if not os.path.exists(f'{sh}/bins.tsv'): continue
        B.append(pd.read_csv(f'{sh}/bins.tsv', sep='\t').assign(genus=g))
        C.append(pd.read_csv(f'{sh}/contigs.tsv.gz', sep='\t').assign(genus=g))
        for f in glob.glob(f'{sh}/units_*.tsv.gz'):
            U.append(pd.read_csv(f, sep='\t').assign(genus=g))
B = pd.concat(B); C = pd.concat(C); U = pd.concat(U)
print(f'# metapanclust v0 held-out evaluation: {B.genus.nunique()} species, {B.held.nunique()} held-out genomes '
      f'(per species up to 20), {len(B)} simulated bins\n')

print('## 1. Unit prediction on clean held-out bins (bits per unit; lower is better)\n')
U['y'] = U.observed.astype(int)
rows = []
for (g, c), d in U.groupby(['genus', 'target_comp']):
    N = d.ref_carriers.max() + 1
    pp = (d.ref_carriers + 1) / (N + 2)
    rows.append(dict(genus=g, comp=c, units=len(d) // d.held.nunique(), model=d.bits.mean(),
                     tree=bt(d.p_tree * d.detect + (1 - d.p_tree) * 0.002, d.y).mean(),
                     prevalence=bt(pp * d.detect + (1 - pp) * 0.002, d.y).mean()))
R = pd.DataFrame(rows)
R['gain_vs_prev_%'] = 100 * (1 - R.model / R.prevalence)
R['context_gain_%'] = 100 * (1 - R.model / R.tree)
print(R.round(4).to_string(index=False))
print('\nmedian gain vs prevalence: %.1f%%, context over tree: %.1f%%' % (R['gain_vs_prev_%'].median(),
                                                                         R['context_gain_%'].median()))
U['qb'] = pd.cut(U.q, [0, .001, .01, .05, .2, .5, .8, .95, 1])
print('\ncalibration (all species, clean bins):')
print(U.groupby('qb', observed=True).agg(n=('y', 'size'), mean_q=('q', 'mean'), observed=('y', 'mean'),
                                         excess=('excess', 'mean')).round(4).to_string())
cl = B[B.cond == 'clean']
print(f'\nZ+ on clean bins: median {cl.Zplus.median():.2f} (IQR {cl.Zplus.quantile(.25):.2f} .. '
      f'{cl.Zplus.quantile(.75):.2f}); Z-: median {cl.Zminus.median():.2f}  (0 if calibrated per genome)')

print('\n## 2. Bin level: does the score separate contaminated or strain-merged bins from clean ones?\n')
scores = ['novel_bits', 'Zplus', 'rare_present', 'excess_presence']
rows = []
for c in sorted(B.target_comp.unique()):
    b = B[B.target_comp == c]
    clean = b[b.cond == 'clean']
    for cond, x in [('foreign', 0.05), ('foreign', 0.1), ('foreign', 0.2), ('strain', 0.1)]:
        t = b[(b.cond == cond) & (b.contam == x)]
        r = dict(comp=c, cond=cond, contam=x, n=len(t))
        for s in scores: r[f'AUC_{s}'] = auc(t[s], clean[s])
        rows.append(r)
print(pd.DataFrame(rows).round(3).to_string(index=False))
print('\nforeign vs strain merges (10%):')
for c in sorted(B.target_comp.unique()):
    b = B[B.target_comp == c]
    f = b[(b.cond == 'foreign') & (b.contam == 0.1)]; s = b[b.cond == 'strain']
    print(f'  comp {c}: novel_bits AUC {auc(f.novel_bits, s.novel_bits):.3f}; rare_present (strain > foreign) '
          f'AUC {auc(s.rare_present, f.rare_present):.3f}')
print('\nper species, foreign 5% vs clean, novel_bits AUC:')
print(B.groupby('genus').apply(lambda d: auc(d[(d.cond == 'foreign') & (d.contam == 0.05)].novel_bits,
                                             d[d.cond == 'clean'].novel_bits), include_groups=False).round(3).to_string())
print('\ncompleteness estimate (clean bins): target vs estimate')
print(cl.groupby('target_comp').completeness.describe()[['mean', '25%', '75%']].round(3).to_string())

print('\n## 3. Contig level: finding the foreign contigs inside contaminated bins\n')
Cf = C[C.cond == 'foreign'].copy()
Cf['is_foreign'] = Cf.src == 'foreign'
Cf['score_novel'] = Cf.novel_frac
Cf['score_excess'] = Cf.presence_excess
Cf['anchored'] = Cf.persistent > 0
rows = []
for (c, x), d in Cf.groupby(['target_comp', 'contam']):
    r = dict(comp=c, contam=x, contigs=len(d), foreign=int(d.is_foreign.sum()))
    for s in ('score_novel', 'score_excess'):     # identity is not simulated (native genes keep their labels)
        r[f'AUC_{s[6:]}'] = auc(d[d.is_foreign][s], d[~d.is_foreign][s])
    w = d.len.fillna(0)
    r['foreign_bp_unanchored'] = round(w[d.is_foreign & ~d.anchored].sum() / max(w[d.is_foreign].sum(), 1), 3)
    r['native_bp_unanchored'] = round(w[~d.is_foreign & ~d.anchored].sum() / max(w[~d.is_foreign].sum(), 1), 3)
    rows.append(r)
print(pd.DataFrame(rows).round(3).to_string(index=False))
