"""Compare MAGLean certificates on CAMI II marine bins with the gold standard and CheckM2.

usage: compare.py OUTDIR MAGLEAN_BIN BINS.tsv CHECKM2_REPORT.tsv > results.txt
Runs `maglean --json` on OUTDIR/evidence/*.json (cached in OUTDIR/certs/), then reports:
  1. completeness / contamination vs truth (MAE, bias, Spearman) for MAGLean and CheckM2;
  2. detecting truly contaminated bins (AUC) from MAGLean marker contamination, conservative taxonomic impurity
     (1 - purity_cons) per rank, max chimerism, and CheckM2 contamination;
  3. alert frequencies by true contamination class.
"""
import glob, json, os, subprocess, sys
from fractions import Fraction
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

out, mbin, bins_f, ck_f = sys.argv[1:5]
os.makedirs(f'{out}/certs', exist_ok=True)
fr = lambda s: float(Fraction(s)) if s not in (None, '') else np.nan
rows = []
for ev in sorted(glob.glob(f'{out}/evidence/*.json')):
    b = os.path.basename(ev)[:-5]
    cf = f'{out}/certs/{b}.json'
    if not os.path.exists(cf):
        r = subprocess.run([mbin, '--json', ev], capture_output=True, text=True)
        if r.returncode:
            print(f'# {b}: maglean exit {r.returncode}: {r.stderr.strip()[:200]}', file=sys.stderr); continue
        open(cf, 'w').write(r.stdout)
    c = json.load(open(cf))
    row = dict(bin=b, ml_compl=fr(c['completeness_exact']), ml_contam=fr(c['contamination_exact']),
               tier=c['mimag_tier'], tax_coherent=c['tax_coherent'], n_alerts=len(c['alerts']),
               stealth='STEALTH_CHIMERA' in ' '.join(c['alerts']), incoherent='TAXONOMIC_INCOHERENCE' in ' '.join(c['alerts']),
               dilution='UNCLASSIFIED_DILUTION' in ' '.join(c['alerts']), accessory_bp=c['accessory_bp'] / c['total_bp'])
    for ra in c['rank_audits']:
        row[f'imp_{ra["rank"]}'] = 1 - fr(ra['purity_cons'])
        row[f'chim_{ra["rank"]}'] = fr(ra['chimerism'])
        row[f'uncl_{ra["rank"]}'] = fr(ra['unclassified_frac'])
    rows.append(row)
M = pd.DataFrame(rows)
T = pd.read_csv(bins_f, sep='\t')
ck = pd.read_csv(ck_f, sep='\t').rename(columns={'Name': 'bin', 'Completeness': 'ck_compl', 'Contamination': 'ck_contam'})
X = M.merge(T, on='bin').merge(ck[['bin', 'ck_compl', 'ck_contam']], on='bin', how='left')
ranks = [c[4:] for c in X.columns if c.startswith('imp_')]
pd.set_option('display.width', 220)
print(f'# MAGLean on CAMI II marine: {len(X)} bins from {X.binner.nunique()} binners '
      f'(markers: CheckM domain-level sets; contig taxonomy: Kraken2 standard DB)\n')
print('## 1. Estimates vs truth (percentage points)\n')
for nm, cc, ct, sc in (('MAGLean', 'ml_compl', 'ml_contam', 100), ('CheckM2', 'ck_compl', 'ck_contam', 1)):
    ec = X[cc] * sc - X.true_compl * 100; et = X[ct] * sc - X.true_contam * 100
    print(f'{nm:8s} completeness MAE {ec.abs().mean():5.1f} bias {ec.mean():+5.1f} rho {spearmanr(X[cc], X.true_compl).statistic:.3f}; '
          f'contamination MAE {et.abs().mean():5.1f} bias {et.mean():+5.1f} rho {spearmanr(X[ct], X.true_contam).statistic:.3f}')
print('\n## 2. Detecting contaminated bins (AUC)\n')
X['max_chim'] = X[[f'chim_{r}' for r in ranks]].max(axis=1)
scores = ['ck_contam', 'ml_contam'] + [f'imp_{r}' for r in ranks] + ['max_chim']
res = []
for thr in (0.05, 0.10):
    y = (X.true_contam >= thr).astype(int)
    r = dict(true_contam_ge=thr, positives=int(y.sum()), n=len(y))
    for s in scores:
        v = X[s].fillna(0).values
        r[s] = round(roc_auc_score(y, v), 3) if 0 < y.sum() < len(y) else np.nan
    res.append(r)
print(pd.DataFrame(res).to_string(index=False))
print('\n## 3. Alerts by true contamination\n')
X['true_class'] = pd.cut(X.true_contam, [-1, .01, .05, .10, 1.1], labels=['<1%', '1-5%', '5-10%', '>=10%'])
print(X.groupby('true_class', observed=True).agg(n=('bin', 'size'), stealth=('stealth', 'mean'),
      incoherent=('incoherent', 'mean'), dilution=('dilution', 'mean'), ml_contam=('ml_contam', 'median'),
      ck_contam=('ck_contam', 'median')).round(3).to_string())
print('\nMIMAG tiers (marker-only; RNA unknown so HQ is impossible):')
print(pd.crosstab(X.tier, X.true_class).to_string())
X.to_csv(f'{out}/maglean_vs_truth.tsv', sep='\t', index=False)
