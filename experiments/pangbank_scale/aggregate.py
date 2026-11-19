"""Aggregate worker outputs of the PanGBank data-size experiment.

usage: aggregate.py OUTDIR CATALOG.tsv
Reads OUTDIR/A/*.summary.tsv and OUTDIR/B/*.summary.tsv (+ *.kappa.tsv); writes OUTDIR/agg/*.tsv and prints
markdown tables.
  A: gain of the context model over Beta(1,1) prevalence vs subsample size N, in bits per (genome, target)
     and as % of the prevalence loss; C is fixed per row (no per-species tuning); slope per doubling of N;
     by phylum.
  B: held-out bits per genome for species-only (Beta(1,1) and EB Beta(a,b)), genus mean and genus-shrunk
     prevalence, per training size m and prevalence stratum; fraction of focal species where shrinkage wins.
"""
import glob, os, sys
import numpy as np, pandas as pd

out, catf = sys.argv[1], sys.argv[2]
os.makedirs(f'{out}/agg', exist_ok=True)
cat = pd.read_csv(catf, sep='\t').rename(columns={'id': 'pid'})
pd.set_option('display.width', 250)


def md(df, floatfmt='.3f'):
    try: return df.to_markdown(floatfmt=floatfmt)
    except ImportError: return df.to_string()


# ------------------------------------------------------------------------------------------------- A
fa = sorted(glob.glob(f'{out}/A/*.summary.tsv'))
if fa:
    A = pd.concat([pd.read_csv(f, sep='\t') for f in fa], ignore_index=True)
    A = A[A.n_targets > 0].merge(cat[['pid', 'phylum', 'genus', 'species', 'mean_compl']], left_on='id', right_on='pid')
    A['Nlab'] = np.where(A.N == A.n_total, 'all', A.N.astype(str))
    lr = [c for c in A.columns if c.startswith('bits_lr')]
    models = ['bits_prev'] + [c for c in A.columns if c.startswith('bits_knn')] + lr + ['bits_mix0.1']
    # mean over reps, per species x N x kind (N = 'all' rows kept separately)
    P = A.groupby(['id', 'phylum', 'n_total', 'N', 'Nlab', 'kind'])[models + ['n_targets', 'sec']].mean().reset_index()
    for m in models[1:]:
        P['gain_' + m[5:]] = P.bits_prev - P[m]
        P['pct_' + m[5:]] = 100 * (P.bits_prev - P[m]) / P.bits_prev
    P.to_csv(f'{out}/agg/A_per_species.tsv', sep='\t', index=False)
    sub = P[P.Nlab != 'all']
    pcts = [c for c in P.columns if c.startswith('pct_')]
    T = sub.groupby(['kind', 'N']).agg(species=('id', 'nunique'), bits_prev=('bits_prev', 'mean'),
                                       **{c: (c, 'median') for c in pcts}).round(2)
    T.to_csv(f'{out}/agg/A_by_N.tsv', sep='\t')
    print(f'## A: {P.id.nunique()} species, {P.phylum.nunique()} phyla\n\nMedian % of prevalence loss removed, '
          f'by subsample size N (subsamples only):\n')
    print(md(T, '.1f'))
    # best fixed C per N (by median gain across species)
    best = sub.groupby(['kind', 'N'])[[c for c in pcts if c.startswith('pct_lr')]].median().idxmax(axis=1)
    print('\nBest fixed C per (kind, N):', best.to_dict())
    # gain at the best global C (0.1 default) per doubling: per-species regression of pct on log2 N
    sl = []
    for (sp, kind), g in sub.groupby(['id', 'kind']):
        if g.N.nunique() >= 3:
            for c in ('pct_lr0.1', 'pct_mix0.1', 'pct_knn2'):
                if c in g: sl.append(dict(id=sp, kind=kind, model=c, slope=np.polyfit(np.log2(g.N), g[c], 1)[0],
                                          top=g.N.max()))
    if sl:
        SL = pd.DataFrame(sl)
        print('\nPer-species slope of % gain per doubling of N (median, IQR):\n')
        print(md(SL.groupby(['kind', 'model']).slope.describe()[['count', '25%', '50%', '75%']], '.2f'))
        # still rising at the top? slope over the last two sizes for species reaching 400
        last = []
        for (sp, kind), g in sub.groupby(['id', 'kind']):
            g = g.sort_values('N')
            if len(g) >= 2 and g.N.iloc[-1] >= 200:
                last.append(dict(id=sp, kind=kind, N_hi=g.N.iloc[-1],
                                 d=g['pct_lr0.1'].iloc[-1] - g['pct_lr0.1'].iloc[-2]))
        if last:
            L = pd.DataFrame(last)
            print('\nChange in % gain (C = 0.1) between the two largest subsamples (N_hi = 200 or 400):\n')
            print(md(L.groupby(['kind', 'N_hi']).d.describe()[['count', '25%', '50%', '75%']], '.2f'))
    ph = sub[sub.N == 50].groupby(['kind', 'phylum']).agg(species=('id', 'nunique'),
                                                          pct_lr=('pct_lr0.1', 'median'),
                                                          pct_knn=('pct_knn2', 'median')).round(1)
    ph.to_csv(f'{out}/agg/A_by_phylum_N50.tsv', sep='\t')
    rt = A.groupby('N').sec.describe()[['count', '50%', 'max']]
    print('\nWall seconds per (N, rep) subsample, 1 thread:\n'); print(md(rt, '.0f'))

# ------------------------------------------------------------------------------------------------- B
fb = sorted(glob.glob(f'{out}/B/*.summary.tsv'))
if fb:
    B = pd.concat([pd.read_csv(f, sep='\t') for f in fb], ignore_index=True)
    mods = [c for c in ('bits_sp', 'bits_sp_eb', 'bits_genus', 'bits_shrunk', 'bits_mix', 'bits_oracle') if c in B]
    Q = B.groupby(['genus', 'focal', 'n_focal', 'n_species', 'm', 'stratum', 'kappa'])[mods + ['n_clusters']].mean().reset_index()
    Q['gain_vs_eb'] = Q.bits_sp_eb - Q.bits_shrunk
    Q['pct_vs_eb'] = 100 * Q.gain_vs_eb / Q.bits_sp_eb
    Q.to_csv(f'{out}/agg/B_per_focal.tsv', sep='\t', index=False)
    T = Q.groupby(['stratum', 'm']).agg(focal=('focal', 'nunique'), clusters=('n_clusters', 'mean'),
                                         **{m: (m, 'mean') for m in mods},
                                         pct_vs_eb_median=('pct_vs_eb', 'median'),
                                         frac_shrunk_wins=('gain_vs_eb', lambda v: (v > 0).mean())).round(2)
    T.to_csv(f'{out}/agg/B_by_m.tsv', sep='\t')
    print(f'\n## B: {Q.genus.nunique()} genera, {Q.focal.nunique()} focal species\n\nHeld-out bits per genome '
          f'(summed over genus clusters in the stratum), mean over focal species:\n')
    print(md(T, '.2f'))
    K = pd.concat([pd.read_csv(f, sep='\t') for f in sorted(glob.glob(f'{out}/B/*.kappa.tsv'))])
    print('\nFitted kappa (per focal species, from its congeners):', K.kappa.describe()[['count', '25%', '50%', '75%']].to_dict())
