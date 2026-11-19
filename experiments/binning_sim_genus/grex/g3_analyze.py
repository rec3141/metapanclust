"""Grex run summary: merge-rule performance per (focal species, contaminant species, strict family cut).

Per pair and cut:
  ani, af            skani ANI / alignment fraction of the contaminant's genomes vs the focal held-out
                     genomes (mean; empty below skani's ~80% range); rep_ani from the selection
  strict_frac        fraction of the contaminant's genes in a focal family at this cut (strain: same
                     for the focal species' own held-out genomes)
  auc_*              P(contaminant merge scores higher than a same-species strain merge), at equal
                     added fraction and base, averaged over fractions (5-40%) and bases (whole, 70%)
  rej0_*             rejected by the zero-excess rule (excess > 0)
  rejcal_*           rejected by a calibrated rule: added bits/Mb above the 95th percentile of strain
                     merges among the *calibration* half of the held-out genomes (per base), evaluated
                     on the *test* half only
  piece_*            single >= 20-kb contigs offered to the 70% base: AUC vs own/strain contigs, rejection
usage: g3_analyze.py BASE   -> BASE/results/pairs_grex.tsv, prints markdown summary
"""
import glob, os, re, sys
import numpy as np, pandas as pd

B = sys.argv[1]
R = f'{B}/results'


def auc(pos, neg):
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    return (pos[:, None] > neg[None, :]).mean() + 0.5 * (pos[:, None] == neg[None, :]).mean()


def md(t, index=True):
    t = t.reset_index() if index else t
    cols = [' / '.join(map(str, c)) if isinstance(c, tuple) else str(c) for c in t.columns]
    out = ['| ' + ' | '.join(cols) + ' |', '|' + '---|' * len(cols)]
    out += ['| ' + ' | '.join('' if pd.isna(v) else (f'{v:.3g}' if isinstance(v, float) else str(v)) for v in r) + ' |'
            for r in t.itertuples(index=False)]
    return '\n'.join(out)


plan = pd.read_csv(f'{B}/plan.tsv', sep='\t', header=None, names=['acc', 'species', 'genus', 'role', 'tag', 'rep_ani'])
plan['ncbi'] = plan.acc.str[3:]
rows = []
for f in sorted(glob.glob(f'{R}/sim_*_id*.tsv')):
    if f.endswith('_strictfrac.tsv'):
        continue
    m = re.match(r'sim_(.+)_id([0-9.]+)\.tsv', os.path.basename(f))
    tag, ident = m.group(1), float(m.group(2))
    d = pd.read_csv(f, sep='\t')
    sfr = pd.read_csv(f.replace('.tsv', '_strictfrac.tsv'), sep='\t')
    p = plan[plan.tag == tag]
    focal = p[p.role == 'held'].species.iloc[0]
    held = set(p[p.role == 'held'].ncbi)
    a = pd.read_csv(f'{B}/ani/{tag}.tsv', sep='\t', header=None, names=['q', 'r', 'ani', 'af_q', 'af_r']) \
        if os.path.exists(f'{B}/ani/{tag}.tsv') else pd.DataFrame(columns=['q', 'r', 'ani', 'af_q', 'af_r'])
    sp_of = dict(zip(p.ncbi, p.species))
    a = a[a.r.isin(held) & ~a.q.isin(held)].assign(sp=lambda x: x.q.map(sp_of))
    merges = d[(d.kind != 'clean') & ~d.kind.str.startswith('piece_')]
    pieces = d[d.kind.str.startswith('piece_')]
    s = merges[merges.kind == 'strain']
    same_pc = pieces[pieces.kind.isin(['piece_own', 'piece_strain'])]
    T = {b: np.quantile(s[(s.base == b) & (s.half == 'calib')].d_S_perMbAdded, 0.95)
         if (s[(s.base == b) & (s.half == 'calib')]).shape[0] else np.nan for b in ('full', '70')}
    st = s[s.half == 'test']
    strain_rejcal = np.mean([(st[st.base == b].d_S_perMbAdded > T[b]).mean() for b in T])
    for kind, c in merges[merges.kind != 'strain'].groupby('kind'):
        if kind == 'far':
            csp = p[p.role == 'far'].species.iloc[0]
        else:
            ep = kind.split(':', 1)[1]
            csp = p[(p.role == 'contam') & (p.species.str.split(' ', n=1).str[1] == ep)].species.iloc[0]
        ac = a[a.sp == csp]
        ct = c[c.half == 'test']
        pc = pieces[pieces.kind == f'piece_{kind}']
        rows.append(dict(
            tag=tag, focal=focal, contaminant=csp, ident=ident,
            cls='far' if kind == 'far' else 'congener',
            rep_ani=p[(p.species == csp)].rep_ani.astype(float).max() if kind != 'far' else np.nan,
            ani=ac.ani.mean() if len(ac) else np.nan, af=ac.af_q.mean() if len(ac) else np.nan,
            strict_frac=sfr[sfr.src == kind].strict_frac.mean(), strain_strict_frac=sfr[sfr.src == 'self'].strict_frac.mean(),
            auc_excess=np.nanmean([auc(c[(c.frac == fr) & (c.base == b)].excess, s[(s.frac == fr) & (s.base == b)].excess)
                                   for fr in c.frac.unique() for b in ('full', '70')]),
            auc_foreign=np.nanmean([auc(c[(c.frac == fr) & (c.base == b)].d_foreign, s[(s.frac == fr) & (s.base == b)].d_foreign)
                                    for fr in c.frac.unique() for b in ('full', '70')]),
            auc_dup_rel=np.nanmean([auc(c[(c.frac == fr) & (c.base == b)].d_dup_rel, s[(s.frac == fr) & (s.base == b)].d_dup_rel)
                                    for fr in c.frac.unique() for b in ('full', '70')]),
            rej0_full=(c[c.base == 'full'].excess > 0).mean(), rej0_70=(c[c.base == '70'].excess > 0).mean(),
            rej0_strain=(s.excess > 0).mean(),
            rejcal_full=(ct[ct.base == 'full'].d_S_perMbAdded > T['full']).mean(),
            rejcal_70=(ct[ct.base == '70'].d_S_perMbAdded > T['70']).mean(),
            rejcal_strain=strain_rejcal,
            piece_auc=auc(pc.excess, same_pc.excess), piece_rej0=(pc.excess > 0).mean(),
            piece_rej0_own=(pieces[pieces.kind == 'piece_own'].excess > 0).mean()))
P = pd.DataFrame(rows)
if P.empty:
    sys.exit('no simulation results found')
P.to_csv(f'{R}/pairs_grex.tsv', sep='\t', index=False)

bins = [0, 80, 85, 88, 90, 92, 94, 95, 96, 101]
lab = ['no skani (<~80)', '80-85', '85-88', '88-90', '90-92', '92-94', '94-95', '95-96', '>=96']
P['ani_bin'] = pd.cut(P.ani.fillna(0), bins, labels=lab, right=False)
P['sf_bin'] = pd.cut(P.strict_frac, [0, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 1.01], right=False)
cols = ['auc_excess', 'auc_foreign', 'auc_dup_rel', 'rej0_full', 'rej0_70', 'rej0_strain', 'rejcal_full', 'rejcal_70',
        'rejcal_strain', 'piece_auc', 'piece_rej0', 'strict_frac']
print(f'# Same-genus contaminants, grex run\n\n{P.tag.nunique()} focal species, {P.contaminant.nunique()} contaminant '
      f'species, {len(P[P.ident == P.ident.min()])} focal x contaminant pairs; strict cuts {[float(x) for x in sorted(P.ident.unique())]}.\n')
for ident, d in P.groupby('ident'):
    g = d.groupby('ani_bin', observed=True)
    t = g[cols].median().assign(n=g.size(), auc_excess_q10=g.auc_excess.quantile(0.1))
    print(f'\n## Strict cut {ident:.2f}: by genome ANI bin (medians; n = pairs)\n')
    print(md(t.round(3)))
for ident, d in P.groupby('ident'):
    g = d.groupby('sf_bin', observed=True)
    print(f'\n## Strict cut {ident:.2f}: by shared-gene fraction (contaminant genes in focal families)\n')
    print(md(g[['auc_excess', 'rej0_full', 'rej0_70', 'rejcal_full', 'rejcal_70', 'ani']].median().assign(n=g.size()).round(3)))
# where does it break: highest ANI bin whose median AUC >= 0.9, per cut
print('\n## Breakdown point (lowest ANI bin with median AUC < 0.9 / < 0.75)\n')
for ident, d in P[P.ani.notna()].groupby('ident'):
    m = d.groupby('ani_bin', observed=True).auc_excess.median()
    lo9 = next((str(k) for k, v in m.items() if v < 0.9), 'none')
    lo75 = next((str(k) for k, v in m.items() if v < 0.75), 'none')
    print(f'- cut {ident:.2f}: AUC < 0.9 from {lo9}; AUC < 0.75 from {lo75}')
