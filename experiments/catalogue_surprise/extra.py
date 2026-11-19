"""Extra checks: calibration of predicted detection, system-level presence of atp/nuo/rp where the
main family is absent, partial correlations, duplicate copies in flagged bins."""
import numpy as np, pandas as pd, numpy.linalg as la
from scipy.stats import spearmanr, rankdata

D = '/data/scratch/metapanclust_data'; OUT = f'{D}/catalogue_surprise'
z = np.load(f'{OUT}/inputs.npz'); fz = np.load(f'{OUT}/fit_flat.npz')
obs, fams, bins = z['obs'], z['fams'], z['bins']
per = pd.read_csv(f'{OUT}/per_bin.tsv', sep='\t')
print(per[per.bin.isin(['dastool-maxbin_102', 'dastool-maxbin_040', 'dastool-vamb_711_sub'])][['bin', 'compl', 'contam', 'dup_intact_expected', 'S_pres_per_Mb', 'S_abs_per_expMb']].to_string())
for b in ['dastool-maxbin_102', 'dastool-maxbin_040']:
    v = per.set_index('bin').loc[b, 'dup_intact_expected']
    print(b, 'dup percentile', round((per.dup_intact_expected < v).mean(), 2))
X = np.c_[np.ones(len(per)), rankdata(per.compl)]
res = lambda y: y - X @ la.lstsq(X, y, rcond=None)[0]
print('partial rho absence~contam | compl', round(np.corrcoef(res(rankdata(per.S_abs_per_expMb)), res(rankdata(per.contam)))[0, 1], 2),
      '; rho compl~contam', round(spearmanr(per.compl, per.contam)[0], 2))
# calibration: predicted P(seen) vs observed, by model
seen = obs > 0
rows = []
for name, pk in [('tax_shrunk', 'p'), ('flat', 'p_flat')]:
    p = fz[pk]
    q = p * (fz['PI'] + fz['PP']) + (1 - p) * fz['eps']
    for lo, hi in [(0, .01), (.01, .05), (.05, .2), (.2, .5), (.5, .8), (.8, .95), (.95, 1.01)]:
        m = (q >= lo) & (q < hi)
        rows.append({'model': name, 'pred_bin': f'{lo}-{hi}', 'n': int(m.sum()), 'pred_mean': q[m].mean() if m.any() else np.nan, 'obs_rate': seen[m].mean() if m.any() else np.nan})
cal = pd.DataFrame(rows); print(cal.round(3).to_string()); cal.round(4).to_csv(f'{OUT}/calibration.tsv', sep='\t', index=False)
# high-confidence prior on few bins
p = fz['p']; nb = seen.sum(0)
m = (p > 0.99) & (fz['PI'] > 0.5)
print('cells with p_tax>0.99 & P_intact>0.5:', int(m.sum()), 'seen rate', round(seen[m].mean(), 3),
      '| of these, family in <=3 bins:', int((m & (nb[None, :] <= 3)).sum()), 'seen rate', round(seen[m & (nb[None, :] <= 3)].mean(), 3))
# system-level presence where the main family is absent
g = pd.read_csv(f'{D}/catalogue/fam30/genes.tsv', sep='\t', usecols=['bin', 'product'])
mem = pd.read_csv(f'{D}/catalogue/final/chesterfield_members.tsv', sep='\t', usecols=['bin', 'family'])
for fam, pat in [('CF00033', r'ATP synthase subunit'), ('CF00005', r'NADH-quinone oxidoreductase subunit'),
                 ('CF00000', r'ribosomal subunit protein u[LS](?:3|4|5|14|22)$')]:
    has = set(mem[mem.family == fam].bin)
    nsys = g[g['product'].fillna('').str.contains(pat, regex=True)].groupby('bin').size().reindex(bins, fill_value=0)
    ab = [b for b in bins if b not in has]
    print(fam, 'absent in', len(ab), 'bins; of these with >=3 system genes elsewhere:', int((nsys[ab] >= 3).sum()),
          '| present bins with >=3:', int((nsys[list(has)] >= 3).sum()), 'of', len(has))
# within-bin paired comparison: mean presence bits per locus, off-population vs core-like contigs
from scipy.stats import wilcoxon
loci = pd.read_csv(f'{OUT}/loci_span.tsv', sep='\t')
st = pd.read_csv(f'{OUT}/states.tsv', sep='\t', usecols=['bin', 'family', 'state', 'bits'])
st = st[st.state != 'absent']
loci = loci.merge(st[['bin', 'family', 'bits']], on=['bin', 'family'])
cc = pd.read_csv(f'{D}/coverage_test/allbins/allbins_contigs.tsv', sep='\t', usecols=['contigName', 'bin', 'class'])
loci = loci.merge(cc, left_on=['bin', 'contig'], right_on=['bin', 'contigName'], how='left')
loci['cls'] = loci['class'].fillna('').str.replace('+low-r', '', regex=False)
pv = loci[loci.cls.isin(['off-population', 'core-like'])].groupby(['bin', 'cls']).bits.mean().unstack().dropna()
w = wilcoxon(pv['off-population'], pv['core-like'])
print('bins with both classes', len(pv), 'off-pop > core-like in', int((pv['off-population'] > pv['core-like']).sum()),
      'median diff', round(float((pv['off-population'] - pv['core-like']).median()), 2), 'wilcoxon p', w.pvalue)
for c_ in ['no-reads', 'other-bin', 'accessory-like', 'high-copy']:
    print(c_, 'loci', int((loci.cls == c_).sum()), 'mean bits', round(loci[loci.cls == c_].bits.mean(), 2))
