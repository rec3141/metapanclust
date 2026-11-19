"""Follow-ups: absence surprise as a contamination flag at matched completeness; off-population loci
without TMED14; system genes in the bins with the strongest atp/nuo absences."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
D = '/data/scratch/metapanclust_data'; OUT = f'{D}/catalogue_surprise'
per = pd.read_csv(f'{OUT}/per_bin.tsv', sep='\t')
hi = per[per.compl >= 0.7]
print('compl>=0.7 bins', len(hi), 'spearman S_abs_per_expMb~contam', spearmanr(hi.contam, hi.S_abs_per_expMb),
      '| S_pres_per_Mb~contam', spearmanr(hi.contam, hi.S_pres_per_Mb))
c10 = hi.contam >= 0.10
print('compl>=0.7: contam>=10%:', int(c10.sum()), 'median S_abs/expMb', hi[c10].S_abs_per_expMb.median().round(1),
      'vs', hi[~c10].S_abs_per_expMb.median().round(1), 'MWU p', mannwhitneyu(hi[c10].S_abs_per_expMb, hi[~c10].S_abs_per_expMb).pvalue,
      '| AUC', round(mannwhitneyu(hi[c10].S_abs_per_expMb, hi[~c10].S_abs_per_expMb).statistic / (c10.sum() * (~c10).sum()), 2),
      '| presence AUC', round(mannwhitneyu(hi[c10].S_pres_per_Mb, hi[~c10].S_pres_per_Mb).statistic / (c10.sum() * (~c10).sum()), 2))
print('top-10 absence bins with contam>=10%:', int((per.sort_values('S_abs_per_expMb', ascending=False).head(10).contam >= 0.1).sum()),
      'of 10; all bins with contam>=10%:', int((per.contam >= 0.1).sum()), 'of', len(per))
loci = pd.read_csv(f'{OUT}/loci_span.tsv', sep='\t')
st = pd.read_csv(f'{OUT}/states.tsv', sep='\t', usecols=['bin', 'family', 'state', 'bits'])
loci = loci.merge(st[st.state != 'absent'][['bin', 'family', 'bits']], on=['bin', 'family'])
cc = pd.read_csv(f'{D}/coverage_test/allbins/allbins_contigs.tsv', sep='\t', usecols=['contigName', 'bin', 'class'])
loci = loci.merge(cc, left_on=['bin', 'contig'], right_on=['bin', 'contigName'], how='left')
loci['cls'] = loci['class'].fillna('').str.replace('+low-r', '', regex=False)
x = loci[(loci.cls == 'off-population') & (loci.bin != 'dastool-maxbin_102')].bits
y = loci[(loci.cls == 'core-like') & (loci.bin != 'dastool-maxbin_102')].bits
print('without TMED14: off-pop loci', len(x), 'median', round(x.median(), 2), 'mean', round(x.mean(), 2), '| core-like median', round(y.median(), 2),
      'mean', round(y.mean(), 2), 'MWU p', mannwhitneyu(x, y).pvalue)
print('share of loci >5 bits: off-pop', round((loci[loci.cls == 'off-population'].bits > 5).mean(), 3), 'core-like', round((loci[loci.cls == 'core-like'].bits > 5).mean(), 3),
      'no-reads', round((loci[loci.cls == 'no-reads'].bits > 5).mean(), 3), 'accessory-like', round((loci[loci.cls == 'accessory-like'].bits > 5).mean(), 3))
g = pd.read_csv(f'{D}/catalogue/fam30/genes.tsv', sep='\t', usecols=['bin', 'product', 'family'])
for b, pat in [('dastool-vamb_623', 'ATP synthase'), ('dastool-vamb_049', 'NADH-quinone oxidoreductase|NADH dehydrogenase'),
               ('dastool-lorbin_029', 'ribosomal subunit protein')]:
    h = g[(g.bin == b) & g['product'].fillna('').str.contains(pat)]
    print(b, pat, len(h), sorted(h['product'].str.replace(' subunit', '').str.slice(0, 40).unique())[:14])
mem = pd.read_csv(f'{D}/catalogue/final/chesterfield_members.tsv', sep='\t', usecols=['bin', 'family', 'products'])
for b, pat in [('dastool-vamb_623', 'ATP synthase'), ('dastool-vamb_049', 'NADH-quinone')]:
    m = mem[(mem.bin == b) & mem.products.str.contains(pat)]
    print(b, 'loci with', pat, m.family.tolist())
