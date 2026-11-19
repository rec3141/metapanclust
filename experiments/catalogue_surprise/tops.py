"""Summaries of the most surprising family states for the report."""
import pandas as pd
OUT = '/data/scratch/metapanclust_data/catalogue_surprise'
st = pd.read_csv(f'{OUT}/states.tsv', sep='\t')
pr = st[st.state != 'absent']; ab = st[st.state == 'absent']
hyp = lambda s: all(x.strip() == 'hypothetical protein' for x in str(s).split('|'))
print('present states', len(pr), '>5 bits', int((pr.bits > 5).sum()), 'of which singleton families', int(((pr.bits > 5) & (pr.n_bins_family == 1)).sum()),
      'all-hypothetical', int(pr[pr.bits > 5].medoid_products.map(hyp).sum()))
print('multi-bin families >5 bits:', int(((pr.bits > 5) & (pr.n_bins_family >= 2)).sum()))
print('absent states >1 bit', len(ab), '>2 bits', int((ab.bits > 2).sum()), 'max', ab.bits.max().round(2))
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 110)
print(ab.sort_values('bits', ascending=False).head(15)[['bin', 'taxon', 'family', 'p_tax', 'bits', 'n_bins_family', 'medoid_products']].round(2).to_string())
m = pr[pr.n_bins_family >= 3].sort_values('bits', ascending=False)
m = m[~m.medoid_products.map(hyp)]
print(m.head(20)[['bin', 'taxon', 'family', 'state', 'p_tax', 'p_flat', 'bits', 'n_bins_family', 'medoid_products']].round(3).to_string())
