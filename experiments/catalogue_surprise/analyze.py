"""Per-bin surprise, top family states with Bakta products, contamination/misbinning checks,
core-operon checks. Reads fit.npz from surprise.py."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

D = '/data/scratch/metapanclust_data'
OUT = f'{D}/catalogue_surprise'
z = np.load(f'{OUT}/inputs.npz'); fz = np.load(f'{OUT}/fit_flat.npz')
obs, nloc, span, fams, bins = z['obs'], z['nloc'], z['span'], z['fams'], z['bins']
meta = pd.read_csv(f'{OUT}/bins.tsv', sep='\t')
cat = pd.read_csv(f'{D}/catalogue/final/chesterfield_catalogue.tsv', sep='\t', usecols=['family', 'n_bins', 'medoid_products']).set_index('family')
prod = cat.loc[fams, 'medoid_products'].fillna('').str.slice(0, 160).values
nb = cat.loc[fams, 'n_bins'].values
B, F = obs.shape


def states(p, a_, PI, PP, PL, eps, FS):
    """3-state probabilities for every (bin, family) and surprise of the observed state + entropy."""
    q = np.stack([p * PL + (1 - p) * (1 - eps), p * PP + (1 - p) * eps * (1 - FS), p * PI + (1 - p) * eps * FS])
    q = np.clip(q, 1e-12, 1)
    s = -np.log2(np.take_along_axis(q, obs[None].astype(int), 0)[0])
    h = -(q * np.log2(q)).sum(0)
    return s, h


tag = {}
for v, (pk, ek) in {'main': ('p', 'eps'), 'r0': ('p_r0', 'eps_r0')}.items():
    s, h = states(fz[pk], None, fz['PI'], fz['PP'], fz['PL'], fz[ek], fz['FS'])
    tag[v] = (s, h)
s, h = tag['main']; s0, _ = tag['r0']
s_det, h_det = s, h


def content(pk, ak, bk):
    """presence: -log2 P(O|p)/P(O|z=1) = -log2(p + (1-p) b/a); 0 when p = 1, -log2 p when b = 0.
    absence keeps the full surprise -log2 P(A|p) (P(A|z=0) ~ 1, so it is already content-only)."""
    pp = fz[pk]
    return -np.log2(pp + (1 - pp) * fz[bk] / fz[ak])
c_main = content('p', 'a', 'b'); c_r0 = content('p_r0', 'a_r0', 'b_r0')
pf_ = fz['p_flat']; c_flat = -np.log2(pf_ + (1 - pf_) * fz['b'] / fz['a'])
pres_ = obs > 0
s = np.where(pres_, c_main, s_det); s0 = np.where(pres_, c_r0, s0)
p = fz['p']; pf = fz['p_flat']
sf, _ = states(pf, None, fz['PI'], fz['PP'], fz['PL'], fz['eps'], fz['FS'])
sf = np.where(pres_, c_flat, sf)
pres = obs > 0
Mb = meta.size_bp.values / 1e6
exp_Mb = (p * span[None, :]).sum(1) / 1e6
short = meta.lineage.str.split(';').apply(lambda x: next((t for t in reversed(x) if len(t) > 3), '')).values
per = pd.DataFrame({
    'bin': bins, 'taxon': short, 'compl': meta.compl, 'contam': meta.contam, 'Mb': Mb.round(2),
    'n_present': pres.sum(1), 'n_singleton_present': (pres & (nb[None, :] == 1)).sum(1),
    'S_pres_bits': (s * pres).sum(1), 'S_pres_per_Mb': (s * pres).sum(1) / Mb,
    'S_pres_det_per_Mb': (s_det * pres).sum(1) / Mb,
    'S_pres_multi_per_Mb': (s * pres * (nb[None, :] >= 2)).sum(1) / Mb,
    'S_pres_r0_per_Mb': (s0 * pres).sum(1) / Mb,
    'S_pres_flat_per_Mb': (sf * pres).sum(1) / Mb,
    'S_abs_bits': (s * ~pres).sum(1), 'S_abs_per_expMb': (s * ~pres).sum(1) / exp_Mb,
    'n_pres_gt5bits': ((s > 5) & pres).sum(1), 'n_abs_gt2bits': ((s > 2) & ~pres).sum(1)})
per.round(4).to_csv(f'{OUT}/per_bin.tsv', sep='\t', index=False)

# long table of notable states
rows = []
for j in range(B):
    for k in np.where(pres[j])[0]:
        rows.append((bins[j], short[j], fams[k], 'present_intact' if obs[j, k] == 2 else 'present_partial', p[j, k], pf[j, k], s[j, k], nb[k], span[k], prod[k]))
    ab = np.where(~pres[j] & (s[j] > 1))[0]
    for k in ab:
        rows.append((bins[j], short[j], fams[k], 'absent', p[j, k], pf[j, k], s[j, k], nb[k], span[k], prod[k]))
st = pd.DataFrame(rows, columns=['bin', 'taxon', 'family', 'state', 'p_tax', 'p_flat', 'bits', 'n_bins_family', 'span_bp', 'medoid_products'])
st.sort_values(['bin', 'bits'], ascending=[True, False]).round(4).to_csv(f'{OUT}/states.tsv', sep='\t', index=False)
top = st.sort_values('bits', ascending=False).groupby(['bin', st.state.str.startswith('present')]).head(5)
top.sort_values(['bin', 'state', 'bits'], ascending=[True, True, False]).round(4).to_csv(f'{OUT}/top_states_per_bin.tsv', sep='\t', index=False)

print('=== per-bin summary (medians)')
print(per.drop(columns=['bin', 'taxon']).median().round(3).to_string())
print('totals: present bits', round(float((s * pres).sum())), 'absent bits', round(float((s * ~pres).sum())),
      'singleton-family present bits', round(float((s * pres * (nb[None, :] == 1)).sum())))
print('=== top 12 bins by presence bits/Mb')
print(per.sort_values('S_pres_per_Mb', ascending=False).head(12)[['bin', 'taxon', 'compl', 'contam', 'Mb', 'n_present', 'S_pres_per_Mb', 'S_pres_multi_per_Mb', 'S_pres_r0_per_Mb', 'S_abs_per_expMb']].round(2).to_string())
print('=== top 8 bins by absence bits per expected Mb')
print(per.sort_values('S_abs_per_expMb', ascending=False).head(8)[['bin', 'taxon', 'compl', 'contam', 'S_abs_bits', 'S_abs_per_expMb', 'n_abs_gt2bits']].round(2).to_string())
for col in ['S_pres_per_Mb', 'S_pres_multi_per_Mb', 'S_pres_r0_per_Mb', 'S_pres_flat_per_Mb', 'S_abs_per_expMb']:
    r1 = spearmanr(per.contam, per[col]); r2 = spearmanr(per.compl, per[col])
    print(f'spearman {col}: contam rho={r1[0]:.2f} p={r1[1]:.2g}; compl rho={r2[0]:.2f} p={r2[1]:.2g}')
for b in ['dastool-maxbin_102', 'dastool-maxbin_040', 'dastool-vamb_711_sub']:
    r = per[per.bin == b].iloc[0]
    for col in ['S_pres_per_Mb', 'S_pres_multi_per_Mb', 'S_pres_r0_per_Mb']:
        print(b, col, round(r[col], 1), 'rank', int((per[col] >= r[col]).sum()), 'of', B, 'percentile', round((per[col] < r[col]).mean(), 2))

# contig classes (coverage_test): presence bits per Mb by class
cc = pd.read_csv(f'{D}/coverage_test/allbins/allbins_contigs.tsv', sep='\t', usecols=['contigName', 'len', 'bin', 'class'])
loci = pd.read_csv(f'{OUT}/loci_span.tsv', sep='\t')
bi = {b: i for i, b in enumerate(bins)}; fi = {f: i for i, f in enumerate(fams)}
jj = loci.bin.map(bi).values; kk = loci.family.map(fi).values
loci['bits'] = s[jj, kk] / nloc[jj, kk]
loci['bits_r0'] = s0[jj, kk] / nloc[jj, kk]
loci['multi'] = nb[kk] >= 2
lc = loci.groupby(['bin', 'contig']).agg(bits=('bits', 'sum'), bits_r0=('bits_r0', 'sum'), n=('bits', 'size'),
                                         hi=('bits', lambda x: (x > 5).sum())).reset_index()
cc = cc.merge(lc, left_on=['bin', 'contigName'], right_on=['bin', 'contig'], how='left').fillna({'bits': 0, 'bits_r0': 0, 'n': 0, 'hi': 0})
g = cc.groupby('class').agg(contigs=('len', 'size'), Mb=('len', lambda x: x.sum() / 1e6), loci=('n', 'sum'), bits=('bits', 'sum'), hi=('hi', 'sum'))
g['loci_per_Mb'] = g.loci / g.Mb; g['bits_per_Mb'] = g.bits / g.Mb; g['bits_per_locus'] = g.bits / g.loci; g['loci_gt5bits_per_Mb'] = g.hi / g.Mb
print('=== presence surprise by coverage_test contig class (all bins)')
print(g.round(2).to_string())
g.round(3).to_csv(f'{OUT}/contig_class.tsv', sep='\t')
# matched by length class (>=10 kb only) to avoid short-contig confounding
big = cc[cc.len >= 10000]
g2 = big.groupby('class').agg(Mb=('len', lambda x: x.sum() / 1e6), loci=('n', 'sum'), bits=('bits', 'sum'))
g2['bits_per_Mb'] = g2.bits / g2.Mb; g2['bits_per_locus'] = g2.bits / g2.loci
print('contigs >= 10 kb'); print(g2.round(2).to_string())
lb = loci.merge(cc[['bin', 'contigName', 'class']], left_on=['bin', 'contig'], right_on=['bin', 'contigName'], how='left')
x = lb[lb['class'] == 'off-population'].bits; y = lb[lb['class'] == 'core-like'].bits
print('per-locus bits off-population median', x.median().round(2), 'mean', x.mean().round(2), 'n', len(x),
      '| core-like median', y.median().round(2), 'mean', y.mean().round(2), 'n', len(y), 'MWU p', mannwhitneyu(x, y).pvalue)
for b in ['dastool-maxbin_102', 'dastool-maxbin_040']:
    t = lb[lb.bin == b].groupby('class').agg(loci=('bits', 'size'), bits_mean=('bits', 'mean'), gt5=('bits', lambda v: (v > 5).sum()))
    print(b); print(t.round(2).to_string())

# core operons
core = {'CF00000': 'S10-spc-alpha rp', 'CF00009': 'L11-L10-L12-rpoBC', 'CF00011': 'uL13-uS9', 'CF00004': 'uS2-tsf-pyrH-frr',
        'CF00030': 'alpha operon', 'CF00033': 'atp F0F1', 'CF00005': 'nuo', 'CF00017': 'groESL', 'CF00014': 'dcw'}
rows = []
for f, nm in core.items():
    k = fi[f]
    pr, ab = pres[:, k], ~pres[:, k]
    rows.append({'family': f, 'name': nm, 'n_present': int(pr.sum()), 'p_tax_present_median': np.median(p[pr, k]),
                 'bits_present_median': np.median(s[pr, k]), 'bits_present_max': s[pr, k].max(),
                 'bits_present_with_detection_median': np.median(s_det[pr, k]),
                 'n_absent': int(ab.sum()), 'p_tax_absent_median': np.median(p[ab, k]) if ab.any() else np.nan,
                 'P_intact_absent_median': np.median(fz['PI'][ab, k]) if ab.any() else np.nan,
                 'bits_absent_median': np.median(s[ab, k]) if ab.any() else np.nan,
                 'bits_absent_max': s[ab, k].max() if ab.any() else np.nan,
                 'bits_absent_flat_median': np.median(sf[ab, k]) if ab.any() else np.nan,
                 'span_kb': span[k] / 1e3})
co = pd.DataFrame(rows); co.round(3).to_csv(f'{OUT}/core_operons.tsv', sep='\t', index=False)
print('=== core operons'); print(co.round(2).to_string())
# absent-operon scores vs completeness for atp / nuo / CF00000
for f in ['CF00000', 'CF00033', 'CF00005']:
    k = fi[f]; ab = ~pres[:, k]
    d = pd.DataFrame({'bin': bins[ab], 'taxon': short[ab], 'compl': meta.compl.values[ab], 'p': p[ab, k], 'PI': fz['PI'][ab, k], 'bits': s[ab, k]})
    print(f, 'most surprising absences'); print(d.sort_values('bits', ascending=False).head(6).round(2).to_string())
    k = fi[f]; pr = pres[:, k]
    d = pd.DataFrame({'bin': bins[pr], 'taxon': short[pr], 'p': p[pr, k], 'bits': s[pr, k]})
    print(f, 'most surprising presences'); print(d.sort_values('bits', ascending=False).head(4).round(2).to_string())

# duplicated intact copies (two loci of one family, neither touching a contig end) among expected families
loci['intact'] = ~loci.at_contig_end
du = loci[loci.intact].groupby(['bin', 'family']).size().reset_index(name='n_intact')
du = du[du.n_intact >= 2]
du['p'] = p[du.bin.map(bi).values, du.family.map(fi).values]
dd = du[du.p >= 0.8].groupby('bin').size().reindex(bins, fill_value=0)
per['dup_intact_expected'] = dd.values
per.round(4).to_csv(f'{OUT}/per_bin.tsv', sep='\t', index=False)
r = spearmanr(per.contam, per.dup_intact_expected)
print('dup intact copies of expected (p>=0.8) families: median', per.dup_intact_expected.median(), 'spearman vs CheckM2 contam', round(r[0], 2), r[1])
print(per.sort_values('dup_intact_expected', ascending=False).head(8)[['bin', 'taxon', 'compl', 'contam', 'dup_intact_expected']].to_string())
