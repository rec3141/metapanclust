"""Z+ (calibrated presence-direction excess, r_j = 0) on original vs coverage-cleaned bins, and per-contig
presence excess of removed vs kept contigs. Prior p_r0/eps_r0 is the ORIGINAL leave-one-bin-out fit
(fit_flat.npz), not refitted. Detection (PI/PL/FS) is recomputed on the cleaned contig set with the
genome size G_j = sum(L)/c_j held at its original value (no CheckM2 rerun)."""
import sys, numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
sys.path.insert(0, '/data/dev/metapanclust/experiments/fragmentation')
from detection import p_intact, p_lost
D = '/data/scratch/metapanclust_data'; CS = f'{D}/catalogue_surprise'; OUT = f'{D}/calibrated_check'
z = np.load(f'{CS}/inputs.npz'); fz = np.load(f'{CS}/fit_flat.npz')
obs = z['obs'].astype(int); span = z['span']; fams = z['fams']; bins = z['bins']; B, F = obs.shape
meta = pd.read_csv(f'{CS}/bins.tsv', sep='\t'); c = meta.compl.values
p, eps = fz['p_r0'], fz['eps_r0']
dec = pd.read_csv(f'{D}/bin_cleaning/decisions.tsv', sep='\t')
loci = pd.read_csv(f'{CS}/loci_span.tsv', sep='\t')
fi = {f: i for i, f in enumerate(fams)}; bi = {b: i for i, b in enumerate(bins)}
removed = set(zip(dec.bin[dec.action == 'remove'], dec.contigName[dec.action == 'remove']))
clen = dict(zip(zip(dec.bin, dec.contigName), dec.len))
loci['rm'] = [(b, k) in removed for b, k in zip(loci.bin, loci.contig)]

def score(obs, PI, PL, FS):
    q = np.stack([p*PL+(1-p)*(1-eps), p*(np.clip(1-PI-PL, 0, 1))+(1-p)*eps*(1-FS), p*PI+(1-p)*eps*FS])
    q = np.clip(q, 1e-12, 1); q /= q.sum(0)
    L = -np.log2(q); s = np.take_along_axis(L, obs[None], 0)[0]; h = (q*L).sum(0); v = (q*L**2).sum(0)-h**2
    e = s-h; lo = (1-q[0]) < 0.5
    return e, lo, (e*lo).sum(1), (e*lo).sum(1)/np.sqrt((v*lo).sum(1))

def detect(Ls, G):
    PI = np.zeros((B, F)); PL = np.zeros((B, F)); FS = np.zeros((B, F))
    for j in range(B):
        cj = Ls[j].sum()/G[j]
        PI[j] = p_intact(span, Ls[j], cj); PL[j] = p_lost(span, Ls[j], cj); FS[j] = p_intact(span, Ls[j], 1.0)
    return PI, PL, FS

L0 = [z[f'L_{i}'] for i in range(B)]; G = np.array([l.sum() for l in L0])/np.clip(c, 1e-3, 1)
e0, lo0, Ep0, Z0 = score(obs, *detect(L0, G))
# cleaned obs: rebuild from loci on kept contigs
obs1 = np.zeros_like(obs)
for (b, f), d in loci[~loci.rm].groupby(['bin', 'family']):
    obs1[bi[b], fi[f]] = 2 if (~d.at_contig_end).any() else 1
L1 = []
for j, b in enumerate(bins):
    d = dec[dec.bin == b]
    L1.append(d.len[d.action == 'keep'].values.astype(float) if len(d) else L0[j])
L1 = [l if len(l) else np.array([1.0]) for l in L1]
e1, lo1, Ep1, Z1 = score(obs1, *detect(L1, G))
rmMb = np.array([dec.len[(dec.bin == b) & (dec.action == 'remove')].sum()/1e6 for b in bins])
totMb = np.array([l.sum()/1e6 for l in L0])
cc = pd.read_csv(f'{D}/bin_cleaning/checkm2_compare.tsv', sep='\t').set_index('bin')
out = pd.DataFrame({'bin': bins, 'compl': c, 'contam': meta.contam, 'Mb': totMb.round(2), 'removed_Mb': rmMb.round(3),
                    'removed_frac': (rmMb/totMb).round(3), 'Zp_orig': Z0, 'Zp_clean': Z1, 'dZp': Z1-Z0,
                    'Ep_orig': Ep0, 'Ep_clean': Ep1, 'dEp': Ep1-Ep0,
                    'dConta_ckm2': cc.reindex(bins).dConta.values})
out.round(3).to_csv(f'{OUT}/zplus_cleaned.tsv', sep='\t', index=False)
pd.set_option('display.width', 200)
print(out.sort_values('removed_frac', ascending=False).head(15).round(2).to_string(index=False))
ch = out.removed_Mb > 0
print('bins changed', ch.sum(), 'unchanged', (~ch).sum())
print('unchanged bins: max |dZp| %.3g' % out.dZp[~ch].abs().max())
for nm, x in [('removed_frac', out.removed_frac), ('removed_Mb', out.removed_Mb)]:
    r = spearmanr(x[ch], out.dZp[ch]); print('changed bins: rho(%s, dZp) = %.2f (p=%.2g); rho(%s, dEp) = %.2f' % (nm, r[0], r[1], nm, spearmanr(x[ch], out.dEp[ch])[0]))
r = spearmanr(out.removed_frac[ch], -out.dConta_ckm2[ch]); print('changed bins: rho(removed_frac, CheckM2 conta drop) = %.2f' % r[0])
r = spearmanr(out.dZp[ch], out.dConta_ckm2[ch]); print('changed bins: rho(dZp, dConta_ckm2) = %.2f' % r[0])
print('median dZp changed %.2f, Zp_orig median changed %.2f vs unchanged %.2f' % (out.dZp[ch].median(), out.Zp_orig[ch].median(), out.Zp_orig[~ch].median()))
m = out.compl >= 0.7; y = out.contam >= 0.1
def auc(x, y): return mannwhitneyu(x[y], x[~y]).statistic/(y.sum()*(~y).sum())
print('AUC Zp vs CheckM2 contam>=0.1 (compl>=0.7): orig %.3f' % auc(out.Zp_orig[m].values, y[m].values))
print('AUC Zp_orig vs removed_frac>=0.1 (compl>=0.7): %.3f  (n pos %d)' % (auc(out.Zp_orig[m].values, (out.removed_frac >= 0.1)[m].values), (out.removed_frac >= 0.1)[m].sum()))
# ---- per contig: presence excess on original fit, cells with q(seen)<1/2 and seen; split cell e over its loci
cell = e0*lo0*(obs > 0)
nl = loci.groupby(['bin', 'family']).size()
loci['e'] = [cell[bi[b], fi[f]]/nl[(b, f)] for b, f in zip(loci.bin, loci.family)]
loci['eall'] = [e0[bi[b], fi[f]]*(obs[bi[b], fi[f]] > 0)/nl[(b, f)] for b, f in zip(loci.bin, loci.family)]
ce = loci.groupby(['bin', 'contig']).agg(Eplus=('e', 'sum'), Eall=('eall', 'sum'), nloci=('e', 'size')).reset_index()
dd = dec[dec.bin.isin(bins)][['bin', 'contigName', 'len', 'action']].rename(columns={'contigName': 'contig'})
ce = dd.merge(ce, on=['bin', 'contig'], how='left').fillna({'Eplus': 0, 'Eall': 0, 'nloci': 0})
ce['Eplus_perkb'] = ce.Eplus/ce.len*1e3
ce.to_csv(f'{OUT}/contig_excess.tsv', sep='\t', index=False)
# length-matched within-bin comparison: each removed contig vs kept contigs in the same bin within 2x length
wins = ties = n = 0; rows = []
for b, d in ce.groupby('bin'):
    R = d[d.action == 'remove']; K = d[d.action == 'keep']
    for _, r in R.iterrows():
        k = K[(K.len >= r.len/2) & (K.len <= r.len*2)]
        if len(k) == 0: continue
        wins += (r.Eplus_perkb > k.Eplus_perkb).sum(); ties += (r.Eplus_perkb == k.Eplus_perkb).sum(); n += len(k)
        rows.append((b, r.Eplus_perkb - k.Eplus_perkb.median()))
print('per-contig: removed %d, kept %d contigs in catalogue bins' % ((ce.action == 'remove').sum(), (ce.action == 'keep').sum()))
print('length-matched (2x) within-bin pairs %d: P(removed > kept) = %.3f (ties %.3f), AUC=(w+t/2)/n = %.3f' % (n, wins/n, ties/n, (wins+ties/2)/n))
dl = pd.DataFrame(rows, columns=['bin', 'diff'])
print('removed contigs with a match %d; median (removed - median matched kept) E+/kb = %.3f; in %d/%d bins median diff > 0' % (
    len(dl), dl['diff'].median(), (dl.groupby('bin')['diff'].median() > 0).sum(), dl.bin.nunique()))
for lab, sel in [('<5 kb', ce.len < 5e3), ('5-20 kb', (ce.len >= 5e3) & (ce.len < 2e4)), ('>=20 kb', ce.len >= 2e4)]:
    a = ce[sel]
    print('%-8s removed n=%d median E+/kb %.3f mean %.3f | kept n=%d median %.3f mean %.3f' % (lab, (a.action == 'remove').sum(),
          a.Eplus_perkb[a.action == 'remove'].median(), a.Eplus_perkb[a.action == 'remove'].mean(), (a.action == 'keep').sum(),
          a.Eplus_perkb[a.action == 'keep'].median(), a.Eplus_perkb[a.action == 'keep'].mean()))
for b in ['dastool-maxbin_102', 'dastool-vamb_711_sub', 'dastool-maxbin_040', 'dastool-lorbin_015', 'dastool-maxbin_115']:
    d = ce[ce.bin == b]
    print('%-22s E+ removed %.1f bits on %.2f Mb (%.1f/Mb) | kept %.1f on %.2f Mb (%.1f/Mb)' % (b, d.Eplus[d.action == 'remove'].sum(),
          d.len[d.action == 'remove'].sum()/1e6, d.Eplus[d.action == 'remove'].sum()/max(d.len[d.action == 'remove'].sum()/1e6, 1e-9),
          d.Eplus[d.action == 'keep'].sum(), d.len[d.action == 'keep'].sum()/1e6, d.Eplus[d.action == 'keep'].sum()/(d.len[d.action == 'keep'].sum()/1e6)))
