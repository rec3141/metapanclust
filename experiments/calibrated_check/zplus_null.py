"""Null for dZ+ under cleaning: per changed bin, drop random contigs (any status) until the removed bp
matches the real removal, 200 draws; compare the real dZ+ with that size-matched null. Same prior and
detection conventions as zplus_cleaned.py (original r_j=0 fit, G_j held fixed)."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/data/dev/metapanclust/experiments/fragmentation')
from detection import p_intact, p_lost
D = '/data/scratch/metapanclust_data'; CS = f'{D}/catalogue_surprise'; OUT = f'{D}/calibrated_check'
z = np.load(f'{CS}/inputs.npz'); fz = np.load(f'{CS}/fit_flat.npz')
span = z['span']; fams = z['fams']; bins = list(z['bins']); fi = {f: i for i, f in enumerate(fams)}
meta = pd.read_csv(f'{CS}/bins.tsv', sep='\t'); c = meta.compl.values
P, EPS = fz['p_r0'], fz['eps_r0']
dec = pd.read_csv(f'{D}/bin_cleaning/decisions.tsv', sep='\t')
loci = pd.read_csv(f'{CS}/loci_span.tsv', sep='\t')
rng = np.random.default_rng(1)

def zrow(j, keep, d, lb):
    """Z+ for bin j given a boolean keep mask over contigs d (rows of decisions.tsv)."""
    kc = set(d.contigName[keep]); Ls = d.len[keep].values.astype(float)
    G = d.len.sum()/max(c[j], 1e-3); cj = Ls.sum()/G
    o = np.zeros(len(fams), int)
    for (f), g in lb[lb.contig.isin(kc)].groupby('family'):
        o[fi[f]] = 2 if (~g.at_contig_end).any() else 1
    PI = p_intact(span, Ls, cj); PL = p_lost(span, Ls, cj); FS = p_intact(span, Ls, 1.0)
    p, eps = P[j], EPS[j]
    q = np.stack([p*PL+(1-p)*(1-eps), p*np.clip(1-PI-PL, 0, 1)+(1-p)*eps*(1-FS), p*PI+(1-p)*eps*FS])
    q = np.clip(q, 1e-12, 1); q /= q.sum(0); L = -np.log2(q)
    s = L[o, np.arange(len(o))]; h = (q*L).sum(0); v = (q*L**2).sum(0)-h**2; lo = (1-q[0]) < 0.5
    return ((s-h)*lo).sum()/np.sqrt((v*lo).sum())

rows = []
for b in bins:
    d = dec[dec.bin == b].reset_index(drop=True); rm = (d.action == 'remove').values
    if rm.sum() == 0: continue
    j = bins.index(b); lb = loci[loci.bin == b]
    z0 = zrow(j, np.ones(len(d), bool), d, lb); z1 = zrow(j, ~rm, d, lb)
    target = d.len[rm].sum(); nul = []
    for _ in range(200):
        order = rng.permutation(len(d)); cum = np.cumsum(d.len.values[order])
        k = np.searchsorted(cum, target) + 1; keep = np.ones(len(d), bool); keep[order[:k]] = False
        if keep.sum() == 0: keep[order[k-1]] = True
        nul.append(zrow(j, keep, d, lb) - z0)
    nul = np.array(nul)
    rows.append((b, c[j], meta.contam[j], round(rm.mean(), 3), round(target/d.len.sum(), 3), z0, z1, z1-z0, nul.mean(), nul.std(),
                 ((z1-z0)-nul.mean())/nul.std(), (nul <= z1-z0).mean()))
out = pd.DataFrame(rows, columns=['bin', 'compl', 'contam', 'rm_contig_frac', 'rm_bp_frac', 'Zp_orig', 'Zp_clean', 'dZp',
                                  'null_mean', 'null_sd', 'z_vs_null', 'P_null_le'])
out.round(3).to_csv(f'{OUT}/zplus_null.tsv', sep='\t', index=False)
pd.set_option('display.width', 200)
print(out.sort_values('rm_bp_frac', ascending=False).round(2).to_string(index=False))
from scipy.stats import spearmanr, wilcoxon
x = out.dZp - out.null_mean
print('changed bins %d: median (dZp - null mean) %.2f; real < null mean in %d; one-sided P_null_le<0.05 in %d' % (
    len(out), x.median(), (x < 0).sum(), (out.P_null_le < 0.05).sum()))
print('Wilcoxon dZp vs null mean: p=%.2g' % wilcoxon(out.dZp, out.null_mean).pvalue)
print('rho(rm_bp_frac, dZp-null) = %.2f' % spearmanr(out.rm_bp_frac, x)[0])
h = out.compl >= 0.5
print('compl>=0.5: n=%d, median excess drop %.2f, real<null in %d, P<0.05 in %d' % (h.sum(), x[h].median(), (x[h] < 0).sum(), (out.P_null_le[h] < 0.05).sum()))
