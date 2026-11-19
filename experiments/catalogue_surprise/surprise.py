"""Level-1 leave-one-bin-out surprise over catalogue families (Chesterfield, 93 bins).
Hierarchical empirical-Bayes beta-binomial along the GTDB lineage (model_spec.md 3.4 A), observation
model with three states I/P/A (model_spec.md 2, D1; detection.py for P_intact/P_lost)."""
import os, sys, json, numpy as np, pandas as pd
sys.path.insert(0, '/data/dev/metapanclust/experiments/fragmentation')
from detection import p_intact, p_lost

D = '/data/scratch/metapanclust_data'
OUT = f'{D}/catalogue_surprise'
RANKS = ['d', 'p', 'c', 'o', 'f', 'g', 's']
EPS0 = 1e-4
EM_MODE = os.environ.get('EM_MODE', 'flat')
GRID = [0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 1e4]

z = np.load(f'{OUT}/inputs.npz')
obs, nloc, span, fams, bins = z['obs'], z['nloc'], z['span'], z['fams'], z['bins']
B, F = obs.shape
L = [z[f'L_{i}'] for i in range(B)]
meta = pd.read_csv(f'{OUT}/bins.tsv', sep='\t')
c = meta.compl.values

# clade ids per rank; an unnamed rank (e.g. 's__') makes a node of its own for the bin
cid = {}
for r_i, r in enumerate(RANKS):
    keys = []
    for b, lin in zip(bins, meta.lineage):
        parts = lin.split(';')
        keys.append(';'.join(parts[:r_i + 1]) if parts[r_i][3:] else f'{b}@{r}')
    u, inv = np.unique(keys, return_inverse=True)
    cid[r] = inv

# detection given presence, per bin x family
PI = np.zeros((B, F)); PL = np.zeros((B, F)); FS = np.zeros((B, F))
for j in range(B):
    PI[j] = p_intact(span, L[j], c[j])
    PL[j] = p_lost(span, L[j], c[j])
    FS[j] = p_intact(span, L[j], 1.0)
PP = np.clip(1 - PI - PL, 0, 1)
Lsum = np.array([l.sum() for l in L], float)
G = Lsum / np.clip(c, 1e-3, 1)


def false_presence(r):
    """eps_jk = eps0 + 1 - exp(-r_j G_j rho_k), rho_k = loci of k per bp in the other bins."""
    rho = (nloc.sum(0)[None, :] - nloc) / (Lsum.sum() - Lsum)[:, None]
    eps = EPS0 + 1 - np.exp(-(r * G)[:, None] * rho)
    return np.minimum(eps, 0.5)


def obs_probs(eps):
    a = np.where(obs == 2, PI, np.where(obs == 1, PP, PL))
    b = np.where(obs == 2, eps * FS, np.where(obs == 1, eps * (1 - FS), 1 - eps))
    return np.maximum(a, 1e-12), np.maximum(b, 1e-12)


def loo_prior(post, kappa, a0, b0, ranks=RANKS):
    X = post.sum(0)[None, :]
    p = (X - post + a0) / (B - 1 + a0 + b0)
    for r in ranks:
        k = kappa[r]
        ids = cid[r]
        n = ids.max() + 1
        M = np.zeros((n, B)); M[ids, np.arange(B)] = 1
        Xc = (M @ post)[ids] - post
        Nc = M.sum(1)[ids][:, None] - 1
        p = (Xc + k * p) / (Nc + k)
    return np.clip(p, 1e-6, 1 - 1e-6)


def fit(a, b, kappa, ranks, mode=None):
    """Training counts per bin: P(z=1 | O_g). The held-out bin j must not reach its own prediction
    through the other bins' posteriors, so those posteriors may use only priors that do not depend on
    j's clade: 'none' = observed 0/1 (no detection correction of counts), 'flat' = posterior under the
    flat (root) prior, which carries j at weight 1/(B-1). 'full' = iterate with the hierarchical LOO
    prior (model_spec 3.4 A); leaks j through its clade-mates' posteriors, kept only as a contrast."""
    mode = mode or EM_MODE
    m = (obs > 0).mean(0); v = m.var(); mu = m.mean()
    s_ab = max(mu * (1 - mu) / v - 1, 0.05)
    a0, b0 = mu * s_ab, (1 - mu) * s_ab
    post = (obs > 0).astype(float)
    if mode == 'flat':
        for _ in range(3):
            p = loo_prior(post, {r: 1e9 for r in RANKS}, a0, b0, RANKS)
            post = p * a / (p * a + (1 - p) * b)
    elif mode == 'full':
        for _ in range(3):
            p = loo_prior(post, kappa, a0, b0, ranks)
            post = p * a / (p * a + (1 - p) * b)
    return loo_prior(post, kappa, a0, b0, ranks), (a0, b0)


def logloss(p, a, b, mask=None):
    q = p * a + (1 - p) * b
    ll = -np.log2(q)
    return float(ll.mean() if mask is None else ll[:, mask].mean())


def tune(a, b, fmask):
    kappa = {r: 10.0 for r in RANKS}
    for sweep in range(2):
        for r in reversed(RANKS):
            best = None
            for k in GRID:
                kk = dict(kappa); kk[r] = k
                p, _ = fit(a, b, kk, RANKS)
                ll = logloss(p, a, b, fmask)
                if best is None or ll < best[0]:
                    best = (ll, k)
            kappa[r] = best[1]
    return kappa


def collapse_ll(p, a_I, a_P, a_A, eps, mask=None):
    """2-state (seen vs not) log loss, for comparison with a model that ignores detection."""
    q = p * (a_I + a_P) + (1 - p) * eps
    seen = obs > 0
    ll = -np.log2(np.where(seen, q, 1 - q))
    return float(ll.mean() if mask is None else ll[:, mask].mean())


if __name__ == '__main__':
    res = {}
    eps = false_presence(meta.contam.values)
    a, b = obs_probs(eps)
    even = np.arange(F) % 2 == 0
    # kappa tuned on even families, evaluated on odd (and vice versa) -> honest held-out numbers
    k_even = tune(a, b, even); k_odd = tune(a, b, ~even)
    res['kappa_even'] = k_even; res['kappa_odd'] = k_odd
    kappa = tune(a, b, None); res['kappa_all'] = kappa
    flat = {r: 1e9 for r in RANKS}  # infinite kappa = every clade equals the root = flat prior
    hard = {r: 1e-3 for r in RANKS}
    models = {'flat': flat, 'tax_hard': hard, 'tax_shrunk': kappa}
    multi = (obs > 0).sum(0) >= 2
    rows = []
    P = {}
    for name, kk in models.items():
        p, ab = fit(a, b, kk, RANKS)
        P[name] = p
        row = {'model': name, 'll_all': logloss(p, a, b), 'll_multi_bin_fams': logloss(p, a, b, multi)}
        q = p * a + (1 - p) * b
        ll = -np.log2(q)
        row['ll_present'] = float(ll[obs > 0].mean()); row['ll_absent'] = float(ll[obs == 0].mean())
        row['bits_total'] = float(ll.sum())
        rows.append(row)
    # cross-fitted kappa
    pe, _ = fit(a, b, k_even, RANKS); po, _ = fit(a, b, k_odd, RANKS)
    pc = np.where(even[None, :], po, pe)
    q = pc * a + (1 - pc) * b; ll = -np.log2(q)
    rows.append({'model': 'tax_shrunk_crossfit', 'll_all': float(ll.mean()), 'll_multi_bin_fams': float(ll[:, multi].mean()),
                 'll_present': float(ll[obs > 0].mean()), 'll_absent': float(ll[obs == 0].mean()), 'bits_total': float(ll.sum())})
    # flat with no detection model: P(seen)=p, 2-state
    for name in ['flat', 'tax_shrunk']:
        p = P[name]
        seen = obs > 0
        pn = np.clip(seen.sum(0)[None, :] - seen, 0, None)
        q2 = p * (PI + PP) + (1 - p) * eps
        ll = -np.log2(np.where(seen, q2, 1 - q2))
        rows.append({'model': name + '_2state_with_detection', 'll_all': float(ll.mean()), 'll_multi_bin_fams': float(ll[:, multi].mean()),
                     'll_present': float(ll[seen].mean()), 'll_absent': float(ll[~seen].mean())})
    # 2-state model with no detection: LOO beta-binomial on raw seen/not
    seen = (obs > 0).astype(float)
    mu = seen.mean(); s_ab = max(mu * (1 - mu) / seen.mean(0).var() - 1, 0.05)

    def nd_ll(kk, mask=None):
        p = loo_prior(seen, kk, mu * s_ab, (1 - mu) * s_ab)
        ll = -np.log2(np.where(seen > 0, p, 1 - p))
        return ll if mask is None else ll[:, mask]
    k_nd = {r: 10.0 for r in RANKS}
    for sweep in range(2):
        for r in reversed(RANKS):
            k_nd[r] = min(GRID, key=lambda k: nd_ll({**k_nd, r: k}).mean())
    res['kappa_nodetection'] = k_nd
    for name, kk in [('flat', flat), ('tax_shrunk', k_nd)]:
        ll = nd_ll(kk)
        rows.append({'model': name + '_2state_no_detection', 'll_all': float(ll.mean()), 'll_multi_bin_fams': float(ll[:, multi].mean()),
                     'll_present': float(ll[seen > 0].mean()), 'll_absent': float(ll[seen == 0].mean())})
    pd.DataFrame(rows).to_csv(f'{OUT}/logloss_{EM_MODE}.tsv', sep='\t', index=False)
    print(pd.DataFrame(rows).to_string())
    print(json.dumps(res))
    json.dump(res, open(f'{OUT}/kappa_{EM_MODE}.json', 'w'))
    # contamination-agnostic false-presence (r_j = 0) for the contamination check
    eps0 = false_presence(np.zeros(B)); a0_, b0_ = obs_probs(eps0)
    p0, _ = fit(a0_, b0_, kappa, RANKS)
    np.savez_compressed(f'{OUT}/fit_{EM_MODE}.npz', p=P['tax_shrunk'], p_flat=P['flat'], a=a, b=b, eps=eps, PI=PI, PL=PL, PP=PP, FS=FS,
                        p_r0=p0, a_r0=a0_, b_r0=b0_, eps_r0=eps0)
