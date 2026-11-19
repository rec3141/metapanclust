"""P2 of docs/theory/phylo_prior.md on the PanGBank 2-genus test: held-out bits per test genome of
species-only EB vs genus shrinkage vs an ancestral-state-first two-state predictor, m = 1,2,3,5,10.
usage: p2_eval.py NPZ [NPZ ...] OUT.tsv"""
import sys, os
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, logit, logsumexp
W = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../pangbank_scale/worker.py')
src = open(W).read().split("\nif a.mode == 'species':")[0]
sys_argv = sys.argv; sys.argv = ['w', 'genus', 'x']
g = {'__name__': 'w'}; exec(compile(src, W, 'exec'), g); sys.argv = sys_argv
betabin_ll, genus_mean, fit_kappa, fit_beta_eb, bits = (g[k] for k in
    ('betabin_ll', 'genus_mean', 'fit_kappa', 'fit_beta_eb', 'bits'))

MS = (1, 2, 3, 5, 10); REPS = 5; NTEST = 50


def unpack(th):
    a1, b1, a0, b0 = np.exp(th[:4]); return a1, b1, a0, b0, th[4:]


def mix_ll(th, xs, ns, feats):
    """sum over (species, cluster) of log[pi BB(x; n, a1, b1) + (1 - pi) BB(x; n, a0, b0)];
    pi = sigmoid(c . feats). Component 1 = 'present in the species ancestor' (prevalence near 1),
    component 0 = 'absent in the ancestor' (prevalence near 0)."""
    a1, b1, a0, b0, c = unpack(th)
    tot = 0.0
    for x, n, F in zip(xs, ns, feats):
        lp = F @ c
        l1 = betabin_ll(x, n, a1, b1) + np.log(expit(lp)); l0 = betabin_ll(x, n, a0, b0) + np.log(expit(-lp))
        tot += np.logaddexp(l1, l0).sum()
    return -tot


def fit_mix(xs, ns, feats, nf):
    th0 = np.r_[np.log([5, 0.3, 0.05, 5]), np.zeros(nf)]
    r = minimize(mix_ll, th0, args=(xs, ns, feats), method='L-BFGS-B',
                 bounds=[(-6, 6)] * 4 + [(-15, 15)] * nf)
    return r.x


def mix_pred(th, x, m, F):
    a1, b1, a0, b0, c = unpack(th)
    lp = F @ c
    l1 = betabin_ll(x, m, a1, b1) + np.log(expit(lp)); l0 = betabin_ll(x, m, a0, b0) + np.log(expit(-lp))
    w = np.exp(l1 - np.logaddexp(l1, l0))
    return w * (x + a1) / (m + a1 + b1) + (1 - w) * (x + a0) / (m + a0 + b0), w


def gfeat(X, Nn, use):
    pg = np.clip(genus_mean(X, Nn, use), 1e-6, 1 - 1e-6)
    fr = (X[use] > 0).mean(0)
    return np.c_[np.ones(X.shape[1]), logit(pg), fr]


rows = []
rng = np.random.default_rng(1)
for npz in sys.argv[1:-1]:
    d = np.load(npz); Nn = d['Nn']; names = d['names']; nsp = len(Nn)
    G = [d[f'G{i}'] for i in range(nsp)]
    X = np.array([gg.sum(1) for gg in G], float); nC = X.shape[1]
    genus = os.path.basename(npz).split('.')[0]
    for f in range(nsp):
        others = [s for s in range(nsp) if s != f]
        kap, _ = fit_kappa(X, Nn, others)
        pg = np.clip(genus_mean(X, Nn, others), 1e-6, 1 - 1e-6)
        # ancestral-state prior fitted on congeners only, each against features from the remaining ones
        xs, ns, fs = [], [], []
        for s in others:
            rest = [o for o in others if o != s]
            xs.append(X[s]); ns.append(Nn[s]); fs.append(gfeat(X, Nn, rest))
        th_anc = fit_mix(xs, ns, fs, 3)
        Ff = gfeat(X, Nn, others)
        full = X[f] / Nn[f]
        incong = X[others].sum(0) > 0
        strata = {'all': np.ones(nC, bool), 'shared': (full > 0) & incong, 'focal_only': (full > 0) & ~incong,
                  'congener_only': (full == 0) & incong}
        print(genus, names[f], Nn[f], 'kappa', kap, 'anc theta', np.round(th_anc, 2), flush=True)
        for m in MS:
            if Nn[f] < m + 5: continue
            for rep in range(REPS):
                perm = rng.permutation(Nn[f]); tr = perm[:m]; te = perm[m:m + NTEST]
                x = G[f][:, tr].sum(1).astype(float); Y = G[f][:, te].astype(float)
                al, be = fit_beta_eb(x, m)
                th_sp = fit_mix([x], [m], [np.ones((nC, 1))], 1)    # same two-state form, no genus info
                preds = {'sp_eb': (x + al) / (m + al + be), 'shrunk': (x + kap * pg) / (m + kap)}
                preds['mix'] = 0.5 * preds['shrunk'] + 0.5 * preds['sp_eb']
                preds['sp_2state'] = mix_pred(th_sp, x, m, np.ones((nC, 1)))[0]
                preds['anc'], w = mix_pred(th_anc, x, m, Ff)
                for sname, mask in strata.items():
                    if not mask.any(): continue
                    row = dict(genus=genus, focal=names[f], n_focal=int(Nn[f]), m=m, rep=rep, stratum=sname,
                               n_clusters=int(mask.sum()), kappa=kap)
                    for pn, p in preds.items():
                        row[f'bits_{pn}'] = bits(p[mask][:, None], Y[mask]).sum(0).mean()
                    rows.append(row)
pd.DataFrame(rows).to_csv(sys.argv[-1], sep='\t', index=False)
