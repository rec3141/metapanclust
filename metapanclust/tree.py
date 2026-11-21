"""Level 1: presence of a unit given relatedness, by a two-state gain/loss chain on a gene-content tree.

tree_vs_context (200 species): within a species a UPGMA tree on accessory-family Jaccard (targets masked) is
the best relatedness predictor; xsp_tree_vs_content: between species the GTDB tree is marginally better.
This module does the within-species case. The (rate, stationary frequency) of each unit is averaged over a
grid, weighted by the marginal likelihood of the reference tips (Felsenstein pruning); the query is an
uncertain tip and its posterior comes from an up-down pass.
"""
import numpy as np

MUH = np.array([0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000])   # rate x root height; 1000 ~ no tree signal
PIS = np.array([0.003, 0.01, 0.03, 0.1, 0.25, 0.5, 0.75, 0.9, 0.97, 0.99])   # stationary frequency
GMU = np.repeat(MUH, len(PIS)); GPI = np.tile(PIS, len(MUH))
EPS = 0.002      # P(observed | absent) for reference tips: annotation / clustering noise


def jaccard_dist(F):
    F = F.astype(np.float32)
    inter = F @ F.T; sz = F.sum(1)
    D = 1 - inter / np.maximum(sz[:, None] + sz[None, :] - inter, 1)
    np.fill_diagonal(D, 0)
    return D


def upgma(D):
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import squareform
    N = len(D)
    Z = linkage(squareform(np.maximum(D, 0), checks=False), method='average')
    ch = Z[:, :2].astype(int)
    h = np.concatenate([np.zeros(N), Z[:, 2] / 2])
    h[N:] = np.maximum.accumulate(h[N:])
    bl = np.zeros(2 * N - 1)
    for i in range(N - 1):
        for c in ch[i]: bl[c] = h[N + i] - h[c]
    H = max(h[-1], 1e-9)
    bl[:N] = np.maximum(bl[:N], 0.01 * H)
    return ch, bl, H


def tree_posterior(tree, Y, observed, comp, chunk=200):
    """Y: N x T observed states (0/1) for tips; observed: bool N (False = unknown tip, e.g. the query);
    comp: N detection probabilities. Returns N x T posterior P(unit present) (for observed tips this is the
    leave-one-out posterior: outside message only)."""
    N, T = Y.shape
    out = np.zeros((N, T))
    for a in range(0, T, chunk):
        out[:, a:a + chunk] = _posterior(tree, Y[:, a:a + chunk], observed, comp)
    return out


def _posterior(tree, Y, train, comp):
    ch, bl, H = tree
    N, T = Y.shape; G = len(GMU); nn = 2 * N - 1
    c = np.clip(comp, 0.3, 0.99)[:, None]
    L = np.ones((N, T, 2))
    L[train, :, 0] = np.where(Y[train] > 0, EPS, 1 - EPS)
    L[train, :, 1] = np.where(Y[train] > 0, c[train], 1 - c[train])
    mu = GMU[:, None] / H; pi = GPI[:, None]
    up = np.empty((nn, G, T, 2)); msg = np.empty((nn, G, T, 2)); lsc = np.zeros((G, T))
    up[:N] = L[:, None]

    def message(v):
        e = np.exp(-mu * bl[v])
        u0, u1 = up[v, ..., 0], up[v, ..., 1]
        msg[v, ..., 0] = (1 - pi * (1 - e)) * u0 + pi * (1 - e) * u1
        msg[v, ..., 1] = (1 - pi) * (1 - e) * u0 + (pi + (1 - pi) * e) * u1
    for v in range(N): message(v)
    for i in range(N - 1):
        v = N + i; l, r = ch[i]
        u = msg[l] * msg[r]
        s = u.max(-1, keepdims=True); s[s <= 0] = 1
        up[v] = u / s; lsc += np.log(s[..., 0])
        if v < nn - 1: message(v)
    root = nn - 1
    ll = lsc + np.log(np.maximum((1 - pi) * up[root, ..., 0] + pi * up[root, ..., 1], 1e-300))
    w = np.exp(ll - ll.max(0, keepdims=True)); w /= w.sum(0, keepdims=True)
    down = np.empty((nn, G, T, 2))
    down[root, ..., 0] = np.broadcast_to(1 - pi, (G, T)); down[root, ..., 1] = np.broadcast_to(pi, (G, T))
    for i in range(N - 2, -1, -1):
        v = N + i; l, r = ch[i]
        for cnode, sib in ((l, r), (r, l)):
            x = down[v] * msg[sib]
            e = np.exp(-mu * bl[cnode])
            d0 = x[..., 0] * (1 - pi * (1 - e)) + x[..., 1] * (1 - pi) * (1 - e)
            d1 = x[..., 0] * pi * (1 - e) + x[..., 1] * (pi + (1 - pi) * e)
            s = np.maximum(d0 + d1, 1e-300)
            down[cnode, ..., 0] = d0 / s; down[cnode, ..., 1] = d1 / s
    return (down[:N, ..., 1] * w[None]).sum(1)


def level1(ref, q_fam, q_comp, units, folds=5, seed=0):
    """P(present in the query) for the given units (indices) from the gene-content tree on ref + query.
    q_fam: bool vector of families observed in the query. The distance is Jaccard on accessory families.
    Units are split into `folds` groups; each group's posterior uses a tree whose query placement ignores the
    families of that group's units, so a unit never informs its own placement."""
    rng = np.random.default_rng(seed)
    acc = ref.part != 'P'
    grp = rng.permutation(np.arange(len(units)) % folds)
    N = len(ref.genomes)
    p = np.zeros(len(units))
    obs = np.r_[np.ones(N, bool), False]
    comp = np.r_[ref.comp, q_comp]
    for g in range(folds):
        sel = np.flatnonzero(grp == g)
        if not len(sel): continue
        fk = acc.copy()
        for u in units[sel]: fk[ref.unit_members[u]] = False
        F = np.vstack([ref.R[fk].T, q_fam[fk][None]])
        tree = upgma(jaccard_dist(F))
        Y = np.vstack([ref.U[units[sel]].T.astype(float), np.zeros((1, len(sel)))])
        p[sel] = tree_posterior(tree, Y, obs, comp)[-1]
    return p
