"""Cross-validated logistic model of same-TU from intergenic gap plus per-distance-level conservation.
usage: fit_model.py PANEL_TSV MATRIX_TSV
Asks whether close genomes add information once each level gets its own weight (instead of a union)."""
import csv, sys
import numpy as np

panel = {p['name']: int(p['level']) for p in csv.DictReader(open(sys.argv[1]), delimiter='\t')}
rows = list(csv.DictReader(open(sys.argv[2]), delimiter='\t'))
names = [n for n in rows[0] if n in panel]
y = np.array([int(r['same_tu']) for r in rows], float)
gap = np.array([int(r['gap']) for r in rows])
M = np.array([[int(r[n]) for n in names] for r in rows])          # 0/1/2
C = (M >= 1).astype(float)
lv = np.array([panel[n] for n in names])
edges = [-10**9, 0, 21, 51, 101, 151, 301, 10**9]
G = np.stack([(gap >= lo) & (gap < hi) for lo, hi in zip(edges[:-1], edges[1:])], 1).astype(float)

def lev_any(levels): return np.stack([C[:, lv == L].max(1) for L in levels], 1)
def lev_frac(levels): return np.stack([C[:, lv == L].mean(1) for L in levels], 1)
base5 = C[:, [names.index(n) for n in ('Eco_Sakai', 'Pae', 'Mtu', 'Syn', 'Mja')]].max(1, keepdims=True)
dist = lv >= 3
feats = {
    'gap bins only': G,
    'gap + baseline5 any': np.hstack([G, base5]),
    'gap + any(level>=3)': np.hstack([G, C[:, dist].max(1, keepdims=True)]),
    'gap + count(level>=3)': np.hstack([G, C[:, dist].sum(1, keepdims=True) / dist.sum()]),
    'gap + per-level any, levels 3-9': np.hstack([G, lev_any(range(3, 10))]),
    'gap + per-level any, levels 0-2': np.hstack([G, lev_any(range(0, 3))]),
    'gap + per-level frac, levels 0-2': np.hstack([G, lev_frac(range(0, 3))]),
    'gap + per-level any, levels 0-9': np.hstack([G, lev_any(range(0, 10))]),
    'gap + per-level frac, levels 0-9': np.hstack([G, lev_frac(range(0, 10))]),
    'gap + per-genome indicators': np.hstack([G, C]),
}

def fit(X, y, l2=1e-3, it=25):
    """L2-penalised logistic regression by Newton steps."""
    Xb = np.hstack([X, np.ones((len(X), 1))]); w = np.zeros(Xb.shape[1])
    R = l2 * len(y) * np.diag(np.r_[np.ones(X.shape[1]), 0])
    for _ in range(it):
        p = 1 / (1 + np.exp(-Xb @ w))
        H = Xb.T @ (Xb * (p * (1 - p))[:, None]) + R
        w -= np.linalg.solve(H, Xb.T @ (p - y) + R @ w)
    return w
def pred(w, X): return 1 / (1 + np.exp(-np.hstack([X, np.ones((len(X), 1))]) @ w))

# folds by contiguous blocks of query pairs (rows are in genome order) so neighbouring pairs of one operon stay together
nf = 5; fold = np.arange(len(y)) * nf // len(y)
print(f'{len(y)} pairs; 5-fold CV over contiguous genome blocks; threshold 0.5')
print(f'{"features":36s} acc    prec   rec    logloss')
for lab, X in feats.items():
    p = np.zeros(len(y))
    for f in range(nf):
        w = fit(X[fold != f], y[fold != f]); p[fold == f] = pred(w, X[fold == f])
    c = p >= 0.5; tp = (c & (y == 1)).sum()
    ll = -np.mean(y * np.log(p + 1e-9) + (1 - y) * np.log(1 - p + 1e-9))
    print(f'{lab:36s} {np.mean(c == y):.3f}  {tp / c.sum():.3f}  {tp / y.sum():.3f}  {ll:.3f}')
w = fit(feats['gap + per-level any, levels 0-9'], y)
print('full-data weights, per-level any (log-odds added when conserved at that level):')
print('  ' + '  '.join(f'L{L}:{w[G.shape[1] + i]:+.2f}' for i, L in enumerate(range(10))))
