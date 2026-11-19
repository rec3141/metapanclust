"""Test of cluster_model.md P2: per-genome log-likelihood-ratio evidence for same-TU adjacent pairs.
usage: evidence.py TAG PANEL_TSV MATRIX_TSV RUN_DIR BAKTA_TSV OUT_DIR

Score per pair: sum_k w_k [c_k (1-phi) x_k + (1-c_k) ln((1-e^{-phi x_k})/(1-e^{-x_k}))], phi = 0.53.
x_k = -ln P(conserved | different-TU pair) in genome k, estimated in the training folds only.
Variant 'pos' (conserved-only): breaks carry no evidence.
Variant 'pres': genomes where either gene has no homologue (gene absent from every cluster of that
run) give no evidence, and x_k is conditioned on both genes present.
Classifier: logit = gap-bin prior log-odds (training folds) + score, call same-TU if logit > 0.
5-fold CV over contiguous blocks of pairs."""
import csv, math, os, random, sys
import numpy as np

tag, panel_tsv, mat_tsv, run_dir, bakta_tsv, out_dir = sys.argv[1:7]
PHI = 0.53
panel = list(csv.DictReader(open(panel_tsv), delimiter='\t'))
level = {p['name']: int(p['level']) for p in panel}

rows = [l.rstrip('\n').split('\t') for l in open(mat_tsv)]
names = rows[0][4:]
pairs = [(r[0], r[1]) for r in rows[1:]]
gap = np.array([int(r[2]) for r in rows[1:]])
y = np.array([int(r[3]) for r in rows[1:]])
C = np.array([[int(v) >= 1 for v in r[4:]] for r in rows[1:]])
N, K = C.shape
lv = np.array([level[n] for n in names])

# presence: query gene appears in any cluster of the run against genome k (homologue found)
pres_f = f'{out_dir}/{tag}_presence.npy'
if os.path.exists(pres_f):
    Pr = np.load(pres_f)
else:
    loc = {}
    for l in open(bakta_tsv):
        if l.startswith('#'): continue
        f = l.rstrip('\n').split('\t')
        if len(f) > 5 and f[1] == 'cds': loc[(f[0], int(f[2]), int(f[3]))] = f[5]
    Pr = np.zeros((N, K), bool)
    for k, n in enumerate(names):
        seen = set()
        for l in open(f'{run_dir}/{n}/clusters.tsv'):
            if l[0] != '>': continue
            c, i, a, b = l[1:].split('\t')[0].rsplit('_', 3)
            a, b = int(a), int(b)
            g = loc.get((c, min(a, b), max(a, b)))
            if g: seen.add(g)
        Pr[:, k] = [a in seen and b in seen for a, b in pairs]
    np.save(pres_f, Pr)
Pr = Pr | C

folds = np.minimum(np.arange(N) * 5 // N, 4)
EDGES = [-10**9, -20, -1, 10, 20, 30, 40, 50, 60, 80, 100, 130, 160, 200, 250, 300, 400, 600, 10**9]
gb = np.digitize(gap, EDGES[1:-1], right=True)

def gap_prior(tr):
    lo = np.zeros(len(EDGES) - 1)
    for b in range(len(lo)):
        m = tr & (gb == b)
        lo[b] = math.log((y[m].sum() + 1) / ((1 - y[m]).sum() + 1))
    return lo[gb]

def llr_terms(tr, cols, pres):
    """per-genome conserved / broken LLR, x from training folds' different-TU pairs"""
    neg = tr & (y == 0)
    if pres is True:
        den = (Pr[neg][:, cols]).sum(0)
        pn = (C[neg][:, cols]).sum(0) / np.maximum(den, 1)
    else:
        pn = C[neg][:, cols].mean(0)
    pn = np.clip(pn, 1e-3, 1 - 1e-3)
    x = -np.log(pn)
    lc = (1 - PHI) * x
    lb = np.log((1 - np.exp(-PHI * x)) / (1 - np.exp(-x)))
    return lc, lb, x

def score(tr, cols, w, pres):
    lc, lb, _ = llr_terms(tr, cols, pres)
    c = C[:, cols]; m = Pr[:, cols] if pres is True else np.ones_like(c)
    if pres == 'pos': return (c * lc * w).sum(1)
    return ((c * lc + (m & ~c) * lb) * w).sum(1)

def acc(pred, mask=None):
    return (pred == y).mean() if mask is None else (pred[mask] == y[mask]).mean()

def prf(pred):
    tp = (pred & (y == 1)).sum()
    return tp / max(pred.sum(), 1), tp / y.sum()

def cv(cols, w, pres, mode='prior', slope=1.0):
    """CV accuracy. mode prior: gap prior + slope*score > 0 (slope fixed). logit: logistic on [gap bins, score]."""
    pred = np.zeros(N, bool)
    for f in range(5):
        tr = folds != f; te = ~tr
        s = score(tr, cols, w, pres)
        if mode == 'prior':
            lo = gap_prior(tr) + slope * s
            pred[te] = lo[te] > 0
        else:
            from_bins = np.eye(len(EDGES) - 1)[gb]
            X = np.column_stack([from_bins, s])
            beta = logit_fit(X[tr], y[tr]); SLOPES.append(beta[-1])
            pred[te] = X[te] @ beta > 0
    return pred

def logit_fit(X, t, l2=1e-2, it=50):
    b = np.zeros(X.shape[1])
    for _ in range(it):
        p = 1 / (1 + np.exp(-X @ b))
        g = X.T @ (p - t) + l2 * b
        H = X.T @ (X * (p * (1 - p))[:, None]) + l2 * np.eye(len(b))
        b -= np.linalg.solve(H, g)
    return b

def union_rule(cols, T1, T2):
    any_c = C[:, cols].any(1)
    return (gap <= T1) | (any_c & (gap <= T2))

out = []
SLOPES = []
def rep(label, pred):
    p, r = prf(pred)
    sl = f'\tfitted LLR slope {np.mean(SLOPES):.2f}' if SLOPES else ''
    SLOPES.clear()
    out.append(f'{label}\t{acc(pred):.3f}\t{p:.3f}\t{r:.3f}{sl}')

T1, T2 = (50, 300) if tag == 'eco' else (100, 250)
allc = np.arange(K)
out.append(f'## {tag}: {N} pairs, {K} genomes, base rate {y.mean():.3f}, union-rule thresholds T1={T1} T2={T2}')
out.append('# x and P(c|neutral) per level (all pairs, for reference): level, n, mean x, mean x|present, frac pairs with both present')
for L in sorted(set(lv)):
    cols = np.where(lv == L)[0]
    _, _, x = llr_terms(np.ones(N, bool), cols, False)
    _, _, xp = llr_terms(np.ones(N, bool), cols, True)
    out.append(f'#  L{L}\t{len(cols)}\t{x.mean():.2f}\t{xp.mean():.2f}\t{Pr[:, cols].mean():.3f}')
out.append('label\tacc\tprec\trec')
rep('gap prior only (CV)', cv(allc, np.zeros(K), False))
rep(f'union rule, all {K}', union_rule(allc, T1, T2))
best_lo = 6 if tag == 'eco' else 7
rep(f'union rule, >= level {best_lo}', union_rule(np.where(lv >= best_lo)[0], T1, T2))

wlev = np.array([1 / (lv == L).sum() for L in lv])           # each level = one effective genome
for pres in (False, True, 'pos'):
    pv = {False: 'theory', True: 'pres', 'pos': 'conserved-only'}[pres]
    rep(f'[{pv}] prior+LLR all genomes, w=1, slope 1', cv(allc, np.ones(K), pres))
    rep(f'[{pv}] prior+LLR all genomes, w=1/n_level, slope 1', cv(allc, wlev, pres))
    rep(f'[{pv}] logistic[gap bins, LLR] all genomes, w=1', cv(allc, np.ones(K), pres, 'logit'))
    rep(f'[{pv}] logistic[gap bins, LLR] all genomes, w=1/n_level', cv(allc, wlev, pres, 'logit'))
    for lo_ in range(0, 9):
        cols = np.where(lv >= lo_)[0]
        rep(f'[{pv}] prior+LLR levels >= {lo_} (n={len(cols)}), w=1', cv(cols, np.ones(len(cols)), pres))
    for hi in (3, 4):
        cols = np.where(lv <= hi)[0]
        rep(f'[{pv}] prior+LLR levels <= {hi} only (n={len(cols)}), w=1', cv(cols, np.ones(len(cols)), pres))

# where does the error come from: TU pairs vs non-TU pairs, near-gap
pred = cv(allc, np.ones(K), True)
fp = (pred & (y == 0)); fn = (~pred & (y == 1))
out.append(f'# best-guess model [pres, w=1, all]: FP {fp.sum()} (gap<=50: {(fp & (gap <= 50)).sum()}), FN {fn.sum()} (gap>300: {(fn & (gap > 300)).sum()})')

# spread vs single-clade panels of equal size (pres variant, w=1, prior+LLR); 40 random draws
random.seed(1)
if tag == 'eco':
    clades = {'L5 Gamma': lv == 5, 'L6 Beta': lv == 6, 'L7 Alpha': lv == 7,
              'Bacillota': np.isin(names, ['Bsu', 'Sau', 'Lmo', 'Cac'])}
    spread_levels = [5, 6, 7, 8]
else:
    clades = {'L3 Bacillales other orders': lv == 3, 'L7 Pseudomonadota': lv == 7, 'L8 other phyla': lv == 8,
              'L5 other Bacillota': lv == 5}
    spread_levels = [4, 5, 6, 7, 8]
out.append('# panel size n: single-clade vs spread (one genome per level from ' + str(spread_levels) + '); mean over draws; acc = prior+LLR [conserved-only], union acc = union rule')
out.append('panel\tn\tdraws\tacc_LLR\tacc_union\tcons_recall_any')
for n in (3, 4):
    for cname, m in list(clades.items()) + [('spread', None)]:
        if m is not None and m.sum() < n: continue
        res = []
        for d in range(40):
            if m is None:
                Ls = random.sample(spread_levels, n)
                cols = np.array([random.choice(list(np.where(lv == L)[0])) for L in Ls])
            else:
                cols = np.array(random.sample(list(np.where(m)[0]), n))
            a1 = acc(cv(cols, np.ones(n), 'pos')); a2 = acc(union_rule(cols, T1, T2))
            rc = (C[:, cols].any(1) & (y == 1)).sum() / y.sum()
            res.append((a1, a2, rc))
        r = np.array(res).mean(0)
        out.append(f'{cname}\t{n}\t{len(res)}\t{r[0]:.3f}\t{r[1]:.3f}\t{r[2]:.3f}')

txt = '\n'.join(out)
print(txt)
open(f'{out_dir}/{tag}_results.txt', 'w').write(txt + '\n')
