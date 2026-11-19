"""Combined rule 'gap <= T1 OR (conserved AND gap <= T2)' with thresholds re-tuned for the query, per panel.
The E. coli thresholds (50/300) do not transfer: B. subtilis same-operon pairs have longer spacers.
usage: score_tuned.py PANEL_TSV MATRIX_TSV OUT_TSV   (MATRIX_TSV from score_panels.py)
cv_acc: 5-fold CV over contiguous genome blocks, thresholds chosen on the training folds.
far_*: among pairs with gap > T1 of the gap-only optimum (where distance alone says 'different TU')."""
import csv, random, sys
import numpy as np
panel = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
level = {p['name']: int(p['level']) for p in panel}
rows = list(csv.DictReader(open(sys.argv[2]), delimiter='\t'))
names = [n for n in rows[0] if n in level]
y = np.array([int(r['same_tu']) for r in rows], bool); gap = np.array([int(r['gap']) for r in rows])
C = np.array([[int(r[n]) >= 1 for n in names] for r in rows])
lv = np.array([level[n] for n in names])
T1s = list(range(0, 201, 10)); T2s = [None, 100, 150, 200, 250, 300, 400, 500, 700, 1000]
nf = 5; fold = np.arange(len(y)) * nf // len(y)
def acc(c, idx, t1, t2):
    call = (gap[idx] <= t1) | (c[idx] & (gap[idx] <= t2 if t2 else False))
    return np.mean(call == y[idx]), call
def best(c, idx):
    return max(((acc(c, idx, t1, t2)[0], t1, t2) for t1 in T1s for t2 in T2s if t2 is None or t2 > t1), key=lambda x: x[0])
allidx = np.arange(len(y))
_, g1, _ = best(np.zeros(len(y), bool), allidx)
def score(mask):
    c = C[:, mask].any(1) if mask.any() else np.zeros(len(y), bool)
    a, t1, t2 = best(c, allidx); _, call = acc(c, allidx, t1, t2)
    cv = np.zeros(len(y), bool)
    for f in range(nf):
        _, a1, a2 = best(c, allidx[fold != f]); cv[fold == f] = acc(c, allidx[fold == f], a1, a2)[1]
    far = gap > g1
    return dict(n=int(mask.sum()), cons_prec=(c & y).sum() / max(c.sum(), 1), cons_rec=(c & y).sum() / y.sum(),
                far_prec=(c & y & far).sum() / max((c & far).sum(), 1), far_rec=(c & y & far).sum() / (y & far).sum(),
                T1=t1, T2=t2 or '-', acc=a, cv_acc=np.mean(cv == y), prec=(call & y).sum() / call.sum(), rec=(call & y).sum() / y.sum())
panels = {'gap only': lv < -1}
for L in range(10): panels[f'level{L}'] = lv == L
for L in range(1, 10): panels[f'ge_level{L}'] = lv >= L
panels['levels 3-9 minus Pseudomonadota (7)'] = (lv >= 3) & (lv != 7)
panels['levels 4-9 (outside NCBI Bacillales)'] = lv >= 4
panels['levels 5-9 (outside class Bacilli)'] = lv >= 5
panels['all'] = lv >= 0
hdr = ['panel', 'n', 'cons_prec', 'cons_rec', 'far_prec', 'far_rec', 'T1', 'T2', 'acc', 'cv_acc', 'prec', 'rec']
out = open(sys.argv[3], 'w'); out.write('\t'.join(hdr) + '\n')
print(f'gap-only optimum T1={g1}; far = gap > {g1}: {(y & (gap > g1)).sum()} same-TU, {(~y & (gap > g1)).sum()} different-TU')
for k, m in panels.items():
    s = score(m); out.write('\t'.join([k] + [f'{s[h]:.3f}' if isinstance(s[h], float) else str(s[h]) for h in hdr[1:]]) + '\n')
# recall vs number of genomes, random draws (30 per size) from the pool outside GTDB o__Bacillales (levels >= 3) and outside NCBI Bacillales (>= 4)
rng = random.Random(1)
for pool_name, pool in (('L3-9', np.where(lv >= 3)[0]), ('L4-9', np.where(lv >= 4)[0]), ('L6-9', np.where(lv >= 6)[0])):
    for size in (1, 2, 5, 10, 20, 30, 40, 60):
        if size > len(pool): continue
        res = []
        for _ in range(30):
            m = np.zeros(len(names), bool); m[rng.sample(list(pool), size)] = True
            c = C[:, m].any(1); far = gap > g1
            res.append(((c & y).sum() / max(c.sum(), 1), (c & y).sum() / y.sum(), (c & y & far).sum() / max((c & far).sum(), 1),
                        (c & y & far).sum() / (y & far).sum(), best(c, allidx)[0]))
        r = np.mean(res, 0)
        out.write(f'random {pool_name} n={size}\t{size}\t{r[0]:.3f}\t{r[1]:.3f}\t{r[2]:.3f}\t{r[3]:.3f}\t-\t-\t{r[4]:.3f}\t-\t-\t-\n')
out.close(); print(open(sys.argv[3]).read())
