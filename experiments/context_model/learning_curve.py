"""Reduced learning curve for the context model (module targets only, three subsamples per size, four models),
light enough for a shared workstation (~5 CPU-min per species). Reuses the data loading and model code of
context_model.py (everything above its main loop).
usage: learning_curve.py SPECIES OUTDIR"""
import os, sys
os.environ['OMP_NUM_THREADS'] = '1'; os.environ['OPENBLAS_NUM_THREADS'] = '1'
here = os.path.dirname(os.path.abspath(__file__))
src = open(f'{here}/context_model.py').read().split('rows, insample = [], []')[0]
sys.argv = [sys.argv[0], sys.argv[1], sys.argv[2], '0']
exec(src)
GAMMAS[:] = [8]; CS_MOD[:] = [1.0]; CS_FAM[:] = [0.1]
lc = []
sizes = [s for s in (10, 25, 40) if s < len(pub) - 1] + [len(pub) - 1]
for t in tmods:
    y = M.loc[t].values.astype(float)
    Xf, Xm, J = features('module', t)
    for n in sizes:
        for rep in range(3 if n < len(pub) - 1 else 1):
          for h in P:
            pool = P[P != h]; tr = rng.choice(pool, n, replace=False) if n < len(pool) else pool
            pr = predict_all(Xf, Xm, J, y, tr, np.array([h]))
            lc.append(dict(kind='module', target=t, n=n, rep=rep, genome=genomes[h], y=y[h], **{k: v[0] for k, v in pr.items()}))
pd.DataFrame(lc).to_csv(f'{out}/learning_curve.tsv', sep='\t', index=False)
print('learning curve done', flush=True)
