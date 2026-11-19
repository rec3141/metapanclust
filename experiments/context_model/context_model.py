"""Level-2 'context' surprise: predict each module's (or accessory family's) presence in a genome from the rest
of that genome, leave-one-genome-out over the public genomes, and compare with prevalence-only surprise.

Models (probability that the target is observed present in the held-out genome):
  prev      Beta(1,1) prevalence among training genomes.
  prevc     completeness-aware prevalence: true rate q (Beta(1,1) prior, posterior mean on a grid) with
            observed presence q*c, c = genome completeness; P(present) = q*c_g.
  knnG      Jaccard-weighted neighbours on accessory families (phylogeny proxy): weights J^G, shrunk to
            prev with one pseudo-genome.
  lrmodC    L2 logistic regression on the other modules' presence plus completeness (C = inverse penalty).
  lrfamC    L2 logistic regression on accessory families (shell+cloud, >=2 public genomes) plus completeness.
The target's own families are masked from every feature set.

usage: context_model.py SPECIES OUTDIR [N_FAMILY_TARGETS]
"""
import os, sys, glob, time
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression

sp, out = sys.argv[1], sys.argv[2]
NFAM = int(sys.argv[3]) if len(sys.argv) > 3 else 200
d = f'/data/scratch/metapanclust_data/pangenomes/{sp}/ppanggolin'
os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(1)

R = pd.read_csv(f'{d}/gene_presence_absence.Rtab', sep='\t', index_col=0) > 0
genomes = list(R.columns)
ours = [g for g in genomes if g.startswith('chesterfield')]
pub = [g for g in genomes if g not in ours]
st = pd.read_csv(f'{d}/genomes_statistics.tsv', sep='\t', comment='#').set_index('Genome_name')
comp = (st.Completeness / 100).reindex(genomes).values
part = {}
for pname in ('persistent', 'shell', 'cloud'):
    for l in open(f'{d}/partitions/{pname}.txt'): part[l.strip()] = pname
mig = pd.read_csv(f'{d}/modules_in_genomes.tsv', sep='\t')
mig = mig[mig.completion >= 0.5]
fm = pd.read_csv(f'{d}/functional_modules.tsv', sep='\t')
mods = list(fm.module_id.unique())
M = pd.DataFrame(False, index=mods, columns=genomes)
for m, g in zip(mig.module_id, mig.genome): M.loc[m, g] = True
mod_fams = fm.groupby('module_id').family_id.apply(set).to_dict()

acc = [f for f in R.index if part.get(f) in ('shell', 'cloud')]
npub = R.loc[acc, pub].sum(1)
feat_fams = [f for f in acc if 2 <= npub[f] <= len(pub) - 1]
F = R.loc[feat_fams, genomes].values.T.astype(float)        # genomes x families
Mx = M.values.T.astype(float)                               # genomes x modules
gi = {g: i for i, g in enumerate(genomes)}
P = np.array([gi[g] for g in pub]); O = np.array([gi[g] for g in ours])
fam_idx = {f: i for i, f in enumerate(feat_fams)}
inter_all = F @ F.T; size_all = F.sum(1)

# targets: modules with 2..n-2 public carriers; a random sample of accessory families with 5..n-5 carriers
mk = M[pub].sum(1)
tmods = [m for m in mods if 2 <= mk[m] <= len(pub) - 2]
fk = R.loc[feat_fams, pub].sum(1)
cand = [f for f in feat_fams if 5 <= fk[f] <= len(pub) - 5]
tfams = list(rng.choice(cand, min(NFAM, len(cand)), replace=False)) if cand else []
fam2mod = {f: m for m, fs in mod_fams.items() for f in fs}
print(f'{sp}: {len(pub)} public + {len(ours)} ours; {len(feat_fams)} accessory feature families; '
      f'{len(mods)} modules; targets {len(tmods)} modules, {len(tfams)} families (of {len(cand)})', flush=True)

QG = np.linspace(0.0005, 0.9995, 1000)
def prevc(y, c, cg):
    ll = (y[:, None] * np.log(QG * c[:, None]) + (1 - y[:, None]) * np.log1p(-QG * c[:, None])).sum(0)
    w = np.exp(ll - ll.max()); q = (w * QG).sum() / w.sum()
    return q * cg

GAMMAS = [2, 8, 32]; CS_MOD = [0.1, 1.0]; CS_FAM = [0.01, 0.1, 1.0]

def features(kind, tname):
    """feature matrix (all genomes) with the target's own families / module masked"""
    own = mod_fams.get(tname, {tname}) if kind == 'module' else {tname} | mod_fams.get(fam2mod.get(tname, ''), set())
    keepf = np.array([f not in own for f in feat_fams])
    keepm = np.array([(m != tname) and not (mod_fams[m] & own) for m in mods])
    drop = ~keepf
    inter = inter_all - F[:, drop] @ F[:, drop].T; size = size_all - F[:, drop].sum(1)
    J = inter / np.maximum(size[:, None] + size[None, :] - inter, 1)
    return F[:, keepf], Mx[:, keepm], J

def fit_lr(X, y, C):
    if y.min() == y.max(): return None
    return LogisticRegression(C=C, max_iter=500).fit(X, y)

def predict_all(Xf, Xm, J, y, tr, te):
    """probabilities for genomes te from models trained on genomes tr"""
    out = {}
    ytr = y[tr]; k = ytr.sum(); n = len(tr)
    p0 = (k + 1) / (n + 2)
    out['prev'] = np.full(len(te), p0)
    out['prevc'] = np.array([prevc(ytr, comp[tr], comp[t]) for t in te])
    for G in GAMMAS:
        W = J[np.ix_(te, tr)] ** G
        out[f'knn{G}'] = (W @ ytr + p0) / (W.sum(1) + 1)
    cz = ((comp - comp[tr].mean()) / (comp[tr].std() + 1e-9))[:, None]
    for name, X, Cs in (('lrmod', Xm, CS_MOD), ('lrfam', Xf, CS_FAM)):
        Xc = np.hstack([X, cz])
        for C in Cs:
            mdl = fit_lr(Xc[tr], ytr, C)
            out[f'{name}{C:g}'] = mdl.predict_proba(Xc[te])[:, 1] if mdl else np.full(len(te), p0)
    return out

rows, insample = [], []
t0 = time.time()
targets = [('module', m) for m in tmods] + [('family', f) for f in tfams]
for ti, (kind, t) in enumerate(targets):
    y = (M.loc[t].values if kind == 'module' else R.loc[t, genomes].values).astype(float)
    Xf, Xm, J = features(kind, t)
    for h in P:
        tr = P[P != h]
        pr = predict_all(Xf, Xm, J, y, tr, np.array([h]))
        rows.append(dict(kind=kind, target=t, genome=genomes[h], y=y[h], c=comp[h], **{k: v[0] for k, v in pr.items()}))
    # in-sample fit on all public genomes, for the overfitting gap
    pr = predict_all(Xf, Xm, J, y, P, P)
    for j, h in enumerate(P):
        insample.append(dict(kind=kind, target=t, genome=genomes[h], y=y[h], **{k: v[j] for k, v in pr.items()}))
    if ti % 25 == 0: print(f'  {ti}/{len(targets)} targets, {time.time() - t0:.0f}s', flush=True)
pd.DataFrame(rows).to_csv(f'{out}/loo_predictions.tsv', sep='\t', index=False)
pd.DataFrame(insample).to_csv(f'{out}/insample_predictions.tsv', sep='\t', index=False)

# Chesterfield bins: every module, trained on the public genomes (A) and on public + the other bin (B)
ch = []
for m in mods:
    y = M.loc[m].values.astype(float)
    Xf, Xm, J = features('module', m)
    for o in O:
        others = np.array([x for x in O if x != o])
        for setting, tr in (('public', P), ('public+other_bin', np.concatenate([P, others]))):
            pr = predict_all(Xf, Xm, J, y, tr, np.array([o]))
            ch.append(dict(setting=setting, module=m, genome=genomes[o], y=y[o], c=comp[o], n_pub=int(mk[m]),
                           n_fam=len(mod_fams[m]), **{k: v[0] for k, v in pr.items()}))
pd.DataFrame(ch).to_csv(f'{out}/chesterfield_predictions.tsv', sep='\t', index=False)

# learning curve: see learning_curve.py (reduced, for a shared workstation)
print(f'done {time.time() - t0:.0f}s', flush=True)
