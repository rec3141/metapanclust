"""One task of the PanGBank data-size experiment.

  worker.py species ID H5DIR OUTDIR [--threads T] [--keep]
      Experiment A (level 2). Exports the presence/absence matrix, partitions, modules and genome statistics
      of pangenome ID from H5DIR/ID.h5, then for genome subsamples of size N (15, 25, 50, 100, 200, 400, 800;
      3 random subsamples each, 1 for 800; plus all genomes if <= 400) runs 10-fold cross-validation of the
      context model of experiments/context_model/context_model.py (L2 logistic regression on accessory-family presence +
      completeness; target's own module masked; C in 0.03..1) against Beta(1,1) prevalence and against
      Jaccard-weighted nearest neighbours (clade proxy), all trained on the same genomes.
      Targets per subsample: up to 60 modules with 2..N-2 carriers and up to 60 accessory families with
      max(2, N/20)..N-max(2, N/20) carriers.
      Writes OUTDIR/A/ID.summary.tsv (one row per N, rep, target kind) and OUTDIR/A/ID.targets.tsv.gz.

  worker.py genus GENUS ID,ID,... H5DIR OUTDIR [--threads T] [--keep]
      Experiment B (level 1). Exports each congeneric species' presence/absence matrix (subsampled to at most
      200 genomes) and family representative proteins, clusters all representatives with MMseqs2
      (50% identity, 80% coverage of both) into genus-level families, and for each focal species and training
      size m in (5, 10, 20) compares held-out log loss of
        sp      species-only Beta(1,1) prevalence from m training genomes,
        sp_eb   species-only with one Beta(a, b) for all clusters fitted on the m training genomes (empirical
                Bayes; the fair baseline, since most genus clusters are absent from any one species),
        genus   leave-focal-out genus prevalence (species capped at 20 genomes each),
        shrunk  (x + kappa p_genus) / (m + kappa), kappa fitted by empirical Bayes (beta-binomial marginal
                likelihood) on the other species of the genus only,
        mix     0.5 shrunk + 0.5 sp_eb (fixed, untuned),
        oracle  shrunk with the kappa from the grid that is best on the test genomes (upper bound).
      Strata: all clusters; by the focal species' full-data prevalence (core >= 0.95, accessory, rare_or_absent
      <= 0.05); and by where the cluster occurs (shared: focal and >= 1 congener; focal_only; congener_only).
      Writes OUTDIR/B/GENUS.summary.tsv and OUTDIR/B/GENUS.kappa.tsv.

Work files go to OUTDIR/work/<task> and are deleted at the end unless --keep.
Bits are -log2 P(observed state), clipped at P >= 1e-6.
"""
import os, sys, time, subprocess, shutil, argparse, gzip
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np, pandas as pd
from scipy.special import gammaln
from sklearn.linear_model import LogisticRegression

ap = argparse.ArgumentParser()
ap.add_argument('mode', choices=['species', 'genus'])
ap.add_argument('args', nargs='+')
ap.add_argument('--threads', type=int, default=1)
ap.add_argument('--keep', action='store_true')
ap.add_argument('--seed', type=int, default=1)
a = ap.parse_args()
T0 = time.time()
LOG2 = np.log(2)


def log(msg): print(f'[{time.time() - T0:7.0f}s] {msg}', flush=True)


def bits(p, y):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return -(y * np.log2(p) + (1 - y) * np.log2(1 - p))


def sh(cmd, logf):
    with open(logf, 'a') as lf:
        r = subprocess.run(cmd, shell=True, stdout=lf, stderr=subprocess.STDOUT)
    if r.returncode: raise RuntimeError(f'failed ({r.returncode}): {cmd}; see {logf}')


def export(h5, wd, proteins=False):
    """PPanGGOLiN flat tables (and family representative proteins) from a PanGBank .h5"""
    os.makedirs(wd, exist_ok=True)
    if not os.path.exists(f'{wd}/flat/gene_presence_absence.Rtab'):
        sh(f'ppanggolin write_pangenome -p {h5} --Rtab --modules --stats --partitions -o {wd}/flat -f '
           f'--disable_prog_bar -c 1', f'{wd}/ppanggolin.log')
    if proteins and not os.path.exists(f'{wd}/fa/all_protein_families.faa'):
        sh(f'ppanggolin fasta -p {h5} --prot_families all -o {wd}/fa -f --disable_prog_bar',
           f'{wd}/ppanggolin.log')
    d = f'{wd}/flat'
    hdr = open(f'{d}/gene_presence_absence.Rtab').readline().rstrip('\n').split('\t')
    R = pd.read_csv(f'{d}/gene_presence_absence.Rtab', sep='\t', index_col=0,
                    dtype={h: np.uint8 for h in hdr[1:]}) > 0
    st = pd.read_csv(f'{d}/genomes_statistics.tsv', sep='\t', comment='#').set_index('Genome_name')
    comp = (st.Completeness.reindex(R.columns).fillna(st.Completeness.median()) / 100).values
    part = {}
    for pname in ('persistent', 'shell', 'cloud'):
        f = f'{d}/partitions/{pname}.txt'
        if os.path.exists(f):
            for l in open(f): part[l.strip()] = pname
    mods = {}
    if os.path.exists(f'{d}/functional_modules.tsv'):
        fm = pd.read_csv(f'{d}/functional_modules.tsv', sep='\t')
        mods = fm.groupby('module_id').family_id.apply(set).to_dict()
    M = pd.DataFrame(False, index=list(mods), columns=R.columns)
    if os.path.exists(f'{d}/modules_in_genomes.tsv') and mods:
        mig = pd.read_csv(f'{d}/modules_in_genomes.tsv', sep='\t')
        mig = mig[(mig.completion >= 0.5) & mig.genome.isin(R.columns)]
        M.values[M.index.get_indexer(mig.module_id), M.columns.get_indexer(mig.genome)] = True
    return R, comp, part, M, mods


# ------------------------------------------------------------------------------------------------- A
CS = [0.03, 0.1, 0.3, 1.0]
GAMMAS = [2, 8]          # Jaccard-weighted neighbours on accessory families: a phylogeny / clade proxy
MAX_TGT = 60


def run_A(pid, h5dir, out):
    wd = f'{out}/work/A_{pid}'
    R, comp, part, M, mod_fams = export(f'{h5dir}/{pid}.h5', wd)
    n = R.shape[1]
    log(f'{pid}: {n} genomes, {R.shape[0]} families, {len(mod_fams)} modules')
    fams = np.array(R.index); Rv = R.values; Mv = M.values; mods = list(M.index)
    acc_all = np.array([part.get(f) in ('shell', 'cloud') for f in fams])
    fam_pos = {f: i for i, f in enumerate(fams)}
    fam2mod = {f: m for m, fs in mod_fams.items() for f in fs}
    rng = np.random.default_rng(a.seed)
    sizes = [s for s in (15, 25, 50, 100, 200, 400, 800) if s < n]
    if n <= 400: sizes.append(n)
    summ, trows = [], []
    for N in sizes:
        reps = 1 if N == n or N >= 800 else 3
        for rep in range(reps):
            t1 = time.time()
            S = np.sort(rng.choice(n, N, replace=False)) if N < n else np.arange(n)
            RS = Rv[:, S]; k = RS.sum(1)
            feat = acc_all & (k >= 2) & (k <= N - 1)
            fidx = np.where(feat)[0]
            F = RS[fidx].T.astype(np.float32)
            c = comp[S]; cz = ((c - c.mean()) / (c.std() + 1e-9)).astype(np.float32)[:, None]
            MS = Mv[:, S]; mk = MS.sum(1)
            tm = [i for i in range(len(mods)) if 2 <= mk[i] <= N - 2]
            lo = max(2, N // 20)
            tf = [i for i in fidx if lo <= k[i] <= N - lo]
            if len(tm) > MAX_TGT: tm = list(rng.choice(tm, MAX_TGT, replace=False))
            if len(tf) > MAX_TGT: tf = list(rng.choice(tf, MAX_TGT, replace=False))
            K = min(N, 10)
            fold = rng.permutation(np.arange(N) % K)
            fpos_in_feat = {f: j for j, f in enumerate(fidx)}
            for kind, tl in (('module', tm), ('family', tf)):
                for ti in tl:
                    if kind == 'module':
                        y = MS[ti].astype(float); own = mod_fams[mods[ti]]; tname = mods[ti]
                    else:
                        y = RS[ti].astype(float); tname = fams[ti]
                        own = {tname} | mod_fams.get(fam2mod.get(tname, ''), set())
                    drop = [fpos_in_feat[fam_pos[f]] for f in own if f in fam_pos and fam_pos[f] in fpos_in_feat]
                    keep = np.ones(len(fidx), bool); keep[drop] = False
                    Fk = F[:, keep]; X = np.hstack([Fk, cz])
                    inter = Fk @ Fk.T; sz = Fk.sum(1)
                    J = inter / np.maximum(sz[:, None] + sz[None, :] - inter, 1)
                    pr = {'prev': np.zeros(N), **{f'knn{g}': np.zeros(N) for g in GAMMAS},
                          **{f'lr{C:g}': np.zeros(N) for C in CS}}
                    for fo in range(K):
                        te = fold == fo; tr = ~te
                        ytr = y[tr]; p0 = (ytr.sum() + 1) / (tr.sum() + 2)
                        pr['prev'][te] = p0
                        for g in GAMMAS:
                            W = J[np.ix_(te, tr)] ** g
                            pr[f'knn{g}'][te] = (W @ ytr + p0) / (W.sum(1) + 1)
                        for C in CS:
                            if ytr.min() == ytr.max(): pr[f'lr{C:g}'][te] = p0; continue
                            mdl = LogisticRegression(C=C, max_iter=300).fit(X[tr], ytr)
                            pr[f'lr{C:g}'][te] = mdl.predict_proba(X[te])[:, 1]
                    pr['mix0.1'] = 0.9 * pr['lr0.1'] + 0.1 * pr['prev']
                    b = {m: bits(p, y).mean() for m, p in pr.items()}
                    trows.append(dict(id=pid, N=N, rep=rep, kind=kind, target=tname, carriers=int(y.sum()),
                                      n_feat=int(keep.sum()), **{f'bits_{m}': round(v, 5) for m, v in b.items()}))
            dt = time.time() - t1
            T = pd.DataFrame([r for r in trows if r['N'] == N and r['rep'] == rep])
            for kind in ('module', 'family'):
                t = T[T.kind == kind] if len(T) else T
                row = dict(id=pid, n_total=n, N=N, rep=rep, kind=kind, folds=K, n_targets=len(t),
                           n_feat=len(fidx), sec=round(dt, 1))
                if len(t):
                    for m in ['prev'] + [f'knn{g}' for g in GAMMAS] + [f'lr{C:g}' for C in CS] + ['mix0.1']:
                        row[f'bits_{m}'] = round(t[f'bits_{m}'].mean(), 5)
                    row['frac_better'] = round((t['bits_lr0.1'] < t.bits_prev - 0.05).mean(), 3)
                    row['frac_worse'] = round((t['bits_lr0.1'] > t.bits_prev + 0.05).mean(), 3)
                summ.append(row)
            log(f'N={N} rep={rep}: {len(tm)} modules, {len(tf)} families, {len(fidx)} features, {dt:.0f}s')
            # rewritten after every subsample, so a task killed at the time limit keeps its smaller sizes
            os.makedirs(f'{out}/A', exist_ok=True)
            pd.DataFrame(trows).to_csv(f'{out}/A/{pid}.targets.tsv.gz', sep='\t', index=False)
            pd.DataFrame(summ).to_csv(f'{out}/A/{pid}.summary.tsv', sep='\t', index=False)
    open(f'{out}/A/{pid}.done', 'w').close()
    if not a.keep: shutil.rmtree(wd, ignore_errors=True)


# ------------------------------------------------------------------------------------------------- B
NB_MAX = 200
CAP = 20                     # max effective genomes per species in the genus mean (model_spec 3.4 A)
A0 = B0 = 0.5                # root prior of the genus mean
KGRID = np.array([0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000])
MS_TRAIN = (5, 10, 20)
REPS = 5
NTEST = 50


def betabin_ll(x, n, al, be):
    return (gammaln(n + 1) - gammaln(x + 1) - gammaln(n - x + 1)
            + gammaln(x + al) + gammaln(n - x + be) - gammaln(n + al + be)
            - gammaln(al) - gammaln(be) + gammaln(al + be))


def genus_mean(X, Nn, use):
    """capped pooled prevalence over species in `use` (X: species x clusters counts, Nn: genomes)"""
    w = np.minimum(Nn[use], CAP).astype(float)
    frac = X[use] / Nn[use][:, None]
    return ((w[:, None] * frac).sum(0) + A0) / (w.sum() + A0 + B0)


def fit_kappa(X, Nn, others):
    """empirical Bayes kappa from the species in `others`, each scored against the leave-it-out genus mean
    of the remaining ones"""
    ll = np.zeros(len(KGRID))
    for s in others:
        rest = [o for o in others if o != s]
        p = np.clip(genus_mean(X, Nn, rest), 1e-6, 1 - 1e-6)
        for j, kap in enumerate(KGRID):
            ll[j] += betabin_ll(X[s], Nn[s], kap * p, kap * (1 - p)).sum()
    return KGRID[int(np.argmax(ll))], ll


EBG = np.exp(np.linspace(np.log(0.005), np.log(20), 30))


def fit_beta_eb(x, m):
    """species-only empirical Bayes: one Beta(a, b) shared by all clusters, fitted on the m training genomes"""
    best, ab = -np.inf, (1.0, 1.0)
    for al in EBG:
        ll = np.array([betabin_ll(x, m, al, be).sum() for be in EBG])
        j = int(np.argmax(ll))
        if ll[j] > best: best, ab = ll[j], (al, EBG[j])
    return ab


def run_B(genus, ids, h5dir, out):
    wd = f'{out}/work/B_{genus}'
    os.makedirs(wd, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    mats, comps, names = [], [], []
    with open(f'{wd}/reps.faa', 'w') as fo:
        for pid in ids:
            R, comp, part, M, mods = export(f'{h5dir}/{pid}.h5', f'{wd}/{pid}', proteins=True)
            if R.shape[1] > NB_MAX:
                keepg = np.sort(rng.choice(R.shape[1], NB_MAX, replace=False)); R = R.iloc[:, keepg]
            R = R[R.values.any(1)]
            fset = set(R.index)
            seqname, keepseq = None, False
            for l in open(f'{wd}/{pid}/fa/all_protein_families.faa'):
                if l.startswith('>'):
                    nm = l[1:].split()[0]; keepseq = nm in fset
                    if keepseq: fo.write(f'>{pid}|{nm}\n')
                elif keepseq: fo.write(l)
            mats.append(R); names.append(pid)
            log(f'{genus} {pid}: {R.shape[1]} genomes kept, {R.shape[0]} families')
            if not a.keep: shutil.rmtree(f'{wd}/{pid}/flat', ignore_errors=True)
    sh(f'mmseqs easy-cluster {wd}/reps.faa {wd}/mm {wd}/mmtmp --min-seq-id 0.5 -c 0.8 --cov-mode 0 '
       f'--threads {a.threads} -v 1', f'{wd}/mmseqs.log')
    cl = pd.read_csv(f'{wd}/mm_cluster.tsv', sep='\t', header=None, names=['rep', 'mem'])
    cid = {r: i for i, r in enumerate(cl.rep.unique())}
    mem2c = dict(zip(cl.mem, cl.rep.map(cid)))
    nC = len(cid)
    log(f'{genus}: {sum(m.shape[0] for m in mats)} species families -> {nC} genus clusters')
    # per species: clusters x genomes presence (OR over the species' families in a cluster)
    G = []
    for pid, R in zip(names, mats):
        idx = np.array([mem2c[f'{pid}|{f}'] for f in R.index])
        P = np.zeros((nC, R.shape[1]), bool)
        np.logical_or.at(P, idx, R.values)
        G.append(P)
    Nn = np.array([g.shape[1] for g in G])
    X = np.array([g.sum(1) for g in G], float)                 # species x clusters
    nsp = len(G)
    # how many clusters are shared between species (sanity for the mapping)
    shared = (X > 0).sum(0)
    summ, krows = [], []
    for f in range(nsp):
        others = [s for s in range(nsp) if s != f]
        kap, ll = fit_kappa(X, Nn, others)
        pg = np.clip(genus_mean(X, Nn, others), 1e-6, 1 - 1e-6)
        krows.append(dict(genus=genus, focal=names[f], kappa=kap, **{f'll_k{k:g}': round(v, 1) for k, v in zip(KGRID, ll)}))
        full = X[f] / Nn[f]
        strata = {'all': np.ones(nC, bool), 'core': full >= 0.95, 'accessory': (full > 0.05) & (full < 0.95),
                  'rare_or_absent': full <= 0.05}
        incong = X[others].sum(0) > 0
        strata.update({'shared': (full > 0) & incong, 'focal_only': (full > 0) & ~incong,
                       'congener_only': (full == 0) & incong})
        for m in MS_TRAIN:
            if Nn[f] < m + 5: continue
            for rep in range(REPS):
                perm = rng.permutation(Nn[f]); tr = perm[:m]; te = perm[m:m + NTEST]
                x = G[f][:, tr].sum(1).astype(float); Y = G[f][:, te].astype(float)   # clusters x test
                al, be = fit_beta_eb(x, m)
                preds = {'sp': (x + 1) / (m + 2), 'sp_eb': (x + al) / (m + al + be), 'genus': pg,
                         'shrunk': (x + kap * pg) / (m + kap)}
                preds['mix'] = 0.5 * preds['shrunk'] + 0.5 * preds['sp_eb']
                orc = {k: (x + k * pg) / (m + k) for k in KGRID}
                for sname, mask in strata.items():
                    if not mask.any(): continue
                    row = dict(genus=genus, focal=names[f], n_focal=int(Nn[f]), n_species=nsp, m=m, rep=rep,
                               n_test=len(te), stratum=sname, n_clusters=int(mask.sum()), kappa=kap)
                    for pn, p in preds.items():
                        row[f'bits_{pn}'] = round(bits(p[mask][:, None], Y[mask]).sum(0).mean(), 3)
                    ob = {k: bits(p[mask][:, None], Y[mask]).sum(0).mean() for k, p in orc.items()}
                    kb = min(ob, key=ob.get)
                    row['bits_oracle'] = round(ob[kb], 3); row['kappa_oracle'] = kb
                    summ.append(row)
        log(f'{genus} focal {names[f]} (n={Nn[f]}): kappa={kap:g}')
    os.makedirs(f'{out}/B', exist_ok=True)
    S = pd.DataFrame(summ)
    S['clusters_total'] = nC; S['clusters_multi_species'] = int((shared >= 2).sum())
    S.to_csv(f'{out}/B/{genus}.summary.tsv', sep='\t', index=False)
    pd.DataFrame(krows).to_csv(f'{out}/B/{genus}.kappa.tsv', sep='\t', index=False)
    open(f'{out}/B/{genus}.done', 'w').close()
    if not a.keep: shutil.rmtree(wd, ignore_errors=True)


if a.mode == 'species':
    run_A(a.args[0], a.args[1], a.args[2])
else:
    run_B(a.args[0], a.args[1].split(','), a.args[2], a.args[3])
import resource
log(f"done; max RSS self {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6:.2f} GB, "
    f"children {resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1e6:.2f} GB")
