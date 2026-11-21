"""Between-species tree vs gene content: which relatedness predicts a held-out species' gene content?

Tips are species (PanGBank pangenomes) inside one clade (a GTDB family, or class Nitrososphaeria). Two steps:

  worker.py prep ID H5 OUTDIR
      From the PanGBank .h5: per gene family its representative protein, partition and frequency (fraction of the
      species' genomes carrying it); genome GC and CheckM2 completeness (medians). Families in < 2 genomes or
      < 5% of genomes are dropped. Writes OUTDIR/prep/ID.{faa,tsv,json} (json last: marks prep as done).

  worker.py clade CLADE TREE.nwk OUTDIR --ids ID,ID,... [--threads 8]
      1. Orthogroups (OGs): MMseqs2 easy-cluster of all family representatives of the clade (40% identity,
         80% coverage both ways). F[s, og] = highest frequency among species s's families in the OG;
         presence y = F >= 0.5 (the OG is carried by most of the species' genomes).
      2. OG homology neighbours: representatives searched against each other (e 1e-3); a target OG's direct hits
         are masked with it (split families and paralogs would otherwise leak the target into the predictors).
      3. Relatedness between species:
           gtdb   patristic distance on the GTDB r226 bac120/ar53 tree, pruned to the clade (the ML marker tree;
                  its own topology is used for the tree predictor, not a UPGMA rebuild)
           aai    1 - AAI on core proteins (families in >= 90% of genomes): best hit per query and target species
                  (>= 50% query coverage), mean identity, averaged over both directions; UPGMA tree
           jac    1 - Jaccard of OG presence over variable OGs, with every target OG and its homology neighbours
                  removed; UPGMA tree
      4. Leave-one-species-out prediction of up to 400 variable OGs (present in 2 .. S-2 species):
           prev, ker_{gtdb,aai,jac}, tree_{gtdb,aai,jac}, relseq (stacked gtdb + aai), relall (+ jac),
           lr (context: other variable OGs, target and neighbours masked), avg_relseq, avg_relall
         (same models as tree_vs_context/worker.py; the tip model uses completeness 0.99, since a species-level
         majority call is nearly free of detection loss).
      Writes OUTDIR/res/CLADE.{species.tsv, targets.tsv.gz, diag.tsv, matcor.tsv, done}.
Bits are -log2 P(observed state), clipped at P >= 1e-6.
"""
import os, sys, time, json, argparse, subprocess, shutil
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('step', choices=['prep', 'clade'])
ap.add_argument('args', nargs='+')
ap.add_argument('--ids', default='')
ap.add_argument('--clades', default='', help='clades.tsv, for species and genus names')
ap.add_argument('--threads', type=int, default=8)
ap.add_argument('--max-targets', type=int, default=400)
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--keep', action='store_true')
a = ap.parse_args()
T0 = time.time()


def log(msg): print(f'[{time.time() - T0:7.0f}s] {msg}', flush=True)


def sh(cmd, logf):
    with open(logf, 'a') as fo:
        r = subprocess.run(cmd, shell=True, stdout=fo, stderr=subprocess.STDOUT)
    if r.returncode: raise RuntimeError(f'failed ({r.returncode}): {cmd}; see {logf}')


def bits(p, y):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return -(y * np.log2(p) + (1 - y) * np.log2(1 - p))


def logit(p):
    p = np.clip(p, 1e-4, 1 - 1e-4)
    return np.log(p / (1 - p))


# ============================================================================================ prep
def prep(pid, h5, out):
    import tables
    os.makedirs(f'{out}/prep', exist_ok=True)
    H = tables.open_file(h5)
    A = H.root.annotations
    genomes = np.array([x.decode() for x in A.genomes.col('name')])
    gidx = {g: i for i, g in enumerate(genomes)}
    ctg2g = dict(zip(A.contigs.col('ID').tolist(), [gidx[x.decode()] for x in A.contigs.col('genome')]))
    genes = pd.DataFrame({'gene': A.genes.col('ID'), 'contig': A.genes.col('contig')})
    genes['g'] = genes.contig.map(ctg2g)
    gf = pd.DataFrame({'gene': H.root.geneFamilies.col('gene'), 'fam': H.root.geneFamilies.col('geneFam')})
    genes = genes.merge(gf, on='gene')
    FI = H.root.geneFamiliesInfo
    fams = FI.col('name'); part = FI.col('partition'); prot = FI.col('protein')
    fpos = {f: i for i, f in enumerate(fams)}
    genes['f'] = genes.fam.map(fpos)
    pres = genes[['f', 'g']].drop_duplicates()
    cnt = np.bincount(pres.f.values, minlength=len(fams))
    n = len(genomes)
    freq = cnt / n
    md = {}
    try:
        M = H.root.metadata.genomes.pangbank_wf
        gc = np.array(M.col('gc_percentage'), float); ck = np.array(M.col('checkm2_completeness'), float)
        rep = [x.decode() for x, r in zip(M.col('ID'), M.col('gtdb_representative')) if r]
        md = dict(gc=float(np.nanmedian(gc)), gc_iqr=float(np.nanpercentile(gc, 75) - np.nanpercentile(gc, 25)),
                  checkm2=float(np.nanmedian(ck)), gtdb_rep=rep[0] if rep else '')
    except Exception as e:
        log(f'{pid}: no metadata ({e})')
    H.close()
    keep = np.flatnonzero((cnt >= 2) & (freq >= 0.05) | (np.char.startswith(part.astype(str), 'P')))
    with open(f'{out}/prep/{pid}.faa', 'w') as fo:
        for i in keep:
            fo.write(f'>{pid}|{i}\n{prot[i].decode().rstrip("*")}\n')
    pd.DataFrame({'i': keep, 'fam': [fams[i].decode() for i in keep],
                  'part': [part[i].decode()[:1] for i in keep], 'freq': np.round(freq[keep], 4)}
                 ).to_csv(f'{out}/prep/{pid}.tsv', sep='\t', index=False)
    md.update(id=pid, genomes=int(n), families=int(len(fams)), kept=int(len(keep)),
              persistent=int(np.char.startswith(part.astype(str), 'P').sum()))
    json.dump(md, open(f'{out}/prep/{pid}.json.tmp', 'w'))
    os.replace(f'{out}/prep/{pid}.json.tmp', f'{out}/prep/{pid}.json')
    log(f'{pid}: {n} genomes, {len(fams)} families, {len(keep)} kept')


# ============================================================================================ trees
MUH = np.array([0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000])
PIS = np.linspace(0.05, 0.95, 10)
GMU = np.repeat(MUH, len(PIS)); GPI = np.tile(PIS, len(MUH))
EPS = 0.01
QS = np.array([0.03, 0.1, 0.3, 1.0, 3.0])


def read_nwk(path, ids):
    """rooted binary Newick with tips named by id -> (children (N-1, 2), branch lengths, H) in the tip order of
    ids, plus the patristic distance matrix"""
    s = open(path).read().strip()
    parent, bl, name = [-1], [0.0], ['']
    cur = 0; i = 0
    while i < len(s):
        ch = s[i]
        if ch == '(':
            parent.append(cur); bl.append(0.0); name.append(''); cur = len(parent) - 1; i += 1
        elif ch == ',':
            parent.append(parent[cur]); bl.append(0.0); name.append(''); cur = len(parent) - 1; i += 1
        elif ch == ')': cur = parent[cur]; i += 1
        elif ch == ';': break
        elif ch == ':':
            j = i + 1
            while s[j] not in ',();': j += 1
            bl[cur] = float(s[i + 1:j]); i = j
        else:
            j = i
            while s[j] not in ':,();': j += 1
            name[cur] = s[i:j]; i = j
    M = len(parent)
    kids = [[] for _ in range(M)]
    for v in range(1, M): kids[parent[v]].append(v)
    tip_of = {name[v]: v for v in range(M) if not kids[v]}
    N = len(ids)
    assert set(tip_of) == set(map(str, ids)), 'tree tips differ from ids'
    newid = {tip_of[str(x)]: k for k, x in enumerate(ids)}
    order = []
    stack = [(0, False)]
    while stack:
        v, done = stack.pop()
        if done: order.append(v); continue
        stack.append((v, True))
        for c in kids[v]: stack.append((c, False))
    ch = []
    nxt = N
    for v in order:
        if kids[v]:
            assert len(kids[v]) == 2, 'tree must be binary'
            newid[v] = nxt; nxt += 1
            ch.append([newid[kids[v][0]], newid[kids[v][1]]])
    blen = np.zeros(2 * N - 1)
    depth = np.zeros(2 * N - 1)
    for v in range(M):
        blen[newid[v]] = bl[v] if v else 0.0
    ch = np.array(ch)
    for i in range(N - 2, -1, -1):           # root is the last internal node
        for c in ch[i]: depth[c] = depth[N + i] + blen[c]
    H = depth[:N].max()
    # patristic distances via lowest common ancestor (small N)
    par = np.full(2 * N - 1, -1)
    for i in range(N - 1):
        for c in ch[i]: par[c] = N + i
    anc = []
    for t in range(N):
        p = [t]
        while par[p[-1]] >= 0: p.append(par[p[-1]])
        anc.append(p)
    D = np.zeros((N, N))
    for x in range(N):
        sx = set(anc[x])
        for y in range(x + 1, N):
            lca = next(v for v in anc[y] if v in sx)
            D[x, y] = D[y, x] = depth[x] + depth[y] - 2 * depth[lca]
    blen[:N] = np.maximum(blen[:N], 0.01 * H)
    return (ch, blen, H), D, depth[:N]


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


def tree_predict(tree, Y, train, comp):
    """two-state gain/loss pruning with a Bayesian average over the (mu, pi) grid; see tree_vs_context/worker.py.
    Returns N x T P(observed present): posterior for test tips, leave-one-out posterior for training tips."""
    ch, bl, H = tree
    N, T = Y.shape; G = len(GMU); nn = 2 * N - 1
    c = np.clip(comp, 0.5, 0.99)[:, None]
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
    pz = (down[:N, ..., 1] * w[None]).sum(1)
    return np.clip(pz * c + (1 - pz) * EPS, 1e-6, 1 - 1e-6)


def kernel_predict(D, Y, tr, te, p0):
    med = np.median(D[np.triu_indices(len(D), 1)])
    med = med if med > 0 else max(D.max(), 1e-6)
    Ytr = Y[tr]
    best = None
    for q in QS:
        W = np.exp(-D / (q * med))
        Wtt = W[np.ix_(tr, tr)].copy(); np.fill_diagonal(Wtt, 0)
        loo = (Wtt @ Ytr + p0) / (Wtt.sum(1, keepdims=True) + 1)
        b = bits(loo, Ytr).mean(0)
        Wte = W[np.ix_(te, tr)]
        pte = (Wte @ Ytr + p0) / (Wte.sum(1, keepdims=True) + 1)
        if best is None:
            best = [b, loo, pte]
        else:
            m = b < best[0]
            best[0] = np.where(m, b, best[0]); best[1][:, m] = loo[:, m]; best[2][:, m] = pte[:, m]
    out = np.zeros(Y.shape); out[tr] = best[1]; out[te] = best[2]
    return out


def newton_lr(Z, y, C=1.0, it=30):
    X = np.column_stack([np.ones(len(Z)), Z]); w = np.zeros(X.shape[1])
    pen = np.full(X.shape[1], 1 / C); pen[0] = 0
    for _ in range(it):
        p = 1 / (1 + np.exp(-np.clip(X @ w, -30, 30)))
        g = X.T @ (p - y) + pen * w
        Hm = (X * (p * (1 - p))[:, None]).T @ X + np.diag(pen + 1e-9)
        step = np.linalg.solve(Hm, g); w -= step
        if np.abs(step).max() < 1e-6: break
    return w


def newton_predict(w, Z):
    return 1 / (1 + np.exp(-np.clip(w[0] + Z @ w[1:], -30, 30)))


def jaccard_dist(F):
    F = F.astype(np.float32)
    inter = F @ F.T; sz = F.sum(1)
    J = inter / np.maximum(sz[:, None] + sz[None, :] - inter, 1)
    D = 1 - J; np.fill_diagonal(D, 0)
    return D


# ============================================================================================ clade
def read_m8(fn, cols):
    return pd.read_csv(fn, sep='\t', header=None, names=cols)


def clade(cl, nwk, out):
    from sklearn.linear_model import LogisticRegression
    from scipy.stats import spearmanr
    ids = [int(x) for x in a.ids.split(',')]
    wd = f'{out}/work/{cl}'
    os.makedirs(wd, exist_ok=True); os.makedirs(f'{out}/res', exist_ok=True)
    S = len(ids)
    meta = [json.load(open(f'{out}/prep/{i}.json')) for i in ids]
    tabs = {i: pd.read_csv(f'{out}/prep/{i}.tsv', sep='\t') for i in ids}
    sidx = {i: k for k, i in enumerate(ids)}
    # ---- 1. orthogroups
    if not os.path.exists(f'{wd}/og_cluster.tsv'):
        with open(f'{wd}/all.faa', 'w') as fo:
            for i in ids: fo.write(open(f'{out}/prep/{i}.faa').read())
        sh(f'mmseqs easy-cluster {wd}/all.faa {wd}/og {wd}/tmp --min-seq-id 0.4 -c 0.8 --cov-mode 0 -e 1e-5 '
           f'--threads {a.threads} -v 1', f'{wd}/mmseqs.log')
    cl_t = read_m8(f'{wd}/og_cluster.tsv', ['rep', 'mem'])
    reps = pd.unique(cl_t.rep)
    rid = {r: k for k, r in enumerate(reps)}
    cl_t['og'] = cl_t.rep.map(rid)
    sp = cl_t.mem.str.split('|', expand=True)
    cl_t['s'] = sp[0].astype(int).map(sidx); cl_t['fi'] = sp[1].astype(int)
    fq = pd.concat([t.assign(s=sidx[i])[['s', 'i', 'freq', 'part']] for i, t in tabs.items()])
    cl_t = cl_t.merge(fq, left_on=['s', 'fi'], right_on=['s', 'i'])
    nog = len(reps)
    F = np.zeros((S, nog))
    np.maximum.at(F, (cl_t.s.values, cl_t.og.values), cl_t.freq.values)
    Y_all = F >= 0.5
    log(f'{cl}: {S} species, {len(cl_t)} families -> {nog} OGs')
    # ---- 2. homology neighbours of OG representatives
    if not os.path.exists(f'{wd}/rep_hits.m8'):
        sh(f'mmseqs easy-search {wd}/og_rep_seq.fasta {wd}/og_rep_seq.fasta {wd}/rep_hits.m8 {wd}/tmp2 -e 1e-3 '
           f'--format-output query,target --threads {a.threads} -v 1', f'{wd}/mmseqs.log')
    rh = read_m8(f'{wd}/rep_hits.m8', ['q', 't'])
    rh = rh[rh.q != rh.t]
    nbr = {}
    for q, t in zip(rh.q.map(rid).values, rh.t.map(rid).values):
        nbr.setdefault(q, set()).add(t); nbr.setdefault(t, set()).add(q)
    # ---- 3a. core AAI
    if not os.path.exists(f'{wd}/core_hits.m8'):
        with open(f'{wd}/core.faa', 'w') as fo:
            for i in ids:
                keep = set(tabs[i].i[tabs[i].freq >= 0.9].astype(str))
                rec = open(f'{out}/prep/{i}.faa').read().split('>')[1:]
                for r_ in rec:
                    h = r_.split('\n', 1)[0]
                    if h.split('|')[1] in keep: fo.write('>' + r_)
        sh(f'mmseqs easy-search {wd}/core.faa {wd}/core.faa {wd}/core_hits.m8 {wd}/tmp3 -e 1e-5 --max-seqs 1000 '
           f'--format-output query,target,pident,qcov --threads {a.threads} -v 1', f'{wd}/mmseqs.log')
    ch_ = read_m8(f'{wd}/core_hits.m8', ['q', 't', 'pid', 'qcov'])
    ch_ = ch_[ch_.qcov >= 0.5]
    ch_['qs'] = ch_.q.str.split('|').str[0].astype(int).map(sidx)
    ch_['ts'] = ch_.t.str.split('|').str[0].astype(int).map(sidx)
    ch_ = ch_[ch_.qs != ch_.ts]
    best = ch_.sort_values('pid', ascending=False).drop_duplicates(['q', 'ts'])
    aai = best.groupby(['qs', 'ts']).pid.mean()
    A_ = np.full((S, S), np.nan)
    A_[aai.index.get_level_values(0), aai.index.get_level_values(1)] = aai.values
    A_ = np.nanmean(np.stack([A_, A_.T]), 0)
    Daai = 1 - A_ / 100
    np.fill_diagonal(Daai, 0)
    if np.isnan(Daai).any():
        Daai[np.isnan(Daai)] = np.nanmax(Daai) * 1.5
    # ---- 3b. GTDB tree
    tree_g, Dg, r2t = read_nwk(nwk, ids)
    gc = np.array([m.get('gc', np.nan) for m in meta])
    # ---- 4. targets and prediction
    rng = np.random.default_rng(a.seed)
    k = Y_all.sum(0)
    var = np.flatnonzero((k >= 2) & (k <= S - 2))
    tg = var if len(var) <= a.max_targets else np.sort(rng.choice(var, a.max_targets, replace=False))
    mask = set(tg.tolist())
    for t in tg: mask |= nbr.get(t, set())
    jfeat = np.array([o for o in var if o not in mask])
    Djac = jaccard_dist(Y_all[:, jfeat])
    D = {'gtdb': Dg, 'aai': Daai, 'jac': Djac}
    trees = {'gtdb': tree_g, 'aai': upgma(Daai), 'jac': upgma(Djac)}
    meas = list(D)
    rel_seq = [f'{kk}_{m}' for m in ('gtdb', 'aai') for kk in ('ker', 'tree')]
    rel_all = rel_seq + ['ker_jac', 'tree_jac']
    preds = ['prev'] + rel_all + ['relseq', 'relall', 'lr', 'avg_relseq', 'avg_relall']
    Y = Y_all[:, tg].astype(float)
    T_ = len(tg)
    log(f'{cl}: {len(var)} variable OGs, {T_} targets, {len(jfeat)} OGs left for jac; '
        f'median d gtdb {np.median(Dg[np.triu_indices(S, 1)]):.3f} aai {np.median(Daai[np.triu_indices(S, 1)]):.3f} '
        f'jac {np.median(Djac[np.triu_indices(S, 1)]):.3f}')
    Fv = Y_all[:, var].astype(np.float32)
    vpos = {o: j for j, o in enumerate(var)}
    comp = np.full(S, 0.99)
    P = {m: np.zeros((S, T_)) for m in preds}
    for fo in range(S):
        te = np.zeros(S, bool); te[fo] = True; tr = ~te
        Ytr = Y[tr]; p0 = (Ytr.sum(0) + 1) / (tr.sum() + 2)
        P['prev'][te] = p0
        inner = {}
        for m in meas:
            kp = kernel_predict(D[m], Y, tr, te, p0)
            tp = tree_predict(trees[m], Y, tr, comp)
            P[f'ker_{m}'][te] = kp[te]; P[f'tree_{m}'][te] = tp[te]
            inner[f'ker_{m}'] = kp; inner[f'tree_{m}'] = tp
        for t in range(T_):
            ytr = Ytr[:, t]
            if ytr.min() == ytr.max():
                for m in ('relseq', 'relall', 'lr'): P[m][te, t] = p0[t]
                continue
            for nm, rl in (('relseq', rel_seq), ('relall', rel_all)):
                Zm = np.column_stack([logit(inner[r][:, t]) for r in rl])
                P[nm][te, t] = newton_predict(newton_lr(Zm[tr], ytr), Zm[te])
            drop = {tg[t]} | nbr.get(tg[t], set())
            keep = np.ones(len(var), bool)
            keep[[vpos[o] for o in drop if o in vpos]] = False
            X = Fv[:, keep]
            P['lr'][te, t] = LogisticRegression(C=0.1, max_iter=300).fit(X[tr], ytr).predict_proba(X[te])[:, 1]
        if fo % 10 == 0: log(f'{cl}: fold {fo + 1}/{S}')
    for nm in ('relseq', 'relall'):
        P[f'avg_{nm}'] = 0.5 * P['lr'] + 0.5 * P[nm]
    B = {m: bits(P[m], Y) for m in preds}
    names = [m.get('id') for m in meta]
    C = pd.read_csv(a.clades, sep='\t').set_index('id') if a.clades else None
    sprow = []
    for s in range(S):
        r = dict(clade=cl, id=ids[s], species=C.species[ids[s]] if C is not None else '',
                 genus=C.genus[ids[s]] if C is not None else '', genomes=meta[s]['genomes'],
                 gc=round(gc[s], 2), n_present=int(Y[s].sum()))
        r.update({f'bits_{m}': round(B[m][s].mean(), 5) for m in preds})
        sprow.append(r)
    SP = pd.DataFrame(sprow)
    SP.to_csv(f'{out}/res/{cl}.species.tsv', sep='\t', index=False)
    trows = [dict(clade=cl, og=int(tg[t]), carriers=int(Y[:, t].sum()),
                  **{f'bits_{m}': round(B[m][:, t].mean(), 5) for m in preds}) for t in range(T_)]
    pd.DataFrame(trows).to_csv(f'{out}/res/{cl}.targets.tsv.gz', sep='\t', index=False)
    # ---- diagnostics
    gcdev = gc - np.nanmedian(gc)
    dg = dict(SP[['clade', 'id', 'species', 'genus', 'gc']])
    dg['gc_dev'] = np.round(gcdev, 2)
    dg['r2t_gtdb_z'] = np.round((r2t - r2t.mean()) / (r2t.std() + 1e-12), 3)
    nnb = {}
    for m, Dm in D.items():
        Dx = Dm + np.diag(np.full(S, np.inf))
        nn = Dx.argmin(1)
        md_ = Dm.sum(1) / (S - 1)
        dg[f'nn_{m}'] = [SP.species[j] for j in nn]
        dg[f'nnd_{m}'] = np.round(Dx.min(1), 4)
        dg[f'nnpct_{m}'] = np.round(pd.Series(Dx.min(1)).rank(pct=True).values, 3)
        dg[f'meanz_{m}'] = np.round((md_ - md_.mean()) / (md_.std() + 1e-12), 3)
        nnb[m] = np.argsort(Dx, 1, kind='stable')[:, :min(5, S - 1)]
        # how far, in GC, the nearest neighbour is
        dg[f'nn_gcdiff_{m}'] = np.round(np.abs(gc - gc[nn]), 2)
    for x, y in (('gtdb', 'jac'), ('aai', 'jac'), ('gtdb', 'aai')):
        dg[f'knn5_{x}_{y}'] = [len(set(p) & set(q)) / len(p) for p, q in zip(nnb[x], nnb[y])]
        dg[f'rowcor_{x}_{y}'] = [round(spearmanr(np.delete(D[x][s], s), np.delete(D[y][s], s)).statistic, 3)
                                 for s in range(S)]
    DG = pd.DataFrame(dg)
    DG.to_csv(f'{out}/res/{cl}.diag.tsv', sep='\t', index=False)
    iu = np.triu_indices(S, 1)
    dGC = np.abs(gc[:, None] - gc[None, :])
    mc = []
    for x, y in (('gtdb', 'aai'), ('gtdb', 'jac'), ('aai', 'jac')):
        mc.append(dict(clade=cl, a=x, b=y, spearman=round(spearmanr(D[x][iu], D[y][iu]).statistic, 4)))
    for m in meas:
        mc.append(dict(clade=cl, a='dGC', b=m, spearman=round(spearmanr(dGC[iu], D[m][iu]).statistic, 4)))
    # GC signal in each sequence distance beyond gene content: residual of rank(d) on rank(jac) vs dGC
    from scipy.stats import rankdata
    rj = rankdata(Djac[iu])
    for m in ('gtdb', 'aai'):
        rm = rankdata(D[m][iu]); bcoef = np.polyfit(rj, rm, 1)
        res = rm - np.polyval(bcoef, rj)
        mc.append(dict(clade=cl, a='dGC', b=f'{m}|jac_resid', spearman=round(spearmanr(dGC[iu], res).statistic, 4)))
    pd.DataFrame(mc).to_csv(f'{out}/res/{cl}.matcor.tsv', sep='\t', index=False)
    open(f'{out}/res/{cl}.done', 'w').close()
    log(f'{cl}: ' + ' '.join(f'{m}={B[m].mean():.3f}' for m in preds))
    if not a.keep: shutil.rmtree(f'{wd}/tmp', ignore_errors=True); shutil.rmtree(f'{wd}/tmp2', ignore_errors=True); \
        shutil.rmtree(f'{wd}/tmp3', ignore_errors=True)


if a.step == 'prep':
    prep(int(a.args[0]), a.args[1], a.args[2])
else:
    clade(a.args[0], a.args[1], a.args[2])
import resource
log(f"done; max RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6:.2f} GB")
