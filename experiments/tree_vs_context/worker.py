"""Tree vs context: how much of the context model's gain is relatedness, and which relatedness measure.

One species per task, three steps (the middle one runs skani from another container):

  worker.py prep  ID H5 OUTDIR [--pool 400] [--seed 1]
      Reads the PanGBank .h5 directly (PyTables): genomes, gene -> family, partitions, modules, gene DNA,
      genome GC (PanGBank metadata). Draws a pool of at most --pool genomes. Writes into OUTDIR/work/ID/:
        pool.npz      presence/absence (families x pool genomes), partitions, module members, GC (written last;
                      marks prep as done)
        cds/<g>.fna.gz  per genome: its CDS concatenated per contig in gene order (input for skani)
        cds_list.txt  list of those files
        aai.npy       core-protein distance: Mash-style distance on FracMinHash (1/32) amino-acid 5-mers of
                      persistent families, k-mers keyed by family (estimates 1 - AAI over the core)
  (skani)  skani triangle -l cds_list.txt -E --min-af 0 -o skani.tsv   -> ANI on coding sequence
  worker.py model ID OUTDIR [--seed 1]
      Builds the ANI distance from skani.tsv, then for subsamples N in 25, 50, 100, 200, 400 (3 random
      subsamples of the pool each, 1 if N = pool) scores held-out presence of up to 60 modules and 60
      accessory families by 10-fold cross-validation (same folds for every predictor):
        prev           Beta(1,1) prevalence
        ker_{ani,aai,jac}   relatedness-weighted neighbours, w = exp(-d / (q * median d)), q chosen per target
                       and fold by leave-one-out on the training genomes (q in 0.03 .. 3)
        tree_{ani,aai,jac}  UPGMA tree from the distance, two-state gain/loss chain (stationary root),
                       per-target Bayesian average over a (mu, pi) grid weighted by the marginal likelihood
                       of the training tips (Felsenstein pruning; test tips have uncertain tips [1, 1]);
                       up-down pass gives test predictions and leave-one-out predictions for training tips;
                       detection model P(observed | present) = completeness (capped 0.99),
                       P(observed | absent) = 0.01
        relseq         stacked relatedness: logistic (C = 1, Newton) on logits of ker/tree ani+aai
        relall         same plus ker/tree jac
        lr             context only: L2 logistic (C = 0.1) on accessory families (target's own module masked)
                       + completeness, as in pangbank_scale/worker.py
        lr_relseq      context + sequence relatedness: the lr model plus the top min(10, N/10) principal
                       coordinates (PCoA, standardised) of each sequence distance (ani, aai); these carry
                       no outcome information, so training and test rows are treated alike
        lr_relall      same with the jac coordinates added
        avg_relseq, avg_relall   untuned 50/50 probability average of lr and relseq / relall
      'jac' is gene-content Jaccard on accessory families with ALL target families (and their modules)
      removed from the distance for that subsample, so it never sees the targets.
      Writes OUTDIR/res/ID.summary.tsv (per N, rep, kind: mean bits per predictor), ID.targets.tsv.gz,
      ID.genomes.tsv.gz (per subsample genome mean bits per predictor), ID.diag.tsv (per genome
      relatedness diagnostics on the pool), ID.matcor.tsv (correlations between distance matrices),
      ID.done.
Bits are -log2 P(observed state), clipped at P >= 1e-6.
"""
import os, sys, time, gzip, argparse, shutil
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('step', choices=['prep', 'model'])
ap.add_argument('args', nargs='+')
ap.add_argument('--pool', type=int, default=400)
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--keep', action='store_true')
a = ap.parse_args()
T0 = time.time()


def log(msg): print(f'[{time.time() - T0:7.0f}s] {msg}', flush=True)


def bits(p, y):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return -(y * np.log2(p) + (1 - y) * np.log2(1 - p))


def logit(p):
    p = np.clip(p, 1e-4, 1 - 1e-4)
    return np.log(p / (1 - p))


# ============================================================================================ prep
CODON = {}
_b = 'TCAG'
_aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'   # table 11
AAIDX = {c: i for i, c in enumerate('ACDEFGHIKLMNPQRSTVWY')}
LUT = np.full(64, 20, np.int64)
for i, c1 in enumerate(_b):
    for j, c2 in enumerate(_b):
        for k, c3 in enumerate(_b):
            aa = _aa[16 * i + 4 * j + k]
            LUT[16 * 'TCAG'.index(c1) + 4 * 'TCAG'.index(c2) + 'TCAG'.index(c3)] = AAIDX.get(aa, 20)
NT = np.full(256, 4, np.int64)
for c, v in zip(b'TCAG', range(4)): NT[c] = v
for c, v in zip(b'tcag', range(4)): NT[c] = v
KAA = 5
SCALE = 32


def mix64(x):
    x = x.astype(np.uint64)
    x ^= x >> np.uint64(31); x *= np.uint64(0x7FB5D329728EA185)
    x ^= x >> np.uint64(27); x *= np.uint64(0x81DADEF4BC2DD44D)
    x ^= x >> np.uint64(33)
    return x


def prot_kmers(seq, fam):
    """sampled hashed amino-acid 5-mers of one CDS (frame 1), keyed by family"""
    s = NT[np.frombuffer(seq, np.uint8)]
    L = (len(s) // 3) * 3
    if L < 3 * (KAA + 1): return np.empty(0, np.uint64)
    c = s[:L].reshape(-1, 3)
    bad = (c == 4).any(1)
    aa = LUT[(c[:, 0] * 16 + c[:, 1] * 4 + c[:, 2]) % 64]; aa[bad] = 20
    n = len(aa) - KAA + 1
    code = np.zeros(n, np.int64); ok = np.ones(n, bool)
    for j in range(KAA):
        w = aa[j:j + n]; code = code * 21 + w; ok &= w < 20
    h = mix64(code[ok] + np.int64(fam) * np.int64(21 ** KAA + 7))
    return h[h % np.uint64(SCALE) == 0]


def prep(pid, h5, out):
    import tables
    wd = f'{out}/work/{pid}'
    os.makedirs(f'{wd}/cds', exist_ok=True)
    rng = np.random.default_rng(a.seed)
    H = tables.open_file(h5)
    A = H.root.annotations
    genomes_all = np.array([x.decode() for x in A.genomes.col('name')])
    n = len(genomes_all)
    pool = np.sort(rng.choice(n, min(n, a.pool), replace=False))
    gpool = genomes_all[pool]
    gidx = {g: i for i, g in enumerate(gpool)}
    log(f'{pid}: {n} genomes, pool {len(gpool)}')
    ctg_id = A.contigs.col('ID'); ctg_gen = A.contigs.col('genome')
    ctg2g = dict(zip(ctg_id.tolist(), [gidx.get(x.decode(), -1) for x in ctg_gen]))
    genes = pd.DataFrame({'gene': A.genes.col('ID'), 'contig': A.genes.col('contig'),
                          'gd': A.genes.col('genedata_id')})
    genes['g'] = genes.contig.map(ctg2g).astype(np.int64)
    genes = genes[genes.g >= 0]
    pos = A.genedata.col('position'); gdid = A.genedata.col('genedata_id')
    posmap = np.zeros(gdid.max() + 1, np.int64); posmap[gdid] = pos
    genes['pos'] = posmap[genes.gd.values]
    gs = pd.DataFrame({'gene': A.geneSequences.col('gene'), 'seqid': A.geneSequences.col('seqid')})
    gf = pd.DataFrame({'gene': H.root.geneFamilies.col('gene'), 'fam': H.root.geneFamilies.col('geneFam')})
    genes = genes.merge(gs, on='gene').merge(gf, on='gene')
    fi_name = H.root.geneFamiliesInfo.col('name'); fi_part = H.root.geneFamiliesInfo.col('partition')
    fams = np.array([x.decode() for x in fi_name])
    part = np.array([{'P': 'persistent', 'S': 'shell', 'C': 'cloud'}.get(x.decode()[:1], '?') for x in fi_part])
    fpos = {f.encode(): i for i, f in enumerate(fams)}
    genes['f'] = genes.fam.map(fpos).astype(np.int64)
    log(f'{pid}: {len(genes)} pool genes, {len(fams)} families')
    # presence / absence
    R = np.zeros((len(fams), len(gpool)), bool)
    R[genes.f.values, genes.g.values] = True
    keepf = R.any(1)
    # modules
    mods = {}
    if '/modules' in H:
        for f, m in zip(H.root.modules.col('geneFam'), H.root.modules.col('module')):
            if f in fpos: mods.setdefault(int(m), []).append(fpos[f])
    # genome GC (PanGBank metadata) and checkm2 completeness
    gc = np.full(len(gpool), np.nan); ck = np.full(len(gpool), np.nan)
    try:
        md = H.root.metadata.genomes.pangbank_wf
        for g, v, c in zip(md.col('ID'), md.col('gc_percentage'), md.col('checkm2_completeness')):
            i = gidx.get(g.decode())
            if i is not None: gc[i] = v; ck[i] = c
    except Exception as e:
        log(f'no GC metadata ({e}); using CDS GC')
    # sequences: stream the deduplicated DNA table, keep what the pool needs
    need = np.zeros(int(genes.seqid.max()) + 1, bool); need[genes.seqid.values] = True
    S = A.sequences
    seqs = {}
    step = 20000
    for st in range(0, S.nrows, step):
        blk = S.read(st, min(st + step, S.nrows))
        for dna, sid in zip(blk['dna'], blk['seqid']):
            if sid < len(need) and need[sid]: seqs[int(sid)] = dna
    H.close()
    log(f'{pid}: {len(seqs)} unique sequences loaded')
    # CDS fasta per genome (contigs in gene order, CDS concatenated) and CDS GC as fallback
    genes = genes.sort_values(['g', 'contig', 'pos'])
    files = []
    gc_cds = np.zeros(len(gpool))
    for g, sub in genes.groupby('g', sort=True):
        fn = f'{wd}/cds/{g}.fna.gz'
        nb = gcb = 0
        with gzip.open(fn, 'wb', compresslevel=1) as fo:
            for ctg, s2 in sub.groupby('contig', sort=False):
                seq = b''.join(seqs[i] for i in s2.seqid.values)
                fo.write(b'>c%d\n' % ctg); fo.write(seq); fo.write(b'\n')
                nb += len(seq); gcb += seq.count(b'G') + seq.count(b'C')
        gc_cds[g] = 100 * gcb / max(nb, 1)
        files.append(os.path.abspath(fn))
    with open(f'{wd}/cds_list.txt', 'w') as fo: fo.write('\n'.join(files) + '\n')
    gc = np.where(np.isnan(gc), gc_cds, gc)
    # core-protein k-mer distance (persistent families)
    pers = genes[part[genes.f.values] == 'persistent']
    sid_u = np.unique(pers.seqid.values)
    sid_fam = pers.drop_duplicates('seqid').set_index('seqid').f
    kms = {int(s): prot_kmers(seqs[int(s)], int(sid_fam[s])) for s in sid_u}
    del seqs
    rows, cols = [], []
    for g, sub in pers.groupby('g'):
        h = np.unique(np.concatenate([kms[int(s)] for s in np.unique(sub.seqid.values)]))
        rows.append(np.full(len(h), g, np.int64)); cols.append(h)
    rows = np.concatenate(rows); cols = np.concatenate(cols)
    _, ci = np.unique(cols, return_inverse=True)
    import scipy.sparse as sp
    X = sp.csr_matrix((np.ones(len(rows), np.float32), (rows, ci)), shape=(len(gpool), ci.max() + 1))
    inter = (X @ X.T).toarray(); sz = np.asarray(X.sum(1)).ravel()
    J = inter / np.maximum(sz[:, None] + sz[None, :] - inter, 1)
    D = -np.log(np.clip(2 * J / (1 + J), 1e-9, 1)) / KAA
    np.fill_diagonal(D, 0)
    np.save(f'{wd}/aai.npy', D.astype(np.float32))
    log(f'{pid}: core-protein distance from {X.shape[1]} sampled k-mers (median {np.median(sz):.0f}/genome); '
        f'median d {np.median(D[np.triu_indices(len(gpool), 1)]):.4f}')
    mods = {m: v for m, v in mods.items() if keepf[v].any()}
    newidx = np.cumsum(keepf) - 1
    np.savez_compressed(f'{wd}/pool.tmp.npz', R=R[keepf], fams=fams[keepf], part=part[keepf], genomes=gpool,
                        gc=gc, checkm2=ck, n_total=n,
                        mod_ids=np.array(list(mods), dtype=np.int64),
                        mod_members=np.array([','.join(str(newidx[x]) for x in v if keepf[x]) for v in mods.values()],
                                             dtype=object) if mods else np.array([], object))
    os.replace(f'{wd}/pool.tmp.npz', f'{wd}/pool.npz')        # written last: its presence marks prep as done
    log(f'{pid}: prep done')


# ============================================================================================ model
SIZES = (25, 50, 100, 200, 400)
MAX_TGT = 60
QS = np.array([0.03, 0.1, 0.3, 1.0, 3.0])
MUH = np.array([0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000])   # turnover rate x root height (1000 ~ no tree signal)
PIS = np.linspace(0.05, 0.95, 10)
GMU = np.repeat(MUH, len(PIS)); GPI = np.tile(PIS, len(MUH))
EPS = 0.01                                            # P(observed present | absent): annotation/clustering noise


def read_skani(fn, files):
    idx = {os.path.basename(f): i for i, f in enumerate(files)}
    n = len(files)
    D = np.full((n, n), np.nan); np.fill_diagonal(D, 0)
    t = pd.read_csv(fn, sep='\t')
    r = t.Ref_file.map(lambda x: idx[os.path.basename(x)]).values
    q = t.Query_file.map(lambda x: idx[os.path.basename(x)]).values
    d = 1 - t.ANI.values / 100
    D[r, q] = d; D[q, r] = d
    miss = np.isnan(D)
    if miss.any():
        fill = np.nanmax(D) * 1.5 if np.isfinite(np.nanmax(D)) and np.nanmax(D) > 0 else 0.05
        D[miss] = fill
    return D, int(miss.sum() // 2)


def upgma(D):
    """tips 0..N-1, internal N..2N-2 (root last); returns children (N-1, 2), branch length per node"""
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
    # floor on tip branches: identical sequences (d = 0) would otherwise force identical content at every rate,
    # leaving no way for the rate grid to say that relatedness does not inform a target
    bl[:N] = np.maximum(bl[:N], 0.01 * H)
    return ch, bl, H


def tree_predict(tree, Y, train, comp):
    """two-state gain/loss pruning, Bayesian average over (mu, pi) grid.
    Y: N x T observed (0/1); train: bool N; comp: N detection prob. Returns N x T predicted P(observed present):
    for test tips the posterior given all training tips; for training tips the leave-one-out posterior (outside
    message only; grid weights include the tip)."""
    ch, bl, H = tree
    N, T = Y.shape; G = len(GMU); nn = 2 * N - 1
    c = np.clip(comp, 0.5, 0.99)[:, None]
    L = np.ones((N, T, 2))
    L[train, :, 0] = np.where(Y[train] > 0, EPS, 1 - EPS)
    L[train, :, 1] = np.where(Y[train] > 0, c[train], 1 - c[train])
    mu = GMU[:, None] / H; pi = GPI[:, None]                     # G x 1
    up = np.empty((nn, G, T, 2)); msg = np.empty((nn, G, T, 2)); lsc = np.zeros((G, T))
    up[:N] = L[:, None]
    def message(v):
        e = np.exp(-mu * bl[v])                                   # G x 1
        u0, u1 = up[v, ..., 0], up[v, ..., 1]
        m0 = (1 - pi * (1 - e)) * u0 + pi * (1 - e) * u1
        m1 = (1 - pi) * (1 - e) * u0 + (pi + (1 - pi) * e) * u1
        msg[v, ..., 0] = m0; msg[v, ..., 1] = m1
    for v in range(N): message(v)
    for i in range(N - 1):
        v = N + i; l, r = ch[i]
        u = msg[l] * msg[r]
        s = u.max(-1, keepdims=True); s[s <= 0] = 1
        up[v] = u / s; lsc += np.log(s[..., 0])
        if v < nn - 1: message(v)
    root = nn - 1
    ll = lsc + np.log(np.maximum((1 - pi) * up[root, ..., 0] + pi * up[root, ..., 1], 1e-300))
    w = np.exp(ll - ll.max(0, keepdims=True)); w /= w.sum(0, keepdims=True)      # G x T
    down = np.empty((nn, G, T, 2))
    down[root, ..., 0] = np.broadcast_to(1 - pi, (G, T)); down[root, ..., 1] = np.broadcast_to(pi, (G, T))
    for i in range(N - 2, -1, -1):
        v = N + i; l, r = ch[i]
        for cnode, sib in ((l, r), (r, l)):
            x = down[v] * msg[sib]                                # state at v, data outside cnode
            e = np.exp(-mu * bl[cnode])
            d0 = x[..., 0] * (1 - pi * (1 - e)) + x[..., 1] * (1 - pi) * (1 - e)
            d1 = x[..., 0] * pi * (1 - e) + x[..., 1] * (pi + (1 - pi) * e)
            s = np.maximum(d0 + d1, 1e-300)
            down[cnode, ..., 0] = d0 / s; down[cnode, ..., 1] = d1 / s
    pz = (down[:N, ..., 1] * w[None]).sum(1)                       # N x T
    return np.clip(pz * c + (1 - pz) * EPS, 1e-6, 1 - 1e-6), w


def kernel_predict(D, Y, tr, te, p0):
    """returns test predictions and training leave-one-out predictions (N x T), q chosen per target by
    training LOO"""
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
    """small L2 logistic regression (intercept unpenalised) by Newton's method; returns coefficients"""
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


def pcoa(D, k=10):
    """classical MDS: top k principal coordinates of a distance matrix, standardised"""
    n = len(D); J = np.eye(n) - 1 / n
    B = -0.5 * J @ (D ** 2) @ J
    ev, V = np.linalg.eigh(B)
    o = np.argsort(ev)[::-1][:k]; ev = np.maximum(ev[o], 0)
    X = V[:, o] * np.sqrt(ev)
    return ((X - X.mean(0)) / (X.std(0) + 1e-12)).astype(np.float32)


def jaccard_dist(F):
    F = F.astype(np.float32)
    inter = F @ F.T; sz = F.sum(1)
    J = inter / np.maximum(sz[:, None] + sz[None, :] - inter, 1)
    D = 1 - J; np.fill_diagonal(D, 0)
    return D


def spearman_pairs(A, B):
    from scipy.stats import spearmanr
    iu = np.triu_indices(len(A), 1)
    return spearmanr(A[iu], B[iu]).statistic


def diagnostics(pid, Ds, gc, genomes, out):
    """per genome: nearest-neighbour distance percentile and mean-distance z under each measure, k=5
    neighbour overlap between measures, GC deviation from the species median"""
    n = len(genomes)
    rows = {'genome': genomes, 'gc': np.round(gc, 3), 'gc_dev': np.round(gc - np.nanmedian(gc), 3)}
    nbr = {}
    for m, D in Ds.items():
        Dm = D + np.diag(np.full(n, np.inf))
        nn = Dm.min(1)
        rows[f'nn_{m}'] = np.round(nn, 6)
        rows[f'nnpct_{m}'] = np.round(pd.Series(nn).rank(pct=True).values, 3)
        md = D.sum(1) / (n - 1)
        rows[f'meanz_{m}'] = np.round((md - md.mean()) / (md.std() + 1e-12), 3)
        k = min(5, n - 1)
        nbr[m] = np.argsort(Dm, 1, kind='stable')[:, :k]
    ms = list(Ds)
    for i in range(len(ms)):
        for j in range(i + 1, len(ms)):
            rows[f'knn5_{ms[i]}_{ms[j]}'] = [len(set(x) & set(y)) / len(x) for x, y in zip(nbr[ms[i]], nbr[ms[j]])]
    T = pd.DataFrame(rows); T.insert(0, 'id', pid)
    # discordance flags: outlying under a sequence measure but typical under gene content, and vice versa
    for s in ('ani', 'aai'):
        if f'meanz_{s}' not in T: continue
        seq_out = (T[f'meanz_{s}'] > 2) | (T[f'nnpct_{s}'] >= 0.95)
        con_out = (T['meanz_jac'] > 2) | (T['nnpct_jac'] >= 0.95)
        T[f'flag_{s}_out_jac_in'] = seq_out & (T.nnpct_jac <= 0.5) & (T.meanz_jac < 1)
        T[f'flag_jac_out_{s}_in'] = con_out & (T[f'nnpct_{s}'] <= 0.5) & (T[f'meanz_{s}'] < 1)
    T.to_csv(f'{out}/res/{pid}.diag.tsv', sep='\t', index=False)
    mc = []
    for i in range(len(ms)):
        for j in range(i + 1, len(ms)):
            mc.append(dict(id=pid, n=n, a=ms[i], b=ms[j], spearman=round(spearman_pairs(Ds[ms[i]], Ds[ms[j]]), 4)))
    gcd = np.abs(T.gc_dev.values)
    for m in ms:
        from scipy.stats import spearmanr
        mc.append(dict(id=pid, n=n, a='absGCdev', b=f'meanz_{m}',
                       spearman=round(spearmanr(gcd, T[f'meanz_{m}']).statistic, 4) if np.nanstd(gcd) > 0 else np.nan))
    pd.DataFrame(mc).to_csv(f'{out}/res/{pid}.matcor.tsv', sep='\t', index=False)


def model(pid, out):
    from sklearn.linear_model import LogisticRegression
    wd = f'{out}/work/{pid}'
    Z = np.load(f'{wd}/pool.npz', allow_pickle=True)
    R = Z['R']; fams = Z['fams']; part = Z['part']; genomes = Z['genomes']; gc = Z['gc']
    n = R.shape[1]; n_total = int(Z['n_total'])
    mod_ids = list(Z['mod_ids']); mod_mem = [np.array([int(x) for x in s.split(',')]) for s in Z['mod_members']]
    fam2mod = {}
    for mi, mem in enumerate(mod_mem):
        for f in mem: fam2mod[f] = mi
    M = np.array([R[mem].mean(0) >= 0.5 for mem in mod_mem]) if mod_mem else np.zeros((0, n), bool)
    pers = part == 'persistent'
    comp = R[pers].mean(0) if pers.any() else np.ones(n)
    files = [l.strip() for l in open(f'{wd}/cds_list.txt') if l.strip()]
    Dseq = {}
    if os.path.exists(f'{wd}/skani.tsv'):
        Dseq['ani'], nmiss = read_skani(f'{wd}/skani.tsv', files)
        log(f'{pid}: ANI from skani, {nmiss} missing pairs filled')
    else:
        log(f'{pid}: no skani.tsv, ANI skipped')
    Dseq['aai'] = np.load(f'{wd}/aai.npy').astype(float)
    acc = (part == 'shell') | (part == 'cloud')
    Dfull = dict(Dseq); Dfull['jac'] = jaccard_dist(R[acc].T)
    diagnostics(pid, Dfull, gc, genomes, out)
    log(f'{pid}: diagnostics written; ' + ', '.join(f'{m} median d {np.median(D[np.triu_indices(n, 1)]):.4g}'
                                                  for m, D in Dfull.items()))
    meas = list(Dseq) + ['jac']
    rel_seq = [f'{k}_{m}' for m in Dseq for k in ('ker', 'tree')]
    rel_all = rel_seq + ['ker_jac', 'tree_jac']
    preds_order = ['prev'] + rel_all + ['relseq', 'relall', 'lr', 'lr_relseq', 'lr_relall', 'avg_relseq', 'avg_relall']
    rng = np.random.default_rng(a.seed + 7)
    sizes = [s for s in SIZES if s < n] + ([n] if n <= max(SIZES) else [])
    summ, trows, grows = [], [], []
    for N in sizes:
        reps = 1 if N == n else 3
        for rep in range(reps):
            t1 = time.time()
            S = np.sort(rng.choice(n, N, replace=False)) if N < n else np.arange(n)
            RS = R[:, S]; k = RS.sum(1)
            fidx = np.flatnonzero(acc & (k >= 2) & (k <= N - 1))
            MS = M[:, S] if len(M) else M; mk = MS.sum(1) if len(M) else np.array([])
            tm = [i for i in range(len(mod_mem)) if 2 <= mk[i] <= N - 2]
            lo = max(2, N // 20)
            tf = [i for i in fidx if lo <= k[i] <= N - lo]
            if len(tm) > MAX_TGT: tm = list(rng.choice(tm, MAX_TGT, replace=False))
            if len(tf) > MAX_TGT: tf = list(rng.choice(tf, MAX_TGT, replace=False))
            tgts = [('module', i) for i in tm] + [('family', i) for i in tf]
            if not tgts: continue
            Y = np.array([MS[i] if kd == 'module' else RS[i] for kd, i in tgts], float).T      # N x T
            own = []
            for kd, i in tgts:
                s_ = set(mod_mem[i]) if kd == 'module' else ({i} | set(mod_mem[fam2mod[i]] if i in fam2mod else []))
                own.append(s_)
            allown = set().union(*own)
            jkeep = np.array([f not in allown for f in fidx])
            D = {m: Dseq[m][np.ix_(S, S)] for m in Dseq}
            D['jac'] = jaccard_dist(RS[fidx[jkeep]].T)
            trees = {m: upgma(D[m]) for m in meas}
            pcs = {m: pcoa(D[m], k=min(10, max(2, N // 10))) for m in meas}
            PCseq = np.hstack([pcs[m] for m in Dseq]); PCall = np.hstack([pcs[m] for m in meas])
            F = RS[fidx].T.astype(np.float32)
            c = comp[S]; cz = ((c - c.mean()) / (c.std() + 1e-9)).astype(np.float32)[:, None]
            fpos = {f: j for j, f in enumerate(fidx)}
            K = min(N, 10)
            fold = rng.permutation(np.arange(N) % K)
            T_ = len(tgts)
            P = {m: np.zeros((N, T_)) for m in preds_order}
            for fo in range(K):
                te = fold == fo; tr = ~te
                Ytr = Y[tr]; p0 = (Ytr.sum(0) + 1) / (tr.sum() + 2)
                P['prev'][te] = p0
                inner = {}
                for m in meas:
                    kp = kernel_predict(D[m], Y, tr, te, p0)
                    tp, _ = tree_predict(trees[m], Y, tr, c)
                    P[f'ker_{m}'][te] = kp[te]; P[f'tree_{m}'][te] = tp[te]
                    inner[f'ker_{m}'] = kp; inner[f'tree_{m}'] = tp
                for t in range(T_):
                    ytr = Ytr[:, t]
                    if ytr.min() == ytr.max():
                        for m in ('relseq', 'relall', 'lr', 'lr_relseq', 'lr_relall'): P[m][te, t] = p0[t]
                        continue
                    Zs = np.column_stack([logit(inner[r][:, t]) for r in rel_seq])
                    Za = np.column_stack([logit(inner[r][:, t]) for r in rel_all])
                    for nm, Zm in (('relseq', Zs), ('relall', Za)):
                        P[nm][te, t] = newton_predict(newton_lr(Zm[tr], ytr), Zm[te])
                    keep = np.ones(len(fidx), bool)
                    keep[[fpos[f] for f in own[t] if f in fpos]] = False
                    X = np.hstack([F[:, keep], cz])
                    mdl = LogisticRegression(C=0.1, max_iter=300).fit(X[tr], ytr)
                    P['lr'][te, t] = mdl.predict_proba(X[te])[:, 1]
                    # context + relatedness: the same penalised logistic plus the top principal coordinates of
                    # the relatedness distances (no outcome information, so train and test rows are symmetric),
                    # standardised and penalised like the family features
                    for nm, PC in (('relseq', PCseq), ('relall', PCall)):
                        X2 = np.hstack([X, PC])
                        mdl = LogisticRegression(C=0.1, max_iter=500).fit(X2[tr], ytr)
                        P[f'lr_{nm}'][te, t] = mdl.predict_proba(X2[te])[:, 1]
            for nm in ('relseq', 'relall'):                     # untuned 50/50 average of context and relatedness
                P[f'avg_{nm}'] = 0.5 * P['lr'] + 0.5 * P[nm]
            B = {m: bits(P[m], Y) for m in preds_order}              # N x T
            for t, (kd, i) in enumerate(tgts):
                trows.append(dict(id=pid, N=N, rep=rep, kind=kd, target=str(mod_ids[i]) if kd == 'module' else fams[i],
                                  carriers=int(Y[:, t].sum()),
                                  **{f'bits_{m}': round(B[m][:, t].mean(), 5) for m in preds_order}))
            kinds = np.array([kd for kd, _ in tgts])
            for gi, g in enumerate(S):
                grows.append(dict(id=pid, N=N, rep=rep, genome=genomes[g],
                                  **{f'bits_{m}': round(B[m][gi].mean(), 5) for m in preds_order}))
            dt = time.time() - t1
            for kd in ('module', 'family'):
                sel = kinds == kd
                row = dict(id=pid, n_total=n_total, pool=n, N=N, rep=rep, kind=kd, folds=K, n_targets=int(sel.sum()),
                           n_feat=len(fidx), n_feat_jac=int(jkeep.sum()), sec=round(dt, 1))
                if sel.any():
                    for m in preds_order: row[f'bits_{m}'] = round(B[m][:, sel].mean(), 5)
                summ.append(row)
            log(f'N={N} rep={rep}: {len(tm)} modules, {len(tf)} families, {len(fidx)} features, {dt:.0f}s; '
                + ' '.join(f'{m}={B[m].mean():.3f}' for m in preds_order))
            os.makedirs(f'{out}/res', exist_ok=True)
            pd.DataFrame(trows).to_csv(f'{out}/res/{pid}.targets.tsv.gz', sep='\t', index=False)
            pd.DataFrame(grows).to_csv(f'{out}/res/{pid}.genomes.tsv.gz', sep='\t', index=False)
            pd.DataFrame(summ).to_csv(f'{out}/res/{pid}.summary.tsv', sep='\t', index=False)
    open(f'{out}/res/{pid}.done', 'w').close()
    if not a.keep: shutil.rmtree(wd, ignore_errors=True)


if a.step == 'prep':
    prep(a.args[0], a.args[1], a.args[2])
else:
    os.makedirs(f'{a.args[1]}/res', exist_ok=True)
    model(a.args[0], a.args[1])
import resource
log(f"done; max RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6:.2f} GB")
