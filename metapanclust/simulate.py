"""Held-out evaluation at gene level: leave a reference genome out, fragment it, drop pieces to a target
completeness, optionally add pieces of a foreign genome (another species, genes mapped to this reference by
MMseqs2) or of another genome of the same species (a strain merge, which should not look surprising), and
score it against the remaining genomes.
"""
import os, subprocess
import numpy as np, pandas as pd
from .ref import Ref
from .score import score


class Genes:
    def __init__(self, refdir):
        z = np.load(f'{refdir}/genes.npz')
        self.g = z['g']; self.contig = z['contig']; self.start = z['start']; self.stop = z['stop']; self.f = z['f']
        self.ctg_len = dict(zip(z['ctg_id'].tolist(), z['ctg_len'].tolist()))
        self.idx = {}
        order = np.argsort(self.g, kind='stable')
        bounds = np.searchsorted(self.g[order], np.arange(self.g.max() + 2))
        for j in range(self.g.max() + 1): self.idx[j] = order[bounds[j]:bounds[j + 1]]

    def genome(self, j):
        i = self.idx[j]
        return pd.DataFrame({'contig': self.contig[i], 'start': self.start[i], 'stop': self.stop[i], 'f': self.f[i]})


def fragment(gdf, ctg_len, rng, median=15000, sigma=1.0, prefix=''):
    """cut each contig into lognormal pieces; genes crossing a break are lost. Returns (genes, piece lengths)"""
    rows, pieces = [], {}
    for ctg, d in gdf.groupby('contig'):
        Lc = ctg_len.get(ctg, int(d.stop.max()) + 100)
        cuts = [0]
        while cuts[-1] < Lc:
            cuts.append(cuts[-1] + max(1000, int(rng.lognormal(np.log(median), sigma))))
        cuts[-1] = Lc
        cuts = np.array(cuts)
        for k in range(len(cuts) - 1):
            pieces[f'{prefix}{ctg}_{k}'] = int(cuts[k + 1] - cuts[k])
        a = np.searchsorted(cuts, d.start.values, side='right') - 1
        b = np.searchsorted(cuts, d.stop.values, side='left') - 1
        ok = a == b
        rows.append(d[ok].assign(contig=[f'{prefix}{ctg}_{k}' for k in a[ok]]))
    return pd.concat(rows), pieces


def drop_to(genes, pieces, frac, rng):
    """keep a random subset of pieces totalling about frac of the length"""
    names = list(pieces); rng.shuffle(names)
    tot = sum(pieces.values()); keep, s = [], 0
    for n in names:
        if s >= frac * tot: break
        keep.append(n); s += pieces[n]
    ks = set(keep)
    return genes[genes.contig.isin(ks)], {n: pieces[n] for n in keep}


def foreign_map(refdir, frefdir, work, threads=4):
    """foreign family -> (this reference's family, pident), from the foreign families' representative proteins"""
    from .annotate import map_to_ref
    os.makedirs(work, exist_ok=True)
    out = f'{work}/frep.faa'
    seen = set()
    with open(f'{frefdir}/ref_prot.faa') as fi, open(out, 'w') as fo:
        h = None
        for line in fi:
            if line.startswith('>'):
                f = int(line[1:].split('|')[1]); h = f if f not in seen else None
                if h is not None: seen.add(f); fo.write(f'>{f}\n')
            elif h is not None:
                fo.write(line)
    m = map_to_ref(out, refdir, work, threads)
    return {int(q): (int(r.fam), float(r.pident)) for q, r in m.iterrows()}


def run(refdir, frefdir, outdir, n=20, comps=(0.7, 0.9), contam=(0.05, 0.1, 0.2), seed=1, threads=4,
        use_context=True):
    rng = np.random.default_rng(seed)
    ref = Ref(refdir); G = Genes(refdir)
    fref = Ref(frefdir); FG = Genes(frefdir)
    fmap = foreign_map(refdir, frefdir, f'{outdir}/foreign_map', threads)
    N = len(ref.genomes)
    cand = np.flatnonzero(ref.comp >= 0.95)
    held = rng.choice(cand, min(n, len(cand)), replace=False)
    rows, crow = [], []
    for j in held:
        keep = np.ones(N, bool); keep[j] = False
        r = ref.drop(keep)
        alive = r.R.sum(1) > 0
        base = G.genome(j)
        for c in comps:
            gq, pq = fragment(base, G.ctg_len, rng)
            gq, pq = drop_to(gq, pq, c, rng)
            native = gq.assign(fam=np.where(alive[gq.f.values], gq.f.values, -1), pident=100.0, src='native')
            conds = [('clean', 0.0, None)] + [('foreign', x, 'f') for x in contam] + [('strain', 0.1, 's')]
            for kind, x, tag in conds:
                genes, lens = native, dict(pq)
                if tag:
                    if tag == 'f':
                        k = rng.integers(len(fref.genomes)); src = FG.genome(k); clen = FG.ctg_len
                    else:
                        k = rng.choice(np.flatnonzero(keep)); src = G.genome(k); clen = G.ctg_len
                    fg, fp = fragment(src, clen, rng, prefix=f'{tag}{k}_')
                    want = x * sum(pq.values()) / (1 - x)
                    fg, fp = drop_to(fg, fp, want / sum(fp.values()), rng)
                    if tag == 'f':
                        mp = [fmap.get(int(f), (-1, np.nan)) for f in fg.f.values]
                        fg = fg.assign(fam=[m[0] for m in mp], pident=[m[1] for m in mp], src='foreign')
                    else:
                        fg = fg.assign(fam=np.where(alive[fg.f.values], fg.f.values, -1), pident=100.0,
                                       src='strain')
                    genes = pd.concat([native, fg]); lens.update(fp)
                s, U, C = score(r, genes, lens.values(), comp=None, use_context=use_context, seed=int(j))
                s.update(held=int(j), genome=ref.genomes[j], target_comp=c, cond=kind, contam=x,
                         true_contam_Mb=round(sum(v for k_, v in lens.items() if not k_ in pq) / 1e6, 3))
                rows.append(s)
                src = genes.groupby('contig').src.first()
                C['src'] = C.contig.map(src); C['held'] = int(j); C['cond'] = kind; C['contam'] = x
                C['target_comp'] = c; C['len'] = C.contig.map(lens)
                crow.append(C)
                if kind == 'clean':
                    U.assign(held=int(j), target_comp=c).to_csv(f'{outdir}/units_{j}_{c}.tsv.gz', sep='\t',
                                                                index=False)
        pd.DataFrame(rows).to_csv(f'{outdir}/bins.tsv', sep='\t', index=False)
        pd.concat(crow).to_csv(f'{outdir}/contigs.tsv.gz', sep='\t', index=False)
        print(f'held-out genome {j} done ({len(rows)} bins)', flush=True)
    return pd.DataFrame(rows)
