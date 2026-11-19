"""(1) Fragment the six complete reference genomes with contig lengths drawn from each Chesterfield
bin, drop pieces to the bin's CheckM2 completeness, and score curated E. coli operons and six-genome
Spacedust cluster instances as intact / partial / lost. Writes one row per cluster x bin x replicate
summary (aggregated) to OUT/sim_clusters.tsv.gz.
usage: simulate.py [REPS]"""
import sys, gzip
import numpy as np
import common as C

REPS = int(sys.argv[1]) if len(sys.argv) > 1 else 10
rng = np.random.default_rng(1)

# --- reference replicons, CDS and clusters -----------------------------------------------------
reps = {}  # replicon -> dict(G, starts, ends, sets: {name: (member_idx flat, cluster_id flat, k, span)})
clusters = {'operon': C.eco_operons(), 'sd6': C.sd_instances(f'{C.D}/bench/sd6_mmseqs/clusters.tsv')}
for g, fa in C.REFS.items():
    lens = C.fasta_lengths(f'{C.BENCH}/{g}/{fa}.fasta')
    cds, _ = C.gff_cds(f'{C.BENCH}/{g}/mag_out/annotation/bakta/basic/annotation.gff3')
    for r, G in lens.items():
        if r in cds:
            a = cds[r]
            reps[r] = dict(genome=g, G=G, s=a[:, 0], e=a[:, 1],
                           idx={(int(x), int(y)): i for i, (x, y) in enumerate(a)})
meta = []  # (set, replicon, local cluster id, k, span)
for name, cl in clusters.items():
    for r in reps:
        reps[r].setdefault('sets', {})
    per = {}
    for c, mem in cl:
        R = reps.get(c)
        if R is None:
            continue
        ix = [R['idx'][m] for m in mem if m in R['idx']]
        if len(ix) < 2:
            continue
        per.setdefault(c, []).append(sorted(ix))
    for c, lst in per.items():
        R = reps[c]
        flat = np.concatenate(lst); cid = np.repeat(np.arange(len(lst)), [len(x) for x in lst])
        k = np.array([len(x) for x in lst])
        span = np.array([R['e'][x].max() - R['s'][x].min() + 1 for x in lst])
        R['sets'][name] = (flat, cid, k, span, len(meta))
        meta += [(name, c, i, k[i], span[i]) for i in range(len(lst))]
print(f'{sum(1 for m in meta if m[0]=="operon")} operons, {sum(1 for m in meta if m[0]=="sd6")} sd6 instances', file=sys.stderr)

q = C.checkm2()
bins = C.bin_names()
out = gzip.open(f'{C.OUT}/sim_clusters.tsv.gz', 'wt')
out.write('bin\tN50\tcompleteness\tset\treplicon\tk\tspan\tfrag_intact\tintact\tdetected\tlost\tend_touch\n')
for b in bins:
    L = np.array(list(C.fasta_lengths(f'{C.BINS}/{b}.fna').values()))
    comp = min(q[b][0], 1.0)
    N50 = C.n50(L)
    acc = {}
    for r, R in reps.items():
        G = R['G']
        for rep in range(REPS):
            # stationary laydown: first piece length-biased, random phase
            L1 = rng.choice(L, p=L / L.sum())
            x0 = -rng.uniform(0, L1)
            n_need = int(G / L.mean() * 1.5) + 10
            while True:
                pieces = np.concatenate([[L1], rng.choice(L, n_need)])
                cuts = x0 + np.concatenate([[0], np.cumsum(pieces)])
                if cuts[-1] > G:
                    break
                n_need *= 2
            cuts = cuts[: np.searchsorted(cuts, G) + 1]
            cuts[0], cuts[-1] = 0, G + 1  # positions are 1-based inclusive; pieces [cuts[j], cuts[j+1])
            npc = len(cuts) - 1
            keep = rng.random(npc) < comp
            ps = np.searchsorted(cuts, R['s'], side='right') - 1
            pe = np.searchsorted(cuts, R['e'], side='right') - 1
            whole = ps == pe
            obs = whole & keep[ps]
            # first/last observed (whole) CDS per piece
            wi = np.flatnonzero(whole)
            first = np.full(npc, -1); last = np.full(npc, -1)
            first[ps[wi][::-1]] = wi[::-1]; last[ps[wi]] = wi
            for name, (flat, cid, k, span, off) in R.get('sets', {}).items():
                nc = len(k)
                p = ps[flat]
                fi = np.bincount(cid, weights=whole[flat], minlength=nc) == k
                same = np.ones(nc, bool)
                pmin = np.full(nc, np.inf); pmax = np.full(nc, -np.inf)
                np.minimum.at(pmin, cid, p); np.maximum.at(pmax, cid, p)
                same = pmin == pmax
                frag_intact = fi & same
                nobs = np.bincount(cid, weights=obs[flat], minlength=nc)
                intact = frag_intact & (nobs == k)
                lost = nobs == 0
                # largest observed group of members on one piece
                o = obs[flat]
                key = cid[o].astype(np.int64) * (npc + 1) + p[o]
                uk, cnt = np.unique(key, return_counts=True)
                best = np.zeros(nc, int); bestp = np.full(nc, -1)
                order = np.argsort(cnt, kind='stable')
                best[uk[order] // (npc + 1)] = cnt[order]; bestp[uk[order] // (npc + 1)] = uk[order] % (npc + 1)
                detected = best >= 2
                # end touch for the detected group
                sel = o & (p == bestp[cid])
                mi = np.full(nc, 1 << 30); ma = np.full(nc, -1)
                np.minimum.at(mi, cid[sel], flat[sel]); np.maximum.at(ma, cid[sel], flat[sel])
                bp = np.maximum(bestp, 0)
                touch = detected & ((mi == first[bp]) | (ma == last[bp]))
                a = acc.setdefault((name, r), np.zeros((nc, 5)))
                a += np.stack([frag_intact, intact, detected, lost, touch], 1)
    for (name, r), a in acc.items():
        flat, cid, k, span, off = reps[r]['sets'][name]
        for i in range(len(k)):
            out.write(f'{b}\t{N50}\t{comp:.4f}\t{name}\t{r}\t{k[i]}\t{span[i]}\t' + '\t'.join(f'{v / REPS:.2f}' for v in a[i]) + '\n')
    print(b, N50, comp, file=sys.stderr)
out.close()
