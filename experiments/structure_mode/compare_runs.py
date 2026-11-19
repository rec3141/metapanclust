"""Compare Spacedust cluster outputs (e.g. MMseqs2 vs Foldseek/ProstT5 search mode) on the same setDB.
usage: compare_runs.py SETDB_LOOKUP NAMES(comma) label1=clusters1.tsv label2=clusters2.tsv [...]
Reports clusters, hits, cluster sizes, genes covered per genome, per genome-pair cluster counts, and
the overlap of covered genes and of conserved gene pairs (two genes in one cluster) between runs."""
import sys, itertools
from collections import Counter, defaultdict
lookup, names = sys.argv[1], sys.argv[2].split(',')
runs = dict(a.split('=', 1) for a in sys.argv[3:])
pid_set = {}; ntot = Counter()
for l in open(lookup):
    k, pid, s = l.split('\t'); pid_set[pid] = int(s); ntot[int(s)] += 1
def load(path):
    clusters = []; cur = None
    for l in open(path):
        if l.startswith('#'):
            cur = {'q': set(), 't': set()}; clusters.append(cur)
        elif l.startswith('>'):
            q, t = l[1:].split('\t')[:2]; cur['q'].add(q); cur['t'].add(t)
    return clusters
res = {lab: load(p) for lab, p in runs.items()}
cov, pairs = {}, {}
for lab, cl in res.items():
    genes = set(); gp = set(); pc = Counter(); sizes = []
    for c in cl:
        genes |= c['q'] | c['t']; sizes.append(len(c['q']))
        sq = pid_set[next(iter(c['q']))]; st = pid_set[next(iter(c['t']))]
        pc[tuple(sorted((names[sq], names[st])))] += 1
        for side in (c['q'], c['t']):
            for a, b in itertools.combinations(sorted(side), 2): gp.add((a, b))
    cov[lab], pairs[lab] = genes, gp
    per = Counter(pid_set[g] for g in genes)
    sizes.sort()
    print(f'== {lab}: {len(cl)} clusters, {sum(sizes)} hits (query genes); median genes/cluster side {sizes[len(sizes)//2] if sizes else 0}, max {max(sizes) if sizes else 0}')
    print(f'   genes covered {len(genes)} of {sum(ntot.values())} ({len(genes)/sum(ntot.values()):.1%}); co-clustered gene pairs {len(gp)}')
    print('   per genome: ' + '  '.join(f'{names[s]} {per[s]}/{ntot[s]} ({per[s]/ntot[s]:.0%})' for s in sorted(ntot)))
    print('   per genome pair: ' + '  '.join(f'{a}-{b} {n}' for (a, b), n in sorted(pc.items())))
for a, b in itertools.combinations(res, 2):
    ga, gb = cov[a], cov[b]; pa, pb = pairs[a], pairs[b]
    print(f'== {a} vs {b}: genes both {len(ga & gb)}, only {a} {len(ga - gb)}, only {b} {len(gb - ga)}; '
          f'gene pairs both {len(pa & pb)}, only {a} {len(pa - pb)}, only {b} {len(pb - pa)}')
