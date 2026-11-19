"""Score E. coli operon-pair recovery for reference panels built from pairwise Spacedust runs.
usage: score_panels.py PANEL_TSV RUN_DIR PAIRS_TSV BAKTA_TSV OUT_PREFIX

A pair is 'conserved' against genome g when both E. coli genes sit in one Spacedust cluster of the
E. coli-vs-g run (the eval_pairs.py definition). 'adjacent' is stricter: both genes hit genes of g that
are neighbours (index differs by 1 on one contig) in the same relative orientation. A panel calls a pair
when at least k of its genomes do."""
import csv, itertools, random, sys
from collections import defaultdict

panel_tsv, run_dir, pairs_tsv, bakta_tsv, out = sys.argv[1:6]
panel = list(csv.DictReader(open(panel_tsv), delimiter='\t'))
names = [p['name'] for p in panel]
level = {p['name']: int(p['level']) for p in panel}
group = {int(p['level']): p['group'] for p in panel}

loc = {}
for l in open(bakta_tsv):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if len(f) > 5 and f[1] == 'cds': loc[(f[0], int(f[2]), int(f[3]))] = f[5]

def parse(pid):
    c, i, a, b = pid.rsplit('_', 3)
    a, b = int(a), int(b)
    return c, int(i), min(a, b), max(a, b), (1 if b >= a else -1)

rows = [r for r in csv.DictReader(open(pairs_tsv), delimiter='\t') if r['same_tu'] != '']
pair_keys = [(r['locus_a'], r['locus_b']) for r in rows]
lab = [r['same_tu'] == '1' for r in rows]
gap = [int(r['gap']) for r in rows]

# per genome: set of pair indices conserved (co-membership) and adjacent-in-target
cons, adj = {}, {}
for n in names:
    member = defaultdict(set); hit = defaultdict(set)   # locus -> cluster ids ; locus -> (cid, contig, idx, strand)
    cid = None
    try: fh = open(f'{run_dir}/{n}/clusters.tsv')
    except FileNotFoundError: print('missing', n, file=sys.stderr); continue
    for l in fh:
        if l[0] == '#': cid = l[1:].split('\t')[0]; continue
        if l[0] != '>': continue
        q, t = l[1:].split('\t')[:2]
        qc, qi, qa, qb, qs = parse(q); tc, ti, ta, tb, ts = parse(t)
        g = loc.get((qc, qa, qb))
        if g is None: continue
        member[g].add(cid); hit[g].add((cid, tc, ti, qs * ts))
    cs, ad = set(), set()
    for k, (a, b) in enumerate(pair_keys):
        if member[a] & member[b]:
            cs.add(k)
            if any(x[0] == y[0] and x[1] == y[1] and abs(x[2] - y[2]) == 1 and x[3] == y[3]
                   for x in hit[a] for y in hit[b]): ad.add(k)
    cons[n], adj[n] = cs, ad
names = [n for n in names if n in cons]

def metrics(called):
    tp = sum(1 for k in called if lab[k]); fp = len(called) - tp; P = sum(lab)
    comb = {k for k in range(len(rows)) if gap[k] <= 50 or (k in called and gap[k] <= 300)}
    ctp = sum(1 for k in comb if lab[k]); cfp = len(comb) - ctp
    acc = (ctp + (len(rows) - P - cfp)) / len(rows)
    far = [k for k in range(len(rows)) if gap[k] > 50]
    return dict(called=len(called), prec=tp / max(len(called), 1), rec=tp / P,
                far_tp=sum(1 for k in far if k in called and lab[k]), far_fp=sum(1 for k in far if k in called and not lab[k] and gap[k] <= 300),
                near_diff_cons=sum(1 for k in called if gap[k] <= 50 and not lab[k]),
                c_acc=acc, c_prec=ctp / len(comb), c_rec=ctp / P)

def panel_call(genomes, k=1, which=cons):
    cnt = defaultdict(int)
    for g in [g for g in genomes if g in which]:
        for p in which[g]: cnt[p] += 1
    return {p for p, c in cnt.items() if c >= k}

hdr = ['panel', 'n', 'k', 'mode', 'called', 'prec', 'rec', 'far_tp', 'far_fp', 'near_diff_cons', 'c_acc', 'c_prec', 'c_rec']
res = []
def add(label, genomes, k=1, mode='cons'):
    m = metrics(panel_call(genomes, k, cons if mode == 'cons' else adj))
    res.append([label, len([g for g in genomes if g in cons]), k, mode] + [m[h] for h in hdr[4:]])

P = sum(lab); far_same = sum(1 for k in range(len(rows)) if gap[k] > 50 and lab[k])
far_diff300 = sum(1 for k in range(len(rows)) if 50 < gap[k] <= 300 and not lab[k])
print(f'{len(rows)} pairs, {P} same TU; gap>50 same-TU {far_same}; 50<gap<=300 different-TU {far_diff300}')
d = metrics(set()); print(f'distance rule alone: acc {d["c_acc"]:.3f} prec {d["c_prec"]:.3f} rec {d["c_rec"]:.3f}')

for n in names: add(f'single:{n}', [n]); add(f'single:{n}', [n], mode='adj')
levels = sorted(set(level.values()))
for L in levels:
    gs = [n for n in names if level[n] == L]
    for mode in ('cons', 'adj'):
        add(f'level{L}', gs, mode=mode)
        if len(gs) >= 3: add(f'level{L}', gs, k=2, mode=mode)
for L in levels:   # nested: everything at distance >= L, and everything <= L
    for mode in ('cons', 'adj'):
        add(f'ge_level{L}', [n for n in names if level[n] >= L], mode=mode)
        add(f'le_level{L}', [n for n in names if level[n] <= L], mode=mode)
add('baseline5', ['Bsu', 'Pae', 'Mtu', 'Syn', 'Mja'])
add('all', names); add('all', names, mode='adj')
for k in (2, 3, 5, 8):
    add('all', names, k=k); add('ge_level2', [n for n in names if level[n] >= 2], k=k)

# number of genomes: random panels of size n from pools; mean over draws
rng = random.Random(1)
curve = []
pools = {'L0-1 (Escherichia)': [n for n in names if level[n] <= 1],
         'L2-3 (Enterobacterales)': [n for n in names if level[n] in (2, 3)],
         'L4-6 (other Gamma+Beta)': [n for n in names if level[n] in (4, 5, 6)],
         'L7-9 (Alpha, other phyla, Archaea)': [n for n in names if level[n] >= 7],
         'L2-9 (all but Escherichia)': [n for n in names if level[n] >= 2],
         'L0-9 (all)': names,
         'one per level 2-8 stratified': None}
for pool, gs in pools.items():
    for size in (1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40, 50, 60, 68):
        draws = []
        for _ in range(30):
            if gs is None:
                if size > 7: break
                lv = rng.sample([2, 3, 4, 5, 6, 7, 8], size)
                pick = [rng.choice([n for n in names if level[n] == L]) for L in lv]
            else:
                if size > len(gs): break
                pick = rng.sample(gs, size)
            draws.append(metrics(panel_call(pick)))
        if not draws: continue
        mean = {h: sum(d[h] for d in draws) / len(draws) for h in ('prec', 'rec', 'c_acc', 'c_prec', 'c_rec', 'far_tp', 'far_fp')}
        curve.append([pool, size] + [mean[h] for h in ('prec', 'rec', 'c_acc', 'c_prec', 'c_rec', 'far_tp', 'far_fp')])

with open(out + '_panels.tsv', 'w') as f:
    f.write('\t'.join(hdr) + '\n')
    for r in res: f.write('\t'.join(f'{x:.3f}' if isinstance(x, float) else str(x) for x in r) + '\n')
with open(out + '_curve.tsv', 'w') as f:
    f.write('pool\tsize\tprec\trec\tc_acc\tc_prec\tc_rec\tfar_tp\tfar_fp\n')
    for r in curve: f.write('\t'.join(f'{x:.3f}' if isinstance(x, float) else str(x) for x in r) + '\n')
print('wrote', out + '_panels.tsv', out + '_curve.tsv')

# per-pair matrix for model fitting: gap, label, then conserved (1) / adjacent (2) / neither (0) per genome
with open(out + '_matrix.tsv', 'w') as f:
    f.write('\t'.join(['locus_a', 'locus_b', 'gap', 'same_tu'] + names) + '\n')
    for k, (a, b) in enumerate(pair_keys):
        f.write('\t'.join([a, b, str(gap[k]), str(int(lab[k]))] + [str(2 if k in adj[n] else 1 if k in cons[n] else 0) for n in names]) + '\n')
print('wrote', out + '_matrix.tsv')
