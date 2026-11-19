"""P1 of docs/theory/phylo_prior.md on the Chesterfield 93 bins: within a deepest-shared GTDB rank, is the
closer half of bin pairs more concordant than the farther half? Distance = patristic distance in the GTDB
r226 bac120 tree between each bin's GTDB-Tk closest reference (ANI-screen reference, else placement
reference); bins with neither are dropped, as are pairs with the same reference.
Concordance = completeness-corrected Pr(z_j=1 | seen_i), as in theory_phylo/rank_decay.py."""
import numpy as np, csv, re, sys
D = '/data/scratch/metapanclust_data/catalogue_surprise/'
TAX = '/data/scratch/metapanclust_data/chesterfield/mag_chesterfield/taxonomy/gtdbtk/gtdbtk_taxonomy.tsv'
TREE = '/data/scratch/refdbs/GTDB/release226/metadata/bac120_r226.tree'
d = np.load(D + 'inputs.npz', allow_pickle=True)
obs = d['obs']; bins = list(d['bins'])
meta = {r['bin']: r for r in csv.DictReader(open(D + 'bins.tsv'), delimiter='\t')}
ref = {}
for r in csv.DictReader(open(TAX), delimiter='\t'):
    x = r['closest_genome_reference'] if r['closest_genome_reference'] != 'N/A' else r['closest_placement_reference']
    if x != 'N/A': ref[r['user_genome']] = x
# --- newick: parent pointers and branch lengths
s = open(TREE).read()
par, bl, name = [-1], [0.0], ['']
stack, cur, last = [], 0, None
for tok in re.finditer(r"\(|\)|,|;|'[^']*'|:[-0-9.eE]+|[^(),;:']+", s):
    t = tok.group()
    if t == '(':
        par.append(cur); bl.append(0.0); name.append(''); stack.append(cur); cur = len(par) - 1; last = cur
    elif t == ',':
        cur = stack[-1]; par.append(cur); bl.append(0.0); name.append(''); cur = len(par) - 1; last = cur
    elif t == ')':
        cur = stack.pop(); last = cur
    elif t == ';': break
    elif t[0] == ':': bl[last] = float(t[1:])
    else: name[last] = t
# the root '(' created node 1 under the dummy 0; leaves are named nodes created at '(' or ','
leafidx = {n[3:] if n[:3] in ('RS_', 'GB_') else n: i for i, n in enumerate(name) if n.startswith(('RS_', 'GB_'))}
par = np.array(par); bl = np.array(bl)
def anc(i):
    out = {}; dist = 0.0
    while i >= 0: out[i] = dist; dist += bl[i]; i = par[i]
    return out
use = [b for b in bins if b in ref and ref[b] in leafidx]
print(f'bins with a reference in the tree: {len(use)} / {len(bins)}', file=sys.stderr)
A = {b: anc(leafidx[ref[b]]) for b in use}
def pdist(b1, b2):
    a1, a2 = A[b1], A[b2]
    i = leafidx[ref[b2]]; dd = 0.0
    while i not in a1: dd += bl[i]; i = par[i]
    return dd + a1[i]
ranks = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
def shared(a, b):
    r = -1
    for i in range(7):
        if a[i] == b[i] and len(a[i]) > 3: r = i
        else: break
    return r
bi = {b: i for i, b in enumerate(bins)}
seen = (obs > 0).astype(float); comp = np.array([float(meta[b]['compl']) for b in bins])
nb = seen.sum(0)
strata = {'2-3': (nb >= 2) & (nb <= 3), '4-9': (nb >= 4) & (nb <= 9), '10-29': (nb >= 10) & (nb <= 29), '>=30': nb >= 30}
pairs = {}
for x in range(len(use)):
    for y in range(x + 1, len(use)):
        b1, b2 = use[x], use[y]
        if ref[b1] == ref[b2]: continue
        r = shared(meta[b1]['lineage'].split(';'), meta[b2]['lineage'].split(';'))
        pairs.setdefault(r, []).append((bi[b1], bi[b2], pdist(b1, b2)))
def cond(P, m):
    I = np.array([p[0] for p in P]); J = np.array([p[1] for p in P])
    a = seen[I][:, m]; b = seen[J][:, m]; ca = comp[I][:, None]; cb = comp[J][:, None]
    return 2 * (a * b).sum() / ((a * cb).sum() + (b * ca).sum())
print('rank\tpairs\td_min\td_median\td_max\tstratum\tcond_near\tcond_far\tcond_all\td_near_med\td_far_med')
for r in sorted(pairs):
    if r < 1: continue
    P = sorted(pairs[r], key=lambda p: p[2]); n = len(P); dd = np.array([p[2] for p in P])
    near, far = P[:n // 2], P[n - n // 2:]
    for sn, m in strata.items():
        print(f'{ranks[r]}\t{n}\t{dd.min():.3f}\t{np.median(dd):.3f}\t{dd.max():.3f}\t{sn}\t{cond(near, m):.2f}\t'
              f'{cond(far, m):.2f}\t{cond(P, m):.2f}\t{np.median(dd[:n//2]):.3f}\t{np.median(dd[n-n//2:]):.3f}')
