"""Spacedust pairwise clusters -> per-bin conserved loci, split at weakly supported gene-gene links.

A link joins consecutive genes i, i+1 on a contig (Spacedust indices count CDSs across the whole genome). A pairwise cluster (bin q vs bin t) supports every link
between two of its q-side members that sit at most MAX_GAP genes apart (and likewise on the t side).
Link support counts distinct partner *genera* other than the bin's own genus: synteny shared with a close
relative spans whole contigs and carries no information about functional units, so it cannot hold a
locus together on its own.
A locus is a maximal run of links with support >= MIN_SUP, cut again where support dips below
DIP * (smaller of the maxima on either side) -- two conserved modules that abut and are chained only by
a few partners become two loci.
usage: 02_loci.py CLUSTERS_TSV GENES_TSV BIN_TAXONOMY OUT_PREFIX [min_sup=1] [dip=0.34] [max_gap=2]
Genomes are taken from GENES_TSV (contig -> bin), not from the Spacedust file names."""
import sys
from collections import defaultdict
import pandas as pd

sd, genes_tsv, tax_tsv, out = sys.argv[1:5]
MIN_SUP = int(sys.argv[5]) if len(sys.argv) > 5 else 1
DIP = float(sys.argv[6]) if len(sys.argv) > 6 else 0.34
MAX_GAP = int(sys.argv[7]) if len(sys.argv) > 7 else 2

genes = pd.read_csv(genes_tsv, sep='\t', dtype={'idx': int}, keep_default_na=False)
contig_bin = dict(zip(genes.contig, genes.bin))
lo_idx = genes.groupby('contig').idx.min().to_dict()   # Spacedust indices run across contigs
ngenes = genes.groupby('contig').idx.max().to_dict()   # last index on the contig
genus = {}
for l in open(tax_tsv):
    b, t = l.rstrip('\n').split('\t')[:2]
    r = t.split(';'); genus[b] = r[5] if len(r) > 5 and r[5] != 'g__' else b   # unnamed genus -> own group

def gid(p):
    c, i, _, _ = p.rsplit('_', 3); return c, int(i)

sup = defaultdict(set)            # (contig, i) link i->i+1 -> partner bins
def add_side(members, partner):
    byc = defaultdict(list)
    for c, i in members: byc[c].append(i)
    for c, idx in byc.items():
        idx = sorted(set(idx))
        for a, z in zip(idx, idx[1:]):
            if z - a <= MAX_GAP + 1:
                for k in range(a, z): sup[(c, k)].add(partner)

cur = None
def flush():
    if cur: add_side(cur[1], cur[2]); add_side(cur[3], cur[0])
for l in open(sd):
    if l[0] == '#':
        flush()
        f = l[1:].split('\t'); cur = [None, [], None, []]
    elif l[0] == '>':
        q, t = l[1:].split('\t', 2)[:2]
        q, t = gid(q), gid(t); cur[1].append(q); cur[3].append(t)
        if cur[0] is None: cur[0], cur[2] = contig_bin[q[0]], contig_bin[t[0]]
flush()

def score(c, i):
    b = contig_bin[c]; g = genus.get(b, b)
    return len({genus.get(p, p) for p in sup.get((c, i), ())} - {g})

link_rows, loci = [], []
for c in sorted({c for c, _ in sup}):
    o = lo_idx[c]; n = ngenes[c] - o + 1             # local position p <-> Spacedust index o + p
    s = [score(c, o + i) for i in range(n - 1)]
    nb = [len(sup.get((c, o + i), ())) for i in range(n - 1)]
    for i in range(n - 1):
        if nb[i]: link_rows.append((c, o + i, s[i], nb[i]))
    i = 0
    while i < n - 1:                                  # runs of links with s >= MIN_SUP
        if s[i] < MIN_SUP: i += 1; continue
        j = i
        while j + 1 < n - 1 and s[j + 1] >= MIN_SUP: j += 1
        segs = [(i, j)]; done = []
        while segs:                                   # recursive split at the deepest dip
            a, z = segs.pop()
            best = None
            for k in range(a + 1, z):
                lm, rm = max(s[a:k]), max(s[k + 1:z + 1])
                if s[k] < DIP * min(lm, rm) and (best is None or s[k] / min(lm, rm) < best[0]):
                    best = (s[k] / min(lm, rm), k)
            if best: k = best[1]; segs += [(a, k - 1), (k + 1, z)]
            else: done.append((a, z))
        for a, z in sorted(done):
            loci.append((contig_bin[c], c, o + a, o + z + 1, min(s[a:z + 1]), max(s[a:z + 1])))
        i = j + 1

pd.DataFrame(link_rows, columns=['contig', 'idx', 'support_genera', 'support_bins']).to_csv(f'{out}_links.tsv', sep='\t', index=False)
L = pd.DataFrame(loci, columns=['bin', 'contig', 'first', 'last', 'min_sup', 'max_sup'])
L.insert(0, 'locus', [f'L{k:06d}' for k in range(len(L))])
L['n_genes'] = L['last'] - L['first'] + 1
L['at_contig_end'] = (L['first'] == L.contig.map(lo_idx)) | (L['last'] == L.contig.map(ngenes))
L.to_csv(f'{out}_loci.tsv', sep='\t', index=False)
q = L.n_genes.quantile([.5, .9, .99]).tolist()
cov = L.n_genes.sum() / len(genes)
print(f'{len(link_rows)} supported links; {len(L)} loci; genes/locus median {q[0]:.0f}, 90% {q[1]:.0f}, 99% {q[2]:.0f}, max {L.n_genes.max()}; '
      f'{cov:.1%} of genes in a locus', file=sys.stderr)
