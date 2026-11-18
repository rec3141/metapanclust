"""Summarise an all-vs-all Spacedust run over bins: genes covered per bin, cluster sizes, and how
cluster sharing depends on the taxonomic relation of the two bins (GTDB-Tk lineages).
usage: sd_summary.py CLUSTERS_TSV BINS_DIR GTDBTK_TSV"""
import csv, glob, os, sys
from collections import Counter, defaultdict
sd, bdir, gt = sys.argv[1:4]
contig2bin, ncds = {}, Counter()
for g in glob.glob(f'{bdir}/*.gff3'):
    b = os.path.basename(g)[:-5]
    for l in open(g):
        if not l.startswith('#'):
            c = l.split('\t')[0]; contig2bin[c] = b; ncds[b] += 1
lin = {r['user_genome']: r['classification'] for r in csv.DictReader(open(gt), delimiter='\t')}
def rel(a, b):
    la, lb = lin.get(a, '').split(';'), lin.get(b, '').split(';')
    for i, rk in ((6, 'same species'), (5, 'same genus'), (4, 'same family'), (3, 'same order'), (2, 'same class'), (1, 'same phylum'), (0, 'same domain')):
        if len(la) > i and len(lb) > i and la[i] == lb[i] and len(la[i]) > 3:
            return rk
    return 'different domain / unknown'
covered = defaultdict(set); sizes = Counter(); pairs = Counter(); pair_clusters = Counter()
cur = None
for l in open(sd):
    if l.startswith('#'):
        cur = l; continue
    if l.startswith('>'):
        q, t = l[1:].split('\t')[:2]
        for p in (q, t):
            c = p.rsplit('_', 3)[0]
            b = contig2bin.get(c)
            if b: covered[b].add(p)
for l in open(sd):
    if l.startswith('#'):
        f = l[1:].rstrip('\n').split('\t'); sizes[min(int(f[5]), 10)] += 1
        qb, tb = f[1].replace('.gff3', ''), f[2].replace('.gff3', '')
        pair_clusters[rel(qb, tb)] += 1
bins = sorted(ncds)
frac = sorted(len(covered[b]) / ncds[b] for b in bins)
print(f'{len(bins)} bins, {sum(ncds.values())} CDS; share of each bin\'s genes in >= 1 conserved cluster: '
      f'median {frac[len(frac) // 2]:.2f} (range {frac[0]:.2f}-{frac[-1]:.2f})')
print('clusters by number of gene pairs (10 = 10+):', dict(sorted(sizes.items())))
print('clusters by relation of the two bins:', dict(pair_clusters.most_common()))
