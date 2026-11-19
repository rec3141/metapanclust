"""Build a small benchmark from an existing all-vs-all Spacedust run: take the clusters between two genomes,
expand each cluster's genes by +-FLANK neighbours, and write per-genome GFF3 files with only those CDS.
Structure-mode and MMseqs2-mode Spacedust are then rerun on the identical subset (fair comparison of
cluster extension within known conserved neighbourhoods; it cannot find wholly new neighbourhoods).
usage: subset_regions.py SETDB_LOOKUP CLUSTERS_TSV SET_A SET_B FLANK GFF_A GFF_B OUTDIR"""
import sys, os
lookup, clu, sa, sb, flank, gffa, gffb, out = sys.argv[1:9]
sa, sb, flank = int(sa), int(sb), int(flank)
genes = {}; order = {}                         # pid -> (set, contig, idx); per set ordered pid list per contig
for l in open(lookup):
    k, pid, s = l.rstrip('\n').split('\t'); s = int(s)
    acc, idx, a, b = pid.rsplit('_', 3)
    genes[pid] = (s, acc, int(idx)); order.setdefault((s, acc), {})[int(idx)] = pid
keep = set(); cur_ok = False
for l in open(clu):
    if l.startswith('#'): cur_ok = None; continue
    q, t = l[1:].split('\t')[:2]
    if {genes[q][0], genes[t][0]} != {sa, sb}: continue
    for p in (q, t):
        s, acc, i = genes[p]
        for j in range(i - flank, i + flank + 1):
            if j in order[(s, acc)]: keep.add(order[(s, acc)][j])
os.makedirs(out, exist_ok=True)
for s, gff in ((sa, gffa), (sb, gffb)):
    coords = set()
    for p in keep:
        if genes[p][0] == s:
            acc, idx, a, b = p.rsplit('_', 3); coords.add((acc, min(int(a), int(b)), max(int(a), int(b))))
    n = 0
    with open(f'{out}/set{s}.gff3', 'w') as fh:
        for l in open(gff):
            f = l.split('\t')
            if l.startswith('#') or (len(f) > 4 and f[2] == 'CDS' and (f[0], int(f[3]) - 1, int(f[4]) - 1) in coords) \
               or (len(f) > 4 and f[2] == 'CDS' and (f[0], int(f[3]), int(f[4])) in coords):
                fh.write(l); n += not l.startswith('#')
    print(f'set {s}: {len(coords)} genes selected, {n} CDS written')
