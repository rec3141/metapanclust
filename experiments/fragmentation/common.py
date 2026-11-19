"""Shared loaders for the fragmentation experiment: CDS tables, contig lengths, Spacedust cluster
instances, curated E. coli operons, CheckM2 completeness."""
import csv, glob, os
from collections import defaultdict
import numpy as np

D = '/data/scratch/metapanclust_data'
BENCH = '/data/dev/metacyc_benchmark'
REFS = {'ecoli': 'ecoli_k12', 'bsub': 'bsub168', 'pae': 'pae', 'mtb': 'mtb', 'syn': 'syn', 'mja': 'mja'}
BINS = f'{D}/chesterfield/bins'
OUT = f'{D}/fragmentation'


def fasta_lengths(path):
    lens, name, n = {}, None, 0
    for l in open(path):
        if l.startswith('>'):
            if name:
                lens[name] = n
            name, n = l[1:].split()[0], 0
        else:
            n += len(l.strip())
    if name:
        lens[name] = n
    return lens


def gff_cds(path):
    """contig -> sorted array of (start, end) CDS coordinates, 1-based inclusive; plus locus_tag map."""
    by, tags = defaultdict(list), {}
    for l in open(path):
        if l.startswith('#') or '\t' not in l:
            continue
        f = l.split('\t')
        if len(f) < 9 or f[2] != 'CDS':
            continue
        s, e = int(f[3]), int(f[4])
        by[f[0]].append((s, e))
        for kv in f[8].split(';'):
            if kv.startswith('locus_tag='):
                tags[kv[10:].strip()] = (f[0], s, e)
    return {c: np.array(sorted(v)) for c, v in by.items()}, tags


def parse_gene(name):
    """Spacedust gene id contig_idx_start_end (start > end on minus strand)."""
    c, i, a, b = name.rsplit('_', 3)
    a, b = int(a), int(b)
    return c, int(i), min(a, b), max(a, b)


def sd_instances(path):
    """Unique cluster instances (one side of one pairwise Spacedust cluster): list of
    (contig, members) with members a sorted tuple of (start, end). Deduplicated by gene set."""
    seen, out, cur = set(), [], None

    def flush(cur):
        if not cur:
            return
        for side in cur:
            by = defaultdict(list)
            for g in side:
                by[g[0]].append(g)
            for c, gs in by.items():
                key = (c, tuple(sorted((g[2], g[3]) for g in gs)))
                if key in seen or len(gs) < 2:
                    continue
                seen.add(key)
                out.append(key)

    for l in open(path):
        if l.startswith('#'):
            flush(cur)
            cur = (set(), set())
        elif l.startswith('>'):
            q, t = l[1:].split('\t')[:2]
            cur[0].add(parse_gene(q)); cur[1].add(parse_gene(t))
    flush(cur)
    return out


def eco_operons():
    """Curated E. coli multi-gene transcription units from adjacent same-TU pairs (chains)."""
    _, tags = gff_cds(f'{BENCH}/ecoli/mag_out/annotation/bakta/basic/annotation.gff3')
    nxt = {}
    rows = list(csv.DictReader(open(f'{D}/bench/eco_pairs.tsv'), delimiter='\t'))
    for r in rows:
        if r['same_tu'] == '1':
            nxt[r['locus_a']] = r['locus_b']
    has_prev = set(nxt.values())
    ops = []
    for a in nxt:
        if a in has_prev:
            continue
        chain = [a]
        while chain[-1] in nxt:
            chain.append(nxt[chain[-1]])
        co = [tags[t] for t in chain if t in tags]
        ops.append((co[0][0], tuple(sorted((x[1], x[2]) for x in co))))
    return ops  # (contig, members)


def checkm2():
    q = {}
    for r in csv.DictReader(open(f'{D}/chesterfield/mag_chesterfield/binning/checkm2/quality_report.tsv'), delimiter='\t'):
        q[r['Name']] = (float(r['Completeness']) / 100, float(r['Contamination']) / 100)
    return q


def bin_names():
    return sorted(os.path.basename(f)[:-4] for f in glob.glob(f'{BINS}/*.fna'))


def n50(lens):
    l = np.sort(np.asarray(lens))[::-1]
    return int(l[np.searchsorted(np.cumsum(l), l.sum() / 2)])
