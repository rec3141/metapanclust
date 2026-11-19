"""E. coli adjacent pairs called conserved by run B but not run A (same co-membership rule as eval_pairs.py).
usage: gained_pairs.py PAIRS_TSV A.tsv B.tsv BAKTA_TSV"""
import csv, sys
from collections import defaultdict
pairs, ta, tb, bakta = sys.argv[1:5]
loc = {}
for l in open(bakta):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if len(f) > 5 and f[1] == 'cds': loc[(f[0], int(f[2]), int(f[3]))] = (f[5], f[7] if len(f) > 7 else '')
def members(tsv):
    m = defaultdict(set); cid = None
    for l in open(tsv):
        if l.startswith('#'): cid = l[1:].split('\t')[0]; continue
        for p in l[1:].split('\t')[:2]:
            acc, i, a, b = p.rsplit('_', 3); g = loc.get((acc, min(int(a), int(b)), max(int(a), int(b))))
            if g: m[g[0]].add(cid)
    return m
ma, mb = members(ta), members(tb)
prod = {v[0]: v[1] for v in loc.values()}
for r in csv.DictReader(open(pairs), delimiter='\t'):
    if r['same_tu'] == '': continue
    ca = bool(ma[r['locus_a']] & ma[r['locus_b']]); cb = bool(mb[r['locus_a']] & mb[r['locus_b']])
    if ca != cb:
        print(f"{'gained' if cb else 'lost'}\tsame_tu={r['same_tu']}\tgap={r['gap']}\t{prod[r['locus_a']][:35]} / {prod[r['locus_b']][:35]}")
