"""Score adjacent gene pairs against curated operons: distance rule, Spacedust co-membership, and both.
usage: eval_pairs.py PAIRS_TSV SPACEDUST_TSV CONTIG_PREFIX [BAKTA_TSV]
A pair is 'conserved' when both genes are members of one Spacedust cluster (query or target side)."""
import csv, sys
from collections import defaultdict
pairs_tsv, sd_tsv = sys.argv[1:3]
bakta = sys.argv[3]
loc = {}                                         # (contig, lo, hi) -> locus
for l in open(bakta):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if len(f) > 5 and f[1] == 'cds': loc[(f[0], int(f[2]), int(f[3]))] = f[5]
def to_locus(pid):
    parts = pid.rsplit('_', 3)                  # accession, index, start, end
    a, b = int(parts[2]), int(parts[3])
    return loc.get((parts[0], min(a, b), max(a, b)))
member = defaultdict(set)                        # locus -> cluster ids (genes of this genome)
cid = None; miss = 0; n = 0
for l in open(sd_tsv):
    if l.startswith('#'): cid = l[1:].split('\t')[0]; continue
    if l.startswith('>'):
        q, t = l[1:].split('\t')[:2]
        for p in (q, t):
            n += 1; g = to_locus(p)
            if g: member[g].add(cid)
            elif p.startswith(tuple(c for c, _, _ in list(loc)[:1])): miss += 1
rows = [r for r in csv.DictReader(open(pairs_tsv), delimiter='\t') if r['same_tu'] != '']
def report(lab, pred):
    tp = sum(1 for r in rows if pred(r) and r['same_tu'] == '1'); fp = sum(1 for r in rows if pred(r) and r['same_tu'] == '0')
    fn = sum(1 for r in rows if not pred(r) and r['same_tu'] == '1'); tn = len(rows) - tp - fp - fn
    print(f'  {lab:40s} called {tp + fp:5d}  accuracy {(tp + tn) / len(rows):.3f}  precision {tp / max(tp + fp, 1):.3f}  recall {tp / max(tp + fn, 1):.3f}')
cons = lambda r: bool(member[r['locus_a']] & member[r['locus_b']])
print(f'{len(rows)} labelled pairs ({sum(r["same_tu"] == "1" for r in rows)} same TU); genes in any cluster: {sum(1 for g in member if member[g])}; unmapped ids: {miss}')
report('gap <= 50 bp', lambda r: int(r['gap']) <= 50)
report('conserved (Spacedust)', cons)
report('conserved AND gap <= 150', lambda r: cons(r) and int(r['gap']) <= 150)
report('gap <= 50 OR (conserved AND gap <= 300)', lambda r: int(r['gap']) <= 50 or (cons(r) and int(r['gap']) <= 300))
# among pairs the distance rule gets wrong
far_same = [r for r in rows if int(r['gap']) > 50 and r['same_tu'] == '1']
near_diff = [r for r in rows if int(r['gap']) <= 50 and r['same_tu'] == '0']
print(f'  distance-rule misses: {len(far_same)} far same-TU pairs, {sum(map(cons, far_same))} conserved; '
      f'{len(near_diff)} close different-TU pairs, {sum(map(cons, near_diff))} conserved')
