"""List hits whose query gene appears only in run B (not in run A), with Bakta products of both genes.
usage: new_hits.py A.tsv B.tsv GFF3 [GFF3 ...]"""
import sys
a, b, *gffs = sys.argv[1:]
prod = {}
for g in gffs:
    for l in open(g):
        f = l.rstrip('\n').split('\t')
        if len(f) > 8 and f[2] == 'CDS':
            p = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x).get('product', '?')
            prod[(f[0], int(f[3]), int(f[4]))] = p
def name(pid):
    acc, i, s, e = pid.rsplit('_', 3); s, e = int(s), int(e)
    return prod.get((acc, min(s, e), max(s, e)), '?')
seen = set()
for l in open(a):
    if l.startswith('>'): seen.update(l[1:].split('\t')[:2])
cid = None
for l in open(b):
    if l.startswith('#'): cid = l.split('\t')[0]; continue
    f = l[1:].rstrip('\n').split('\t')
    if f[0] not in seen or f[1] not in seen:
        print(f'{cid}\tid {f[3]}\tE {f[4]}\t{name(f[0])[:45]:45s} | {name(f[1])[:45]}')
