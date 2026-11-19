"""Adjacent same-strand gene pairs of the Bakta-annotated B. subtilis 168 genome, labelled from SubtiWiki operons
(same operon = same transcription unit; adapted from scripts/operon_pairs.py, which reads EcoCyc).
usage: bsub_pairs.py BAKTA_TSV SUBTIWIKI_DIR OUT_TSV
SUBTIWIKI_DIR holds operon.json (v5 API /operon/) and gene_coords.tsv (fetch_subtiwiki.sh).
Label 1: the two genes share a SubtiWiki operon; 0: both genes are in some operon but never in the same one;
empty: at least one gene unmapped or in no operon."""
import bisect, json, sys
from collections import defaultdict
tsv, swd, out = sys.argv[1:4]
genes = defaultdict(list)
for l in open(tsv):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if len(f) > 5 and f[1] == 'cds':
        genes[f[0]].append((int(f[2]), int(f[3]), f[4], f[5]))
pos = {}
for l in list(open(f'{swd}/gene_coords.tsv'))[1:]:
    gid, name, lt, s, e, st = l.rstrip('\n').split('\t')
    if s and e: pos[gid] = (min(int(s), int(e)), max(int(s), int(e)))
tus = defaultdict(set)
ops = json.load(open(f'{swd}/operon.json'))['data']
for o in ops:
    for g in o['genes']:
        if str(g['id']) in pos: tus[str(g['id'])].add(o['id'])
# Bakta CDS -> SubtiWiki gene by best reciprocal overlap (both on NC_000964.3)
starts = sorted((s, e, g) for g, (s, e) in pos.items())
keys = [s for s, _, _ in starts]
b2s = {}
for c, gl in genes.items():
    for s, e, st, loc in gl:
        best = (0, None)
        for s2, e2, g in starts[bisect.bisect_left(keys, s - 20000):]:
            if s2 > e: break
            ov = min(e, e2) - max(s, s2) + 1
            if ov > 0:
                fr = ov / max(e - s + 1, e2 - s2 + 1)
                if fr > best[0]: best = (fr, g)
        if best[0] >= 0.8: b2s[loc] = best[1]
n = lab = pos1 = 0
with open(out, 'w') as o:
    o.write('contig\tlocus_a\tlocus_b\tstrand\tgap\tsame_tu\n')
    for c, gl in genes.items():
        gl.sort()
        for (s1, e1, st1, a), (s2, e2, st2, b) in zip(gl, gl[1:]):
            if st1 != st2: continue
            y = ''
            if a in b2s and b in b2s:
                ta, tb = tus.get(b2s[a], set()), tus.get(b2s[b], set())
                if ta and tb: y = int(bool(ta & tb)); lab += 1; pos1 += y
            o.write(f'{c}\t{a}\t{b}\t{st1}\t{s2 - e1 - 1}\t{y}\n'); n += 1
print(f'{n} adjacent same-strand pairs; {len(b2s)} of {sum(map(len, genes.values()))} Bakta CDS mapped to SubtiWiki genes; '
      f'{len(ops)} operons; labelled pairs {lab}, same operon {pos1}', file=sys.stderr)
