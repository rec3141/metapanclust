"""Adjacent same-strand gene pairs of a Bakta-annotated genome, with intergenic distance, labelled from
curated transcription units when given (EcoCyc genes.dat + transunits.dat: same TU = same operon).
usage: operon_pairs.py BAKTA_TSV OUT_TSV [ECOCYC_DIR]"""
import sys
from collections import defaultdict
tsv, out = sys.argv[1:3]
eco = sys.argv[3] if len(sys.argv) > 3 else None
genes = defaultdict(list)                       # contig -> [(start, end, strand, locus)]
for l in open(tsv):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if len(f) > 5 and f[1] == 'cds':
        genes[f[0]].append((int(f[2]), int(f[3]), f[4], f[5]))
lab = {}
if eco:
    def recs(path):
        r = defaultdict(list)
        for l in open(path, encoding='latin-1'):
            if l.startswith('//'):
                if r: yield r
                r = defaultdict(list)
            elif ' - ' in l and not l.startswith('#'):
                k, v = l.rstrip('\n').split(' - ', 1); r[k].append(v)
    pos = {}
    for r in recs(f'{eco}/genes.dat'):
        if r.get('LEFT-END-POSITION') and r.get('RIGHT-END-POSITION'):
            pos[r['UNIQUE-ID'][0]] = (int(r['LEFT-END-POSITION'][0]), int(r['RIGHT-END-POSITION'][0]))
    tus = defaultdict(set)                      # EcoCyc gene -> TU ids
    for r in recs(f'{eco}/transunits.dat'):
        for c in r.get('COMPONENTS', []):
            if c in pos: tus[c].add(r['UNIQUE-ID'][0])
    # Bakta CDS -> EcoCyc gene by best reciprocal overlap (same coordinates system: NC_000913.3 / U00096.3)
    starts = sorted((s, e, g) for g, (s, e) in pos.items())
    import bisect
    keys = [s for s, _, _ in starts]
    b2e = {}
    for c, gl in genes.items():
        for s, e, st, loc in gl:
            i = bisect.bisect_left(keys, s - 20000)
            best = (0, None)
            for s2, e2, g in starts[i:]:
                if s2 > e: break
                ov = min(e, e2) - max(s, s2) + 1
                if ov > 0:
                    fr = ov / max(e - s + 1, e2 - s2 + 1)
                    if fr > best[0]: best = (fr, g)
            if best[0] >= 0.8: b2e[loc] = best[1]
n = 0
with open(out, 'w') as o:
    o.write('contig\tlocus_a\tlocus_b\tstrand\tgap\tsame_tu\n')
    for c, gl in genes.items():
        gl.sort()
        for (s1, e1, st1, a), (s2, e2, st2, b) in zip(gl, gl[1:]):
            if st1 != st2: continue
            y = ''
            if eco and a in b2e and b in b2e:
                ta, tb = tus.get(b2e[a], set()), tus.get(b2e[b], set())
                if ta and tb: y = int(bool(ta & tb))
            o.write(f'{c}\t{a}\t{b}\t{st1}\t{s2 - e1 - 1}\t{y}\n'); n += 1
print(f'{n} adjacent same-strand pairs; mapped to EcoCyc: {len(lab) or (len(b2e) if eco else 0)} genes', file=sys.stderr)
