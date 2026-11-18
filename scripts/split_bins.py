"""Write one FASTA + GFF3 (CDS only) per bin from a co-assembly: contig2bin.tsv, bin contigs FASTA, Bakta GFF3.
usage: split_bins.py CONTIG2BIN BINS_FASTA BAKTA_GFF3 OUTDIR  -> OUTDIR/<bin>.fna, OUTDIR/<bin>.gff3, OUTDIR/pairs.txt"""
import os, sys
c2b_path, fa, gff, out = sys.argv[1:5]
os.makedirs(out, exist_ok=True)
c2b = dict(l.rstrip('\n').split('\t')[:2] for l in open(c2b_path) if '\t' in l and not l.startswith(('contig\t', 'contig_id')))
fh = {}
def handle(b, ext):
    k = (b, ext)
    if k not in fh:
        fh[k] = open(f'{out}/{b}.{ext}', 'w')
        if ext == 'gff3': fh[k].write('##gff-version 3\n')
    return fh[k]
cur = None
for l in open(fa):
    if l.startswith('>'):
        cur = c2b.get(l[1:].split()[0])
    if cur: handle(cur, 'fna').write(l)
for l in open(gff):
    if l.startswith('##FASTA'): break
    if l.startswith('#'): continue
    f = l.split('\t')
    if len(f) > 2 and f[2] == 'CDS' and f[0] in c2b:
        handle(c2b[f[0]], 'gff3').write(l)
for h in fh.values(): h.close()
bins = sorted({b for b, _ in fh})
with open(f'{out}/pairs.txt', 'w') as o:
    for b in bins:
        if os.path.exists(f'{out}/{b}.gff3'): o.write(f'{b}:{out}/{b}.fna:{out}/{b}.gff3\n')
print(f'{len(bins)} bins written', file=sys.stderr)
