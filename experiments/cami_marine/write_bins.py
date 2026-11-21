"""Write bin FASTAs from the pooled gold-standard assembly and a bin -> contig list (members.tsv.gz), so the 25 GB of
bin FASTAs need not be copied. Loads the assembly once (about 2.5 GB), then writes each bin in one go.
usage: write_bins.py ASSEMBLY.fasta.gz MEMBERS.tsv.gz OUTDIR
"""
import gzip, os, sys
from collections import defaultdict

asm, mem, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)
bins = defaultdict(list)
with gzip.open(mem, 'rt') as f:
    for line in f:
        b, c = line.rstrip('\n').split('\t'); bins[b].append(c)
need = {c for v in bins.values() for c in v}
seq = {}
with gzip.open(asm, 'rt') as fi:
    cur, buf = None, []
    for line in fi:
        if line.startswith('>'):
            if cur in need: seq[cur] = ''.join(buf)
            cur, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if cur in need: seq[cur] = ''.join(buf)
for b, cs in bins.items():
    if os.path.exists(f'{out}/{b}.fna'): continue
    with open(f'{out}/{b}.fna.tmp', 'w') as fo:
        for c in cs: fo.write(f'>{c}\n{seq[c]}\n')
    os.replace(f'{out}/{b}.fna.tmp', f'{out}/{b}.fna')
print(f'{len(bins)} bins written')
