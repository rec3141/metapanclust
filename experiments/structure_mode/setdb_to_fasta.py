"""Write an MMseqs2/Spacedust amino-acid DB as FASTA with the DB key as header, in key order.
usage: setdb_to_fasta.py SETDB OUT.faa [MAX_SEQS]"""
import sys
db, out = sys.argv[1:3]
limit = int(sys.argv[3]) if len(sys.argv) > 3 else None
data = open(db, 'rb').read()
idx = sorted((int(k), int(o), int(l)) for k, o, l in (x.split('\t') for x in open(db + '.index')))
with open(out, 'w') as fh:
    for i, (k, o, l) in enumerate(idx):
        if limit and i >= limit: break
        fh.write(f'>{k}\n{data[o:o + l - 1].decode().strip()}\n')   # entries end with '\n\0'
