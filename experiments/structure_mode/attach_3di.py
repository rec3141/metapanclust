"""Attach a ProstT5 3Di DB (foldseek createdb on setdb_to_fasta.py output) to a Spacedust setDB as SETDB_ss.
Keys are remapped through the FASTA headers (= setDB keys); 3Di lengths are checked against the setDB.
usage: attach_3di.py SETDB FS_DB"""
import sys
setdb, fs = sys.argv[1:3]
name = {int(k): n for k, n, *_ in (l.rstrip('\n').split('\t') for l in open(fs + '.lookup'))}
fidx = {int(k): (int(o), int(l)) for k, o, l in (x.split('\t') for x in open(fs + '_ss.index'))}
sidx = {int(k): int(l) for k, o, l in (x.split('\t') for x in open(setdb + '.index'))}
data = open(fs + '_ss', 'rb').read()
out = open(setdb + '_ss', 'wb'); idx = []; off = 0; bad = 0
for fk in sorted(fidx, key=lambda k: int(name[k])):
    sk = int(name[fk]); o, l = fidx[fk]; rec = data[o:o + l]
    if l != sidx[sk]: bad += 1
    out.write(rec); idx.append(f'{sk}\t{off}\t{l}\n'); off += l
out.close(); open(setdb + '_ss.index', 'w').writelines(idx)
open(setdb + '_ss.dbtype', 'wb').write(open(fs + '_ss.dbtype', 'rb').read())
print(f'{len(idx)} 3Di records attached ({len(sidx)} in setDB); length mismatches: {bad}')
