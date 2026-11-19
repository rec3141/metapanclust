"""Representative proteins of rare-present module families from both survey tiers, TMED14 left out.
usage: pfam_query.py SURVEY_DIR > query.faa"""
import sys
W = sys.argv[1]
seqs = {}
for tier in ('all', 'small'):
    name = None
    for l in open(f'{W}/{tier}/flagged_families.faa'):
        if l.startswith('>'): name = l[1:].strip()
        else: seqs[name] = seqs.get(name, '') + l.strip()
for k, v in seqs.items():
    if not k.startswith('tmed14|'): print(f'>{k}\n{v}')
