"""Pfam domains (hmmscan --cut_ga) of rare-present module families, next to the Bakta products.
Writes pfam/module_functions.tsv and prints markdown: class counts per family (Bakta vs Pfam) and the
modules with any informative annotation. usage: pfam_summary.py SURVEY_DIR"""
import csv, os, sys
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from classes import fclass  # noqa: E402
W = sys.argv[1]
CONTAM = {'tmed14'}
UNINF = ('hypothetical/DUF', 'no hit', 'other')
dom = defaultdict(list)
for l in open(f'{W}/pfam/query.domtbl'):
    if l.startswith('#'): continue
    f = l.split(None, 22)
    q, name, desc = f[3], f[0], f[22].strip()
    if (name, desc) not in dom[q]: dom[q].append((name, desc))
fam_info, mod_bins = {}, defaultdict(set)
modk = {}
for tier in ('all', 'small'):
    for r in csv.DictReader(open(f'{W}/{tier}/flagged_modules.tsv'), delimiter='\t'):
        mod_bins[(r['species'], r['module'])].add(r['bin']); modk[(r['species'], r['module'])] = (r['k'], r['n'], r['anchored'], r['mge'])
    for r in csv.DictReader(open(f'{W}/{tier}/flagged_genes.tsv'), delimiter='\t'):
        if r['species'] in CONTAM: continue
        k = (r['species'], r['family'])
        main = not r['bin'].startswith(('chesterfield_marine', 'chesterfield_fresh'))
        if k not in fam_info or (main and fam_info[k]['product'] in ('?', 'hypothetical protein')):
            fam_info[k] = dict(species=r['species'], module=r['module'], family=r['family'], product=r['product'])
bak, pf = Counter(), Counter()
rows = []
for (sp, fam), r in fam_info.items():
    d = dom.get(f'{sp}|{fam}', [])
    r['pfam'] = '; '.join(f'{n} ({desc})' for n, desc in d)
    r['bakta_class'] = fclass(r['product'])
    pc = [fclass(desc + ' ' + n) for n, desc in d]
    inf = [c for c in pc if c not in UNINF]
    r['pfam_class'] = inf[0] if inf else ('hypothetical/DUF' if pc else 'no domain')
    r['class'] = r['bakta_class'] if r['bakta_class'] not in UNINF else r['pfam_class']
    bak[r['bakta_class']] += 1; pf[r['class']] += 1
    rows.append(r)
with open(f'{W}/pfam/module_functions.tsv', 'w') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)
n = len(rows)
print(f'{n} families (TMED14 excluded); with a Pfam domain: {sum(1 for r in rows if r["pfam"])}; '
      f'Bakta hypothetical/no hit: {sum(1 for r in rows if r["bakta_class"] in UNINF)}, '
      f'of which Pfam-informative: {sum(1 for r in rows if r["bakta_class"] in UNINF and r["pfam_class"] not in UNINF + ("no domain",))}')
print('\n| Class | Families (Bakta only) | Families (Bakta, else Pfam) |\n|---|---|---|')
for c in sorted(set(bak) | set(pf), key=lambda c: -pf[c]):
    print(f'| {c} | {bak[c]} | {pf[c]} |')
bymod = defaultdict(list)
for r in rows: bymod[(r['species'], r['module'])].append(r)
print('\n| Species | Module | Fam. | k/n | Anchored | Classes | Annotations (Bakta / Pfam) |\n|---|---|---|---|---|---|---|')
for (sp, m), rs in sorted(bymod.items(), key=lambda kv: (kv[0][0], -len(kv[1]))):
    cls = Counter(r['class'] for r in rs)
    if all(c in UNINF + ('no domain',) for c in cls): continue
    k, nn, anch, mge = modk[(sp, m)]
    ann = sorted({r['product'] for r in rs if r['product'] not in ('?', 'hypothetical protein')} |
                 {x.split(' (')[1].rstrip(')') for r in rs if r['bakta_class'] in UNINF for x in r['pfam'].split('; ') if x and 'DUF' not in x})
    print(f"| {sp} | {m} | {len(rs)} | {k}/{nn} | {anch} | " + ', '.join(f'{c} {v}' for c, v in cls.most_common() if c not in UNINF + ('no domain',))
          + f" | {'; '.join(ann)[:220]} |")
