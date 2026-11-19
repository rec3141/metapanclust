"""Markdown tables from survey.py outputs. usage: summarize.py OUTDIR"""
import csv, sys
from collections import Counter, defaultdict
o = sys.argv[1]
rd = lambda f: list(csv.DictReader(open(f'{o}/{f}'), delimiter='\t'))
sp, pb, fm, fg = rd('species.tsv'), rd('per_bin.tsv'), rd('flagged_modules.tsv'), rd('flagged_genes.tsv')

print('## Species\n\n| Species | Public (median compl.) | Ours | Families | Pers./shell/cloud | Modules |\n|---|---|---|---|---|---|')
for r in sp:
    print(f"| {r['species']} | {r['public']} ({float(r['public_completeness_median']):.2f}) | {r['ours']} | {r['families']} | "
          f"{r['persistent']}/{r['shell']}/{r['cloud']} | {r['modules']} |")

print('\n## Per bin\n\n| Species | Bin | Compl. | Modules present / all | Rare-present (k/n <= 0.1) | In no public genome | Bits total | present | absent |\n|---|---|---|---|---|---|---|---|---|')
for r in sorted(pb, key=lambda r: -float(r['bits_total'])):
    print(f"| {r['species']} | {r['bin'].replace('chesterfield_', '')} | {float(r['completeness']):.2f} | {r['present']}/{r['modules']} | "
          f"{r['rare_present']} | {r['absent_from_all_public']} | {r['bits_total']} | {r['bits_present']} | {r['bits_absent']} |")

# unique (species, module): pick one bin's view (prefer main-assembly bin: exact Bakta products and geNomad)
uniq = {}
for r in fm:
    key = (r['species'], r['module'])
    main = not r['bin'].startswith(('chesterfield_marine', 'chesterfield_freshwater'))
    if key not in uniq or main: uniq[key] = r
print(f'\n## Rare-present modules: {len(uniq)} unique (species, module), {len(fm)} (bin, module) pairs\n')
cls = Counter(r['class'] for r in uniq.values())
fams = defaultdict(int)
for r in uniq.values(): fams[r['class']] += int(r['families'])
print('| Leading informative class | Modules | Families |\n|---|---|---|')
for c, n in cls.most_common():
    print(f'| {c} | {n} | {fams[c]} |')
# gene-level classes (unique family per species)
gcls = Counter()
seen = set()
for r in fg:
    k = (r['species'], r['family'])
    if k in seen: continue
    seen.add(k); gcls[r['cls']] += 1
print('\nGene families of rare-present modules by product class: ' + ', '.join(f'{c} {n}' for c, n in gcls.most_common()))
anch = sum(int(r['anchored']) for r in uniq.values())
mge = Counter(r['mge'] or 'none' for r in uniq.values())
print(f'\nAnchored (a persistent gene on a module contig): {anch}/{len(uniq)}; geNomad: {dict(mge)}')
both = sum(1 for r in uniq.values() if int(r['carried_by_ours']) >= 2)
print(f'Carried by >= 2 of our bins: {both}/{len(uniq)}')
print('\n## Rare-present modules per species\n\n| Species | Modules | k=0 | anchored | MGE | classes |\n|---|---|---|---|---|---|')
bysp = defaultdict(list)
for r in uniq.values(): bysp[r['species']].append(r)
for s, rs in sorted(bysp.items(), key=lambda kv: -len(kv[1])):
    c = Counter(r['class'] for r in rs)
    print(f"| {s} | {len(rs)} | {sum(r['k'] == '0' for r in rs)} | {sum(int(r['anchored']) for r in rs)} | "
          f"{sum(1 for r in rs if r['mge'] and r['mge'] != 'n/a')} | " + ', '.join(f'{k} {v}' for k, v in c.most_common()) + ' |')
print('\n## Informative rare-present modules (largest first)\n\n| Species | Module | Fam. | k/n | Bins | Anchored | MGE | Class | Products |\n|---|---|---|---|---|---|---|---|---|')
for r in sorted(uniq.values(), key=lambda r: (-int(r['families']))):
    if r['class'] in ('hypothetical/DUF', 'no hit', 'none', 'other'): continue
    print(f"| {r['species']} | {r['module']} | {r['families']} | {r['k']}/{r['n']} | {r['carried_by_ours']} | {r['anchored']} | "
          f"{r['mge'] or '-'} | {r['class']} | {r['products'][:160]} |")
