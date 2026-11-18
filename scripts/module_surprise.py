"""Prototype surprise for PPanGGOLiN modules in our bins, against the species' public genomes.
p = (n_with + 1) / (n + 2) over public genomes (Beta(1,1)); present: -log2 p; absent: -log2(1 - p*c),
c = the bin's completeness from persistent gene families (PPanGGOLiN genomes_statistics).
usage: module_surprise.py PPANGGOLIN_DIR OUR_GENOME_PREFIX [OUR_GENOME_PREFIX ...]"""
import csv, math, sys
from collections import defaultdict
d, ours = sys.argv[1], tuple(sys.argv[2:])
stats = {}
rows = [l for l in open(f'{d}/genomes_statistics.tsv') if not l.startswith('#')]
for r in csv.DictReader(rows, delimiter='\t'):
    stats[r['Genome_name']] = r
genomes = sorted(stats)
pub = [g for g in genomes if not g.startswith(ours)]; mine = [g for g in genomes if g.startswith(ours)]
comp = {g: float(stats[g].get('Completeness', 'nan') or 'nan') / (100 if float(stats[g].get('Completeness', 0) or 0) > 1.5 else 1) for g in genomes}
present = defaultdict(set)
for r in csv.DictReader(open(f'{d}/modules_in_genomes.tsv'), delimiter='\t'):
    if float(r['completion']) >= 0.5: present[r['module_id']].add(r['genome'])
size = defaultdict(int)
for r in csv.DictReader(open(f'{d}/functional_modules.tsv'), delimiter='\t'): size[r['module_id']] += 1
fams = sum(1 for _ in open(f'{d}/gene_families.tsv'))
part = defaultdict(int)
import glob, os
for f in glob.glob(f'{d}/partitions/*.txt'): part[os.path.basename(f)[:-4]] = sum(1 for _ in open(f))
print(f'{len(genomes)} genomes ({len(pub)} public, {len(mine)} ours); gene families {len(set(l.split()[0] for l in open(d + "/gene_families.tsv")))}; '
      f'partitions {dict(part)}; modules {len(size)} (median {sorted(size.values())[len(size) // 2] if size else 0} families)')
print('  our bins: ' + ', '.join(f'{g} completeness {comp[g]:.2f}' for g in mine))
out = []
for m in size:
    k = len(present[m] & set(pub)); p = (k + 1) / (len(pub) + 2)
    for g in mine:
        if g in present[m]: out.append((-math.log2(p), 'present', m, g, k, size[m]))
        else: out.append((-math.log2(1 - p * comp[g]), 'absent', m, g, k, size[m]))
out.sort(reverse=True)
print(f'  most surprising module states (prevalence among {len(pub)} public genomes):')
for s, st, m, g, k, n in out[:10]:
    print(f'    {s:5.2f} bits  {st:7s} {m:10s} ({n:2d} families, in {k:2d}/{len(pub)} public)  {g}')
tot = defaultdict(float)
for s, st, m, g, k, n in out: tot[g] += s
print('  total surprise per bin over all modules: ' + ', '.join(f'{g} {v:.1f} bits' for g, v in tot.items()))
