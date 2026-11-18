"""Genes of modules found in our bins but in no public genome: Bakta product and whether the contig is a
geNomad plasmid/virus/provirus contig. usage: module_genes.py PPANGGOLIN_DIR OUR_GENOME BAKTA_TSV GENOMAD_DIR"""
import csv, glob, os, sys
from collections import defaultdict
d, ours, bakta, gm = sys.argv[1:5]
pub_with = defaultdict(set)
for r in csv.DictReader(open(f'{d}/modules_in_genomes.tsv'), delimiter='\t'):
    if float(r['completion']) >= 0.5: pub_with[r['module_id']].add(r['genome'])
mods = [m for m, gs in pub_with.items() if ours in gs and not any(not g.startswith('chesterfield') for g in gs)]
fam2mod = {}
for r in csv.DictReader(open(f'{d}/functional_modules.tsv'), delimiter='\t'):
    if r['module_id'] in mods: fam2mod[r['family_id']] = r['module_id']
genes = {}                                      # our gene id -> module
for l in open(f'{d}/gene_families.tsv'):
    f = l.split()
    if f[0] in fam2mod and f[1].startswith(ours): genes[f[1]] = fam2mod[f[0]]
coords = {}
for gf in glob.glob(f'{d}/gff/{ours}.gff'):
    for l in open(gf):
        if l.startswith('#') or '\tCDS\t' not in l: continue
        f = l.split('\t'); gid = dict(x.split('=', 1) for x in f[8].strip().split(';') if '=' in x).get('ID', '')
        if gid in genes: coords[gid] = (f[0], int(f[3]), int(f[4]))
prod = {}
for l in open(bakta):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if len(f) > 7 and f[1] == 'cds': prod[(f[0], int(f[2]), int(f[3]))] = f[7]
mge = {}
for kind, fn in (('plasmid', 'plasmid_summary.tsv'), ('virus', 'virus_summary.tsv')):
    p = f'{gm}/{fn}'
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter='\t'): mge[r['seq_name'].split('|')[0]] = kind
by = defaultdict(list)
for g, m in genes.items():
    c = coords.get(g)
    pr = prod.get(c, '?') if c else '?'
    by[m].append((c[0] if c else '?', pr, mge.get(c[0], '') if c else ''))
print(f'{len(mods)} modules in {ours} and no public genome; {len(genes)} genes')
for m, gl in sorted(by.items(), key=lambda kv: -len(kv[1]))[:12]:
    contigs = {c for c, _, _ in gl}; kinds = {k for _, _, k in gl if k}
    print(f'  {m}: {len(gl)} genes on {len(contigs)} contig(s) {sorted(kinds) or ""}: ' + '; '.join(sorted({p for _, p, _ in gl})[:6]))
