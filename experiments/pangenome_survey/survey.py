"""Module surprise across many species pangenomes (PPanGGOLiN), with function and MGE context for the
modules our Chesterfield bins carry but public genomes rarely do.

Surprise as in scripts/module_surprise.py: p = (k + 1) / (n + 2) over the n public genomes (k carry the
module, completion >= 0.5); present -log2 p, absent -log2(1 - p c), c = bin completeness (persistent
families). A module is "rare-present" in a bin if the bin carries it and k / n <= RARE.

Function: representative protein of each family of a rare-present module, searched (MMseqs2) against the
Bakta proteins of the main mag_chesterfield co-assembly. For genes of main-assembly bins the Bakta
product is taken by exact coordinates instead. MGE: contig in the geNomad plasmid/virus/provirus lists
(main assembly only; marine/freshwater bin genes use their >= 95% identity Bakta hit's contig as a proxy).
Anchoring: number of persistent-partition genes on the same contig (0 = contig carries no core genes,
so binning contamination cannot be excluded).

usage: survey.py OUTDIR SPECIES_TSV   (SPECIES_TSV: short<TAB>gtdb species<TAB>ppanggolin dir[<TAB>family protein faa])"""
import csv, gzip, math, os, re, subprocess, sys
from collections import defaultdict, Counter

RARE = 0.10
out, sp_tsv = sys.argv[1:3]
C = '/data/scratch/metapanclust_data/chesterfield/mag_chesterfield'
BAKTA = f'{C}/annotation/bakta/basic'
GENOMAD = f'{C}/mge/genomad'
OURS = 'chesterfield'

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from classes import fclass  # noqa: E402


def read_stats(d):
    rows = [l for l in open(f'{d}/genomes_statistics.tsv') if not l.startswith('#')]
    stats = {r['Genome_name']: r for r in csv.DictReader(rows, delimiter='\t')}
    comp = {}
    for g, r in stats.items():
        c = float(r.get('Completeness') or 'nan')
        comp[g] = c / 100 if c > 1.5 else c
    return stats, comp


species = [l.rstrip('\n').split('\t') for l in open(sp_tsv) if l.strip() and not l.startswith('#')]
os.makedirs(out, exist_ok=True)

# ---------- pass 1: surprise and rare-present modules per species ----------
per_bin, flagged, sp_rows, fam_reps = [], [], [], {}
sp_data = {}
famfa = {}
for row in species:
    short, spname, d = row[:3]
    famfa[short] = row[3] if len(row) > 3 else f'{d}/fasta/all_protein_families.faa'
    stats, comp = read_stats(d)
    genomes = sorted(stats)
    pub = [g for g in genomes if not g.startswith(OURS)]
    mine = [g for g in genomes if g.startswith(OURS)]
    present = defaultdict(set)
    for r in csv.DictReader(open(f'{d}/modules_in_genomes.tsv'), delimiter='\t'):
        if float(r['completion']) >= 0.5:
            present[r['module_id']].add(r['genome'])
    mfams = defaultdict(list)
    for r in csv.DictReader(open(f'{d}/functional_modules.tsv'), delimiter='\t'):
        mfams[r['module_id']].append(r['family_id'])
    nfam = len({l.split()[0] for l in open(f'{d}/gene_families.tsv')})
    part = {os.path.basename(f)[:-4]: sum(1 for _ in open(f'{d}/partitions/{f}'))
            for f in os.listdir(f'{d}/partitions') if f in ('persistent.txt', 'shell.txt', 'cloud.txt')}
    pc = sorted(comp[g] for g in pub)
    sp_rows.append(dict(species=spname, short=short, public=len(pub), ours=len(mine),
                        public_completeness_median=round(pc[len(pc) // 2], 3), families=nfam,
                        persistent=part.get('persistent', 0), shell=part.get('shell', 0),
                        cloud=part.get('cloud', 0), modules=len(mfams)))
    for g in mine:
        tot = pres_bits = abs_bits = 0.0
        n_pres = n_rare = n_zero = 0
        for m, fams in mfams.items():
            k = len(present[m] & set(pub)); p = (k + 1) / (len(pub) + 2)
            if g in present[m]:
                s = -math.log2(p); pres_bits += s; n_pres += 1
                if k / len(pub) <= RARE:
                    n_rare += 1; n_zero += (k == 0)
                    flagged.append(dict(species=short, module=m, bin=g, k=k, n=len(pub), families=len(fams),
                                        bits=round(s, 2), carried_by_ours=sum(1 for x in mine if x in present[m])))
            else:
                s = -math.log2(1 - p * comp[g]); abs_bits += s
            tot += s
        per_bin.append(dict(species=short, bin=g, completeness=round(comp[g], 3), public=len(pub),
                            modules=len(mfams), present=n_pres, rare_present=n_rare, absent_from_all_public=n_zero,
                            bits_total=round(tot, 1), bits_present=round(pres_bits, 1), bits_absent=round(abs_bits, 1)))
    sp_data[short] = (d, mfams, mine, pub, present)

# ---------- function: representative proteins of flagged module families vs Bakta ----------
flag_mods = {(f['species'], f['module']) for f in flagged}
qfa = f'{out}/flagged_families.faa'
buf = []
for short, m in sorted(flag_mods):
    d, mfams = sp_data[short][0], sp_data[short][1]
    want = set(mfams[m]); keep = False
    for l in open(famfa[short]):
        if l.startswith('>'):
            name = l[1:].split()[0]; keep = name in want
            if keep: buf.append(f'>{short}|{name}\n')
        elif keep:
            buf.append(l.rstrip('*\n') + '\n')
if not os.path.exists(qfa) or open(qfa).read() != ''.join(buf):   # rewrite (and re-search) only on change
    open(qfa, 'w').write(''.join(buf))
hits = f'{out}/flagged_vs_bakta.m8'
if not os.path.exists(hits) or os.path.getmtime(hits) < os.path.getmtime(qfa):
    subprocess.run(['mmseqs', 'easy-search', qfa, f'{BAKTA}/annotation.faa', hits, f'{out}/mmseqs_tmp',
                    '--threads', '2', '-e', '1e-5', '--max-seqs', '50', '--split-memory-limit', '8G',
                    '--format-output', 'query,target,fident,qcov,tcov,evalue,bits'],
                   check=True, stdout=subprocess.DEVNULL)
best = {}
for l in open(hits):
    q, t, fi, qc, tc, ev, b = l.split('\t')
    if q not in best or float(b) > best[q][3]:
        best[q] = (t, float(fi), float(qc), float(b))

lt2prod, lt2contig, coord2prod = {}, {}, {}
for l in open(f'{BAKTA}/annotation.tsv'):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if len(f) > 7 and f[1] == 'cds':
        lt2prod[f[5]] = f[7]; lt2contig[f[5]] = f[0]; coord2prod[(f[0], int(f[2]), int(f[3]))] = f[7]
mge = {}
for kind, fn in (('plasmid', 'plasmid_summary.tsv'), ('virus', 'virus_summary.tsv')):
    for r in csv.DictReader(open(f'{GENOMAD}/{fn}'), delimiter='\t'):
        c = r['seq_name'].split('|')[0]
        mge[c] = 'provirus' if '|provirus' in r['seq_name'] else kind


def fam_product(short, fam):
    h = best.get(f'{short}|{fam}')
    if not h: return '?', None
    t, fi, qc, _ = h
    return lt2prod.get(t, '?') + ('' if fi >= 0.5 else f' (~{fi:.0%})'), (lt2contig.get(t) if fi >= 0.95 and qc >= 0.8 else None)


# ---------- gene-level context in our bins ----------
gene_rows = []
mod_ctx = {}
for short, (d, mfams, mine, pub, present) in sp_data.items():
    fam_of = {}
    want_mods = {f['module'] for f in flagged if f['species'] == short}
    want_fams = {fa: m for m in want_mods for fa in mfams[m]}
    for l in open(f'{d}/gene_families.tsv'):
        f = l.split()
        if f[0] in want_fams and f[1].startswith(OURS): fam_of[f[1]] = f[0]
    for g in mine:
        gff = f'{d}/gff/{g}.gff'
        genes, contig_persist = [], Counter()
        for l in open(gff):
            if l.startswith('#') or '\tCDS\t' not in l: continue
            f = l.rstrip('\n').split('\t'); a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
            if a.get('partition') == 'persistent': contig_persist[f[0]] += 1
            if a.get('ID') in fam_of: genes.append((a['ID'], f[0], int(f[3]), int(f[4])))
        contig_len = {}
        name = None
        for l in open([x.split('\t')[1].strip() for x in open(f'{os.path.dirname(d)}/genomes.tsv') if x.split('\t')[0] == g][0]):
            if l.startswith('>'): name = l[1:].split()[0]; contig_len[name] = 0
            else: contig_len[name] += len(l.strip())
        main = not g.startswith(('chesterfield_marine', 'chesterfield_freshwater'))
        for gid, contig, s, e in genes:
            fam = fam_of[gid]; m = want_fams[fam]
            prod, proxy_contig = fam_product(short, fam)
            if main:
                prod = coord2prod.get((contig, s, e), prod)
                mk = mge.get(contig, '')
            else:
                mk = mge.get(proxy_contig, '') if proxy_contig else 'n/a'
            gene_rows.append(dict(species=short, bin=g, module=m, family=fam, gene=gid, contig=contig,
                                  contig_len=contig_len.get(contig, 0), persistent_on_contig=contig_persist[contig],
                                  product=prod, cls=fclass(prod), mge=mk))

# ---------- module summaries ----------
by_mod = defaultdict(list)
for r in gene_rows: by_mod[(r['species'], r['module'], r['bin'])].append(r)
for f in flagged:
    gl = by_mod.get((f['species'], f['module'], f['bin']), [])
    cls = Counter(r['cls'] for r in gl)
    informative = [c for c, _ in cls.most_common() if c not in ('hypothetical/DUF', 'no hit', 'other')]
    f['class'] = informative[0] if informative else (cls.most_common(1)[0][0] if cls else 'none')
    f['classes'] = '; '.join(f'{c}:{n}' for c, n in cls.most_common())
    f['genes'] = len(gl)
    f['contigs'] = len({r['contig'] for r in gl})
    f['anchored'] = int(any(r['persistent_on_contig'] > 0 for r in gl))
    f['max_contig_len'] = max((r['contig_len'] for r in gl), default=0)
    f['mge'] = ','.join(sorted({r['mge'] for r in gl if r['mge']})) or ''
    f['products'] = '; '.join(sorted({r['product'] for r in gl})[:10])


def write(fn, rows):
    if not rows: return
    with open(f'{out}/{fn}', 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)


write('species.tsv', sp_rows); write('per_bin.tsv', per_bin); write('flagged_modules.tsv', flagged)
write('flagged_genes.tsv', gene_rows)
print(f'{len(sp_rows)} species, {len(per_bin)} bins, {len(flagged)} rare-present (bin, module) pairs')
