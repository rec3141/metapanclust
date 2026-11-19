"""Module and family surprise of Chesterfield bins projected onto PanGBank pangenomes, side by side with
the pangenome survey's own builds (same model as experiments/pangenome_survey/survey.py).

Surprise: p = (k + 1) / (n + 2) over the n PanGBank genomes (k carry the module at completion >= 0.5);
present -log2 p, absent -log2(1 - p c), c = bin completeness reported by the projection (persistent
families). Rare-present: carried and k / n <= 0.1.

Projection can only report modules that exist in PanGBank, so bin content absent from every PanGBank
genome shows up as "specific" genes instead. Specific genes are grouped into islands (>= 3 specific
genes in a row on one contig, in gene order) as the analogue of the survey's k = 0 modules.

usage: compare.py PANGBANK_DIR SURVEY_DIR OUTDIR
"""
import csv, math, os, sys
from collections import defaultdict

PB, SV, OUT = sys.argv[1:4]
RARE = 0.10
TARGETS = [  # tag, survey short name, survey ppanggolin dir
    ('mag120531_r226', 'mag120531', '/data/scratch/metapanclust_data/pangenomes/mag120531/ppanggolin'),
    ('uba9145_r226', 'uba9145', '/data/scratch/metapanclust_data/pangenomes/uba9145/ppanggolin'),
    ('uba9145_r232', 'uba9145', '/data/scratch/metapanclust_data/pangenomes/uba9145/ppanggolin'),
    ('pseudothioglobus_r226', 'pseudothioglobus', f'{SV}/pseudothioglobus/ppanggolin'),
    ('pseudothioglobus_r232', 'pseudothioglobus', f'{SV}/pseudothioglobus/ppanggolin'),
    ('patiriisocius_r226', 'patiriisocius', f'{SV}/patiriisocius/ppanggolin'),
]
OURS = 'chesterfield'
os.makedirs(OUT, exist_ok=True)


def tsv(path, skip_hash=True):
    rows = [l for l in open(path) if not (skip_hash and l.startswith('#'))]
    return list(csv.DictReader(rows, delimiter='\t'))


def rtab(path):
    with open(path) as f:
        hdr = f.readline().rstrip('\n').split('\t')[1:]
        fam = {}
        for l in f:
            x = l.rstrip('\n').split('\t')
            fam[x[0]] = {g for g, v in zip(hdr, x[1:]) if v != '0'}
    return hdr, fam


def gff_coords(path):
    out = {}
    for l in open(path):
        if l.startswith('#') or '\tCDS\t' not in l: continue
        f = l.split('\t'); gid = f[8].split('ID=')[1].split(';')[0]
        out[gid] = (f[0], int(f[3]), int(f[4]))
    return out


sv_bins = {(r['species'], r['bin']): r for r in tsv(f'{SV}/all/per_bin.tsv')}
sv_flag = defaultdict(list)
for r in tsv(f'{SV}/all/flagged_modules.tsv'):
    sv_flag[(r['species'], r['bin'])].append(r)
sv_genes = defaultdict(list)
for r in tsv(f'{SV}/all/flagged_genes.tsv'):
    sv_genes[(r['species'], r['bin'], r['module'])].append(r)

per_bin, fate_rows, island_rows = [], [], []
for tag, short, svd in TARGETS:
    T = f'{PB}/tables/{tag}'
    genomes, fam_genomes = rtab(f'{T}/rtab/gene_presence_absence.Rtab')
    n = len(genomes)
    present = defaultdict(set)
    for r in tsv(f'{T}/modules_in_genomes.tsv'):
        if float(r['completion']) >= 0.5: present[r['module_id']].add(r['genome'])
    mods = {r['module_id'] for r in tsv(f'{T}/functional_modules.tsv')}
    # survey build: public genomes and family prevalence among them
    sv_hdr, sv_fam = rtab(f'{svd}/gene_presence_absence.Rtab')
    sv_pub = [g for g in sv_hdr if not g.startswith(OURS)]
    overlap = len(set(sv_pub) & set(genomes))
    proj_root = f'{PB}/proj/{tag}'
    for b in sorted(os.listdir(proj_root)):
        bd = f'{proj_root}/{b}'
        if not os.path.isdir(bd): continue
        summ = {}
        for l in open(f'{bd}/projection_summary.yaml'):
            if ':' in l and not l.startswith(' '):
                k, v = l.split(':')[-2:] if l.startswith('Projection_summary') else l.split(':', 1)
                summ[k.strip()] = v.strip()
        c = float(summ['Completeness']) / 100
        mine = {r['module_id'] for r in tsv(f'{bd}/modules_in_input_genome.tsv') if float(r['completion']) >= 0.5}
        tot = pb = ab = 0.0; npres = nrare = 0
        for m in mods:
            k = len(present[m]); p = (k + 1) / (n + 2)
            if m in mine:
                s = -math.log2(p); pb += s; npres += 1; nrare += (k / n <= RARE)
            else:
                s = -math.log2(1 - p * c); ab += s
            tot += s
        genes = tsv(f'{bd}/{b}.tsv')
        spec = {l.strip() for l in open(f'{bd}/specific_genes.tsv') if l.strip()}
        # family prevalence: PanGBank vs survey build
        fam_k = {}
        for g in genes:
            fam_k[g['gene']] = None if g['gene'] in spec else len(fam_genomes.get(g['family'], ()))
        pb_fams = {g['family'] for g in genes if g['gene'] not in spec}
        pb_rare_fams = sum(1 for f in pb_fams if len(fam_genomes[f]) / n <= RARE)
        sv_mine = [f for f, gs in sv_fam.items() if b in gs]
        sv_rare = sum(1 for f in sv_mine if len(sv_fam[f] & set(sv_pub)) / len(sv_pub) <= RARE)
        sv_zero = sum(1 for f in sv_mine if not (sv_fam[f] & set(sv_pub)))
        # islands of specific genes
        bycontig = defaultdict(list)
        for g in genes: bycontig[g['contig']].append(g)
        persist_on = {ct: sum(1 for g in gs if g['partition'] == 'persistent' and g['gene'] not in spec)
                      for ct, gs in bycontig.items()}
        islands = []
        for ct, gs in bycontig.items():
            gs.sort(key=lambda g: int(g['start'])); run = []
            for g in gs + [None]:
                if g is not None and g['gene'] in spec:
                    run.append(g); continue
                if len(run) >= 3: islands.append((ct, run))
                run = []
        for ct, run in islands:
            island_rows.append(dict(tag=tag, bin=b, contig=ct, start=run[0]['start'], stop=run[-1]['stop'],
                                    genes=len(run), persistent_on_contig=persist_on[ct]))
        sb = sv_bins.get((short, b), {})
        # fate of the survey's rare-present modules
        coord2gene = {(g['contig'], int(g['start']), int(g['stop'])): g['gene'] for g in genes}
        svgff = gff_coords(f'{svd}/gff/{b}.gff') if os.path.exists(f'{svd}/gff/{b}.gff') else {}
        n_still = n_common = 0
        for fm in sv_flag.get((short, b), []):
            ks = []
            for gr in sv_genes.get((short, b, fm['module']), []):
                if gr['bin'] != b: continue
                co = svgff.get(gr['gene'])
                pg = coord2gene.get(co) if co else None
                ks.append('unmatched' if pg is None else fam_k[pg])
            matched = [k for k in ks if k != 'unmatched']
            nspec = sum(1 for k in matched if k is None)
            kfr = sorted((k or 0) / n for k in matched)
            med = kfr[len(kfr) // 2] if kfr else float('nan')
            fate = 'still rare' if (matched and med <= RARE) else ('now common' if matched else 'unmatched')
            n_still += fate == 'still rare'; n_common += fate == 'now common'
            fate_rows.append(dict(tag=tag, bin=b, module=fm['module'], survey_k=f"{fm['k']}/{fm['n']}",
                                  anchored=fm['anchored'], cls=fm['class'], genes=len(ks), matched=len(matched),
                                  specific=nspec, median_pangbank_freq=round(med, 3), fate=fate,
                                  products=fm['products'][:120]))
        per_bin.append(dict(
            tag=tag, bin=b, n_pangbank=n, n_survey_public=len(sv_pub), shared_public=overlap,
            compl_proj=round(c, 3), compl_survey=sb.get('completeness', ''),
            modules=len(mods), present=npres, rare_present=nrare,
            bits_total=round(tot, 1), bits_present=round(pb, 1), bits_absent=round(ab, 1),
            sv_modules=sb.get('modules', ''), sv_present=sb.get('present', ''), sv_rare_present=sb.get('rare_present', ''),
            sv_bits_total=sb.get('bits_total', ''), sv_bits_present=sb.get('bits_present', ''),
            sv_bits_absent=sb.get('bits_absent', ''),
            genes=len(genes), specific_genes=len(spec), specific_islands=len(islands), islands_anchored=sum(1 for ct, run in islands if persist_on[ct] > 0),
            island_genes=sum(len(run) for ct, run in islands),
            island_bits=round(len(islands) * math.log2(n + 2), 1),
            fam_rare_pb=pb_rare_fams + len(spec), fam_rare_sv=sv_rare, fam_zero_pb=len(spec), fam_zero_sv=sv_zero,
            sv_flagged=len(sv_flag.get((short, b), [])), sv_flagged_still_rare=n_still, sv_flagged_now_common=n_common))


def write(name, rows):
    if not rows: return
    with open(f'{OUT}/{name}', 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)


write('per_bin.tsv', per_bin); write('survey_module_fate.tsv', fate_rows); write('specific_islands.tsv', island_rows)
for r in per_bin:
    print('\t'.join(str(v) for v in r.values()))
