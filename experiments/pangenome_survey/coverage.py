"""Per-sample depth of rare-present module contigs against the bin's core (main mag_chesterfield bins only;
depths from assembly_chesterfield/mapping/depths.txt, MetaBAT jgi format, one column per sample).
Core profile: length-weighted mean depth of the bin's contigs carrying >= 3 persistent genes. For each
module: length-weighted mean depth of its contigs; Pearson r of log2(depth + 0.1) with the core over
samples, and the module/core ratio of depth summed over samples.
Calibration: the same two numbers for each core contig against the core built from the others.
usage: coverage.py OUTDIR SURVEY_SPECIES_TSV"""
import csv, math, os, statistics as st, sys
from collections import defaultdict, Counter
o, sp_tsv = sys.argv[1:3]
DEP = '/data/scratch/metapanclust_data/chesterfield/assembly_chesterfield/mapping/depths.txt'
dirs = {l.split('\t')[0]: l.rstrip('\n').split('\t')[2] for l in open(sp_tsv) if l.strip()}
genes = list(csv.DictReader(open(f'{o}/flagged_genes.tsv'), delimiter='\t'))
mods = list(csv.DictReader(open(f'{o}/flagged_modules.tsv'), delimiter='\t'))
main_bins = sorted({(r['species'], r['bin']) for r in mods if not r['bin'].startswith(('chesterfield_marine', 'chesterfield_fresh'))})

# persistent genes per contig for each main bin
core_contigs = {}
for sp, b in main_bins:
    cnt = Counter()
    for l in open(f"{dirs[sp]}/gff/{b}.gff"):
        if '\tCDS\t' in l and 'partition=persistent' in l: cnt[l.split('\t')[0]] += 1
    core_contigs[(sp, b)] = {c for c, n in cnt.items() if n >= 3}
need = set().union(*core_contigs.values()) | {r['contig'] for r in genes}
dep, ln = {}, {}
with open(DEP) as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    idx = [i for i, h in enumerate(hdr) if h.endswith('.bam')]
    for l in fh:
        f = l.rstrip('\n').split('\t')
        if f[0] in need:
            dep[f[0]] = [float(f[i]) for i in idx]; ln[f[0]] = int(f[1])


def prof(cs):
    cs = [c for c in cs if c in dep]
    if not cs: return None
    L = sum(ln[c] for c in cs)
    return [sum(dep[c][i] * ln[c] for c in cs) / L for i in range(len(idx))]


def compare(m, core):
    x = [math.log2(v + 0.1) for v in m]; y = [math.log2(v + 0.1) for v in core]
    try: r = st.correlation(x, y)
    except st.StatisticsError: r = float('nan')
    # ratio over summed depth: most bins reach >= 1x in only a few samples
    return r, (sum(m) / sum(core) if sum(core) > 0 else float('nan')), sum(1 for c in core if c >= 0.5)


rows, calib = [], []
mod_contigs = defaultdict(set)
for r in genes:
    mod_contigs[(r['species'], r['bin'], r['module'])].add(r['contig'])
for sp, b in main_bins:
    cc = sorted(core_contigs[(sp, b)])
    core = prof(cc)
    if core is None: continue
    for c in cc:                      # calibration: leave-one-out core contigs
        rest = prof([x for x in cc if x != c])
        if rest is None or c not in dep: continue
        r, rat, n = compare(dep[c], rest)
        calib.append(dict(species=sp, bin=b, contig=c, r=r, ratio=rat))
    for m in mods:
        if m['species'] != sp or m['bin'] != b: continue
        cs = mod_contigs[(sp, b, m['module'])]
        p = prof(cs)
        if p is None: continue
        r, rat, n = compare(p, core)
        rows.append(dict(species=sp, bin=b, module=m['module'], families=m['families'], k=m['k'], n=m['n'],
                         cls=m['class'], anchored=m['anchored'], contigs=len(cs), r=round(r, 3), ratio=round(rat, 3),
                         samples_core_ge05x=n, core_total_depth=round(sum(core), 2)))
with open(f'{o}/module_coverage.tsv', 'w') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)
print('| Species | Bin | Core contigs: r median (5th pct) | ratio median (5-95th) | Modules | r < core 5th pct | ratio < 0.5 | ratio > 2 |')
print('|---|---|---|---|---|---|---|---|')
for sp, b in main_bins:
    cr = sorted(x['r'] for x in calib if x['bin'] == b and x['species'] == sp and not math.isnan(x['r']))
    ca = sorted(x['ratio'] for x in calib if x['bin'] == b and x['species'] == sp and not math.isnan(x['ratio']))
    mr = [x for x in rows if x['bin'] == b and x['species'] == sp]
    if not cr or not mr: continue
    q = lambda v, p: v[min(len(v) - 1, int(p * len(v)))]
    lo = q(cr, 0.05)
    print(f"| {sp} | {b.replace('chesterfield_', '')} | {st.median(cr):.2f} ({lo:.2f}) | {st.median(ca):.2f} ({q(ca, .05):.2f}-{q(ca, .95):.2f}) | {len(mr)} | "
          f"{sum(1 for x in mr if x['r'] < lo)} | {sum(1 for x in mr if x['ratio'] < 0.5)} | {sum(1 for x in mr if x['ratio'] > 2)} |")
