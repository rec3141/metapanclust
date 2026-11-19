"""Pick genomes from GTDB r226 (HQ isolates: CheckM2 >= 95% complete, <= 2% contaminated, <= 200 contigs).
Per genus: 2 focal species (NTRAIN training + NHELD held-out genomes each) and 4 extra congeneric
species (2 genomes each) used only as contaminants; plus a far contaminant (another phylum).
usage: 00_select.py HQ_TSV OUTDIR   (HQ_TSV: accession, taxonomy, compl, contam, contigs, size, rep, category)
writes OUTDIR/genomes.tsv: accession  species  genus  role(train|held|contam|far)  focal_tag
"""
import random, sys
from collections import defaultdict

HQ, OUT = sys.argv[1:3]
NTRAIN, NHELD, NCON = 15, 4, 2
PLAN = {  # genus: (focal species, extra contaminant species)
    'Sulfitobacter': (['mediterraneus', 'pontiacus'], ['mediterraneus_A', 'geojensis', 'faviae', 'dubius']),
    'Pseudoalteromonas': (['distincta', 'arctica'], ['atlantica', 'tetraodonis', 'nigrifaciens', 'ruthenica']),
    'Shewanella': (['baltica', 'putrefaciens'], ['hafniensis', 'oncorhynchi', 'vesiculosa', 'algae']),
    'Alteromonas': (['macleodii', 'mediterranea'], ['stellipolaris', 'australica', 'naphthalenivorans', 'hispanica']),
}
FAR = 'Flavobacterium psychrophilum'
rng = random.Random(7)
by_sp = defaultdict(list)
for line in open(HQ):
    acc, tax, comp, cont, ncontig, size, rep, cat = line.rstrip('\n').split('\t')
    if cat != 'none':  # isolates only
        continue
    sp = tax.split(';')[-1][3:]
    by_sp[sp].append((int(ncontig), acc))
rows = []


def pick(sp, n):
    gs = sorted(by_sp[sp]); rng.shuffle(gs)
    gs = sorted(gs[:max(n, len(gs))], key=lambda x: x[0])  # all candidates, fewest contigs first
    sel = [a for _, a in gs[:3 * n]]; rng.shuffle(sel)     # random among the 3n least fragmented
    assert len(sel) >= n, (sp, len(sel))
    return sel[:n]


for genus, (focal, extra) in PLAN.items():
    for s in focal:
        sp = f'{genus} {s}'
        sel = pick(sp, NTRAIN + NHELD)
        tag = s.lower()
        for i, a in enumerate(sel):
            rows.append((a, sp, genus, 'held' if i < NHELD else 'train', tag))
    for s in extra:
        sp = f'{genus} {s}'
        for a in pick(sp, NCON):
            rows.append((a, sp, genus, 'contam', ''))
for a in pick(FAR, NCON):
    rows.append((a, FAR, 'Flavobacterium', 'far', ''))
with open(f'{OUT}/genomes.tsv', 'w') as f:
    for r in rows:
        f.write('\t'.join(r) + '\n')
print(len(rows), 'genomes')
