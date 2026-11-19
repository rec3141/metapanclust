"""LOCAL (run once after p1_rep_ani.sh): choose focal species and their contaminants for the grex run.

Focal species: >= NTRAIN + NHELD HQ isolate genomes and at least one congener at 88-97% ANI
(species representatives, skani). At most MAXGENUS focal species per genus, NFOCAL in all
(genera with more candidates first get their best-supplied species; ties random, seeded).
Contaminants per focal species, 2 genomes each:
  close   up to 3 congeners at 88-97% ANI (highest, lowest, median)
  mid     up to 2 congeners at 80-88% ANI
  distant 1 congener with no skani ANI (below ~80%), if any
  far     Flavobacterium psychrophilum (Escherichia coli for Bacteroidota focal species)
usage: p2_select.py PLAN_DIR   (reads hq_iso.tsv, sp.tsv, rep_ani.tsv; writes plan.tsv, pairs_planned.tsv)
plan.tsv: acc  species  genus  role(train|held|contam|far)  tag  rep_ani(contaminants)
"""
import random, re, sys
from collections import defaultdict
import pandas as pd

D = sys.argv[1]
NTRAIN, NHELD, NCON, NFOCAL, MAXGENUS = 15, 6, 2, 200, 3
rng = random.Random(2026)
hq = pd.read_csv(f'{D}/hq_iso.tsv', sep='\t', header=None, names=['acc', 'rep', 'tax', 'contigs'])
hq['species'] = hq.tax.str.split(';').str[-1].str[3:]
hq['genus'] = hq.tax.str.split(';').str[-2].str[3:]
hq['phylum'] = hq.tax.str.split(';').str[1].str[3:]
by_sp = {s: list(zip(d.contigs, d.acc)) for s, d in hq.groupby('species')}
n_sp = {s: len(v) for s, v in by_sp.items()}
phylum = dict(zip(hq.species, hq.phylum))
ani = pd.read_csv(f'{D}/rep_ani.tsv', sep='\t', header=None, names=['genus', 'focal', 'other', 'ani', 'af_f', 'af_o'])
for c in ('genus', 'focal', 'other'):
    ani[c] = ani[c].str[3:]
ani = ani.groupby(['genus', 'focal', 'other'], as_index=False).ani.max()
need = NTRAIN + NHELD
ani = ani[ani.focal.map(n_sp).fillna(0).ge(need) & ani.other.map(n_sp).fillna(0).ge(NCON)]
close = ani[(ani.ani >= 88) & (ani.ani <= 97)]
cand = close.groupby(['genus', 'focal']).size().rename('n_close').reset_index()
cand['n'] = cand.focal.map(n_sp)
cand['r'] = [rng.random() for _ in range(len(cand))]
# round-robin over genera so no genus dominates: rank within genus by n_close then genome count
cand = cand.sort_values(['genus', 'n_close', 'n', 'r'], ascending=[True, False, False, True])
cand['k'] = cand.groupby('genus').cumcount()
cand = cand[cand.k < MAXGENUS].sort_values(['k', 'r'])
focal = cand.head(NFOCAL)
print(f'{len(cand)} candidate focal species in {cand.genus.nunique()} genera; using {len(focal)} '
      f'in {focal.genus.nunique()} genera')


def pick(sp, n):
    gs = list(by_sp[sp]); rng.shuffle(gs)
    gs = sorted(gs, key=lambda x: x[0])[:3 * n]         # random among the 3n least fragmented
    rng.shuffle(gs)
    return [a for _, a in gs[:n]]


def tag(sp):
    return re.sub(r'[^A-Za-z0-9]+', '_', sp).strip('_')


genus_species = defaultdict(list)
for s, n in n_sp.items():
    if n >= NCON:
        genus_species[s.split(' ')[0]].append(s)   # GTDB genus = first word of the species name
rows, pairs = [], []
for f in focal.itertuples():
    t = tag(f.focal)
    for i, a in enumerate(pick(f.focal, need)):
        rows.append((a, f.focal, f.genus, 'held' if i < NHELD else 'train', t, ''))
    a_f = ani[ani.focal == f.focal].sort_values('ani', ascending=False)
    cl = a_f[(a_f.ani >= 88) & (a_f.ani <= 97)]
    chosen = []
    if len(cl):
        idx = sorted({0, len(cl) - 1, len(cl) // 2})
        chosen += [(r.other, r.ani, 'close') for r in cl.iloc[idx].itertuples()]
    mid = a_f[(a_f.ani >= 80) & (a_f.ani < 88)]
    if len(mid):
        idx = sorted({0, len(mid) - 1})
        chosen += [(r.other, r.ani, 'mid') for r in mid.iloc[idx].itertuples()]
    far_cong = [s for s in genus_species[f.focal.split(' ')[0]] if s != f.focal and s not in set(a_f.other)]
    if far_cong:
        chosen.append((rng.choice(sorted(far_cong)), float('nan'), 'distant'))
    for s, a, cls in chosen:
        for acc in pick(s, NCON):
            rows.append((acc, s, f.genus, 'contam', t, a))
        pairs.append((t, f.focal, s, a, cls))
    far = 'Escherichia coli' if phylum[f.focal] == 'Bacteroidota' else 'Flavobacterium psychrophilum'
    for acc in pick(far, NCON):
        rows.append((acc, far, far.split(' ')[0], 'far', t, ''))
plan = pd.DataFrame(rows, columns=['acc', 'species', 'genus', 'role', 'tag', 'rep_ani'])
plan.to_csv(f'{D}/plan.tsv', sep='\t', index=False, header=False)
pd.DataFrame(pairs, columns=['tag', 'focal', 'contaminant', 'rep_ani', 'class']).to_csv(
    f'{D}/pairs_planned.tsv', sep='\t', index=False)
P = pd.DataFrame(pairs, columns=['tag', 'focal', 'contaminant', 'rep_ani', 'class'])
print(f'{len(plan)} rows, {plan.acc.nunique()} distinct genomes; pairs by class: {P["class"].value_counts().to_dict()}')
print('close pairs by ANI:', pd.cut(P.rep_ani, [88, 90, 92, 94, 95, 96, 97.01], right=False).value_counts().sort_index().to_dict())
