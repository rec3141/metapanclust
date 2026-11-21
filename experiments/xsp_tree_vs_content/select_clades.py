"""Between-species tree vs gene content: pick clades of PanGBank species (GTDB_all 1.0.0 = GTDB r226).

Clades are GTDB families with >= 20 PanGBank species, spread over phyla and chosen to include GC-diverse and
marine groups, plus the archaeal class Nitrososphaeria (incl. Nitrosopumilaceae = Marine Group I). Up to
MAX_SP species per clade, drawn at random among species with <= MAXG genomes (download size).

usage: select_clades.py CATALOG.tsv HAVE_DOWNLOADS.tsv OUTDIR [MAX_SP=40] [MAXG=250] [SEED=1]
writes OUTDIR/clades.tsv (clade, id, name, genomes, genus, species) and OUTDIR/downloads.tsv (id, name, genomes,
md5) for species whose .h5 is not in HAVE_DOWNLOADS.
"""
import sys, os
import pandas as pd

cat, have_f, out = sys.argv[1:4]
MAX_SP = int(sys.argv[4]) if len(sys.argv) > 4 else 40
MAXG = int(sys.argv[5]) if len(sys.argv) > 5 else 250
SEED = int(sys.argv[6]) if len(sys.argv) > 6 else 1
os.makedirs(out, exist_ok=True)

CLADES = [  # (rank, name)
    ('family', 'Rhodobacteraceae'), ('family', 'Flavobacteriaceae'), ('family', 'Burkholderiaceae'),
    ('family', 'Mycobacteriaceae'), ('family', 'Streptomycetaceae'), ('family', 'Lactobacillaceae'),
    ('family', 'Vibrionaceae'), ('family', 'Pseudomonadaceae'), ('family', 'Helicobacteraceae'),
    ('family', 'Campylobacteraceae'), ('family', 'Rikenellaceae'), ('family', 'Lachnospiraceae'),
    ('family', 'Sphingomonadaceae'), ('family', 'Alteromonadaceae'), ('family', 'Gastranaerophilaceae'),
    ('family', 'Bifidobacteriaceae'), ('family', 'Enterobacteriaceae'), ('class', 'Nitrososphaeria'),
]
c = pd.read_csv(cat, sep='\t')
have = set(pd.read_csv(have_f, sep='\t', header=None)[0])
rows = []
for rank, name in CLADES:
    sub = c[(c[rank] == name) & (c.genomes <= MAXG)]
    if len(sub) > MAX_SP:
        # keep already-downloaded species first (no cost), then a random draw of the rest
        h = sub[sub.id.isin(have)]
        rest = sub[~sub.id.isin(have)].sample(frac=1, random_state=SEED)
        sub = pd.concat([h.sample(frac=1, random_state=SEED), rest]).head(MAX_SP)
    for _, r in sub.iterrows():
        rows.append(dict(clade=name, id=r.id, name=r['name'], genomes=r.genomes, genus=r.genus, species=r.species,
                         md5=r.md5))
T = pd.DataFrame(rows)
T.drop(columns='md5').to_csv(f'{out}/clades.tsv', sep='\t', index=False)
D = T[~T.id.isin(have)][['id', 'name', 'genomes', 'md5']].drop_duplicates('id')
D.to_csv(f'{out}/downloads.tsv', sep='\t', index=False, header=False)
print(T.groupby('clade').agg(n=('id', 'size'), genera=('genus', 'nunique'), genomes=('genomes', 'sum')).to_string())
print(f'{len(T)} species, {T.id.isin(have).sum()} already on grex, {len(D)} to download '
      f'({D.genomes.sum()} genomes)')
