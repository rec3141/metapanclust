"""Species for the tree-vs-context run: a stratified subset of part A's species (pangbank_scale sel/species_A.tsv),
whose .h5 files are already downloaded.

usage: select_species.py species_A.tsv OUTDIR [PER_STRATUM=40] [SEED=1]
Strata by genome count: 25-49, 50-99, 100-199, 200-399, >= 400 (species with < 25 genomes are dropped; the pool
is capped at 400 genomes, so the last stratum costs the same as 200-399 apart from reading the .h5).
Within a stratum, species are taken round-robin over phyla (random order inside a phylum), so phyla stay mixed.
Writes OUTDIR/species_tvc.tsv (catalogue rows + stratum) and OUTDIR/tasks_tvc.tsv (one id per line, most
expensive first).
"""
import sys, os
import numpy as np, pandas as pd

src, out = sys.argv[1], sys.argv[2]
per = int(sys.argv[3]) if len(sys.argv) > 3 else 40
seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
rng = np.random.default_rng(seed)
s = pd.read_csv(src, sep='\t')
s = s[s.genomes >= 25].copy()
s['stratum'] = pd.cut(s.genomes, [24, 49, 99, 199, 399, 10 ** 6], labels=['25-49', '50-99', '100-199', '200-399', '400+'])
pick = []
for st, g in s.groupby('stratum', observed=True):
    by = {ph: list(rng.permutation(x.index)) for ph, x in g.groupby('phylum')}
    order = list(rng.permutation(sorted(by)))
    chosen = []
    while len(chosen) < per and any(by.values()):
        for ph in order:
            if by[ph] and len(chosen) < per: chosen.append(by[ph].pop())
    pick += chosen
S = s.loc[pick].copy()
S['cost'] = np.minimum(S.genomes, 400) ** 2 + S.genomes * 50          # skani pairs + .h5 read
S = S.sort_values('cost', ascending=False)
os.makedirs(out, exist_ok=True)
S.drop(columns='cost').to_csv(f'{out}/species_tvc.tsv', sep='\t', index=False)
S[['id']].to_csv(f'{out}/tasks_tvc.tsv', sep='\t', index=False, header=False)
print(f'{len(S)} species, {S.phylum.nunique()} phyla, {S.genomes.sum()} genomes '
      f'({np.minimum(S.genomes, 400).sum()} in pools)')
print(S.groupby('stratum', observed=True).agg(n=('id', 'size'), phyla=('phylum', 'nunique')).to_string())
