"""CAMI II marine source genomes -> GTDB species (GTDB-Tk) -> PanGBank GTDB_all 1.0.0 (r226) pangenome.

usage: map_species.py GTDBTK_OUT_DIR CATALOG.tsv HAVE_IDS.txt OUTDIR
writes OUTDIR/genome_species.tsv (genome, gtdb classification, species, pangbank id, genomes in pangenome) and
OUTDIR/downloads.tsv (id, name, genomes, md5) for pangenomes not already downloaded.
The GTDB-Tk database may be a later release than PanGBank's r226; species whose names changed are unmatched.
"""
import glob, os, sys
import pandas as pd

gdir, cat_f, have_f, out = sys.argv[1:5]
os.makedirs(out, exist_ok=True)
G = pd.concat([pd.read_csv(f, sep='\t') for f in glob.glob(f'{gdir}/**/gtdbtk.*.summary.tsv', recursive=True)])
G['species'] = G.classification.str.split(';s__').str[-1]
G.loc[G.species.isin(['', 'Unclassified']) | ~G.classification.str.contains(';s__'), 'species'] = ''
cat = pd.read_csv(cat_f, sep='\t')
sp2 = cat.set_index('species')
G['pangbank_id'] = G.species.map(sp2.id).astype('Int64')
G['pangbank_genomes'] = G.species.map(sp2.genomes).astype('Int64')
G = G.rename(columns={'user_genome': 'genome'})
G[['genome', 'classification', 'species', 'closest_genome_ani' if 'closest_genome_ani' in G else 'fastani_ani',
   'pangbank_id', 'pangbank_genomes']].to_csv(f'{out}/genome_species.tsv', sep='\t', index=False)
have = set(int(x) for x in open(have_f).read().split()) if os.path.exists(have_f) else set()
need = cat[cat.id.isin(G.pangbank_id.dropna().astype(int)) & ~cat.id.isin(have)]
need[['id', 'name', 'genomes', 'md5']].to_csv(f'{out}/downloads.tsv', sep='\t', index=False, header=False)
print(f'{len(G)} genomes; {(G.species != "").sum()} with a GTDB species; {G.pangbank_id.notna().sum()} with a '
      f'PanGBank pangenome ({G.pangbank_id.nunique()} pangenomes, {len(need)} to download, '
      f'{need.genomes.sum()} genomes)')
