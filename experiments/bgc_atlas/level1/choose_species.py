"""Select GTDB species for the BGC-family level-1 test and extract their Atlas data.

Species must have >= MIN_G Atlas genomes from the GTDB collection (MIMAG medium or high), >= MIN_F distinct primary
families among them, and a PanGBank (GTDB_all 1.0.0, r226) pangenome of 20..MAX_PB genomes. N species are drawn
round-robin over phyla (largest first within a phylum, seed 1).

usage: choose_species.py ATLAS_DERIVED_DIR CATALOG.tsv OUTDIR [N=60] [MIN_G=20] [MIN_F=10] [MAX_PB=400]
writes OUTDIR/species.tsv, OUTDIR/downloads.tsv (PanGBank id, name, genomes, md5), and per species
OUTDIR/atlas/<pangbank_id>.tsv: genome (PanGBank-style accession), completeness, contamination, quality, source_type,
gcf_id, contig, start, end, contig_edge, is_complete (one row per BGC; genomes without BGCs get one row with empty gcf).
"""
import os, sys, random
import pandas as pd

D, cat_f, out = sys.argv[1:4]
N = int(sys.argv[4]) if len(sys.argv) > 4 else 60
MIN_G = int(sys.argv[5]) if len(sys.argv) > 5 else 20
MIN_F = int(sys.argv[6]) if len(sys.argv) > 6 else 10
MAX_PB = int(sys.argv[7]) if len(sys.argv) > 7 else 400
os.makedirs(f'{out}/atlas', exist_ok=True)
G = pd.read_parquet(f'{D}/genomes.parquet')
G = G[(G.collection == 'GTDB') & G.genome_quality.isin(['high', 'medium'])].copy()
G['genome'] = G.accession.str.replace(r'^(GB|RS)_', '', regex=True)
G['sp'] = G.species.fillna('').str.replace('s__', '', regex=False)
X = pd.read_parquet(f'{D}/genome_gcf.parquet')
X = X[X.accession.isin(set(G.accession))].copy()
X['genome'] = X.accession.str.replace(r'^(GB|RS)_', '', regex=True)
cat = pd.read_csv(cat_f, sep='\t')
cat = cat[(cat.genomes >= 20) & (cat.genomes <= MAX_PB)].set_index('species')
ng = G.groupby('sp').size()
nf = X.dropna(subset=['gcf_id']).assign(sp=lambda d: d.species.str.replace('s__', '', regex=False)).groupby('sp').gcf_id.nunique()
cand = pd.DataFrame({'n_genomes': ng, 'n_fams': nf}).dropna()
cand = cand[(cand.n_genomes >= MIN_G) & (cand.n_fams >= MIN_F) & cand.index.isin(cat.index)]
cand = cand.join(cat[['id', 'name', 'genomes', 'md5', 'phylum']])
by = {p: d.sort_values('n_genomes', ascending=False).index.tolist() for p, d in cand.groupby('phylum')}
rng = random.Random(1)
order = sorted(by, key=lambda p: -len(by[p]))
sel = []
while len(sel) < N and any(by.values()):
    for p in order:
        if by[p] and len(sel) < N:
            sel.append(by[p].pop(0))
S = cand.loc[sel].rename_axis('species').reset_index()
S.to_csv(f'{out}/species.tsv', sep='\t', index=False)
S[['id', 'name', 'genomes', 'md5']].to_csv(f'{out}/downloads.tsv', sep='\t', index=False, header=False)
for _, r in S.iterrows():
    g = G[G.sp == r.species][['genome', 'checkm_completeness', 'checkm_contamination', 'genome_quality', 'source_type']]
    x = X[X.species.str.replace('s__', '', regex=False) == r.species]
    # contig names and coordinates come from the region file name: <genome>_<contig>.regionNNN; use bgcs columns
    x = x[['genome', 'gcf_id', 'bgc_id', 'contig', 'start_pos', 'end_pos', 'contig_edge', 'is_complete', 'region_len']]
    t = g.merge(x, on='genome', how='left')
    t.to_csv(f'{out}/atlas/{r.id}.tsv', sep='\t', index=False)
print(f'{len(cand)} candidate species; selected {len(S)} over {S.phylum.nunique()} phyla; '
      f'{int(S.n_genomes.sum())} Atlas genomes; PanGBank genomes to download {int(S.genomes.sum())}')
print(S[['species', 'phylum', 'n_genomes', 'n_fams', 'genomes']].to_string(index=False))
