"""BGC Atlas R2.0 -> genome x gene-cluster-family (GCF) occurrence table for metapanclust's level-1 tests.

usage: build_presence.py ATLAS_DIR OUTDIR
writes OUTDIR/genome_gcf.parquet: one row per BGC in a genome (isolate, MAG or SAG) with its primary family,
  genome accession, source type, collection, CheckM2 completeness/contamination, GTDB lineage, region length,
  contig length, contig_edge and antiSMASH is_complete;
  OUTDIR/genomes.parquet (all genomes, also those without BGCs) and a summary on stdout.
"""
import sys, os
import pandas as pd, pyarrow.csv as pc

A, out = sys.argv[1:3]
os.makedirs(out, exist_ok=True)
db = f'{A}/db'


def read(t, cols):
    return pc.read_csv(f'{db}/{t}.tsv.gz', parse_options=pc.ParseOptions(delimiter='\t', newlines_in_values=True),
                       convert_options=pc.ConvertOptions(include_columns=cols)).to_pandas()


G = read('genomes', ['accession', 'assembly_id', 'source_type', 'collection', 'checkm_completeness',
                     'checkm_contamination', 'genome_quality', 'genome_size', 'domain', 'phylum', 'class', 'order',
                     'family', 'genus', 'species', 'bgc_count'])
B = read('bgcs', ['id', 'assembly', 'start_pos', 'end_pos', 'product_category', 'is_complete'])
F = read('bgc_features', ['bgc_id', 'contig_length', 'contig_edge'])
M = read('bgc_gcf', ['bgc_id', 'gcf_id', 'is_primary'])
M = M[M.is_primary.astype(str).isin(['t', 'True', 'true'])].drop(columns='is_primary')
X = (B.rename(columns={'id': 'bgc_id'}).merge(F, on='bgc_id', how='left').merge(M, on='bgc_id', how='left')
     .merge(G, left_on='assembly', right_on='assembly_id', how='inner'))
X['region_len'] = X.end_pos - X.start_pos
X = X.drop(columns=['assembly', 'assembly_id', 'start_pos', 'end_pos'])
X.to_parquet(f'{out}/genome_gcf.parquet', index=False)
G.to_parquet(f'{out}/genomes.parquet', index=False)

print(f'genomes {len(G):,} ({G.source_type.value_counts().to_dict()}); BGCs in genomes {len(X):,}; '
      f'in a family {X.gcf_id.notna().mean():.1%}')
edge = X.contig_edge.astype(str).isin(['t', 'True', 'true'])
print('BGC regions at a contig edge, by genome source:')
print(X.assign(edge=edge).groupby('source_type').agg(bgcs=('bgc_id', 'size'), at_edge=('edge', 'mean'),
      median_region_kb=('region_len', lambda v: v.median() / 1e3)).round(3).to_string())
hq = G[G.genome_quality == 'high']
sp = hq[hq.species.fillna('').str.len() > 3].groupby('species').size()
print(f'\nGTDB species with >= 15 / >= 50 high-quality genomes: {(sp >= 15).sum():,} / {(sp >= 50).sum():,}')
fam = X[X.gcf_id.notna() & X.species.fillna('').str.len().gt(3)]
per_sp = fam.groupby('species').gcf_id.nunique()
print(f'families per species with a family-assigned BGC: median {per_sp.median():.0f}, '
      f'species with >= 5 families: {(per_sp >= 5).sum():,}')
