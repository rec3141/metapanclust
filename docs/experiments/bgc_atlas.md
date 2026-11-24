# BGC Atlas v2 (R2.0) as a metapanclust and danaSeq resource (2026-11-24)

**Source:** BGC Atlas R2.0, final release of 22 Oct 2026. Data are served from `https://bgc-atlas.ziemertlab.com/downloads`; the
`bgc-atlas.cs.uni-tuebingen.de` host was not answering on 24 Nov. Preprint: Bağcı, Talamas-Tanner & Ziemert,
bioRxiv 2026.09.14.751540. Annotations are CC BY 4.0; the sequences keep their source licences. Site code:
`ZiemertLab/bgc-atlas-web-v2`.

**Downloaded** to `/data/scratch/metapanclust_data/bgc_atlas` (6.5 GB, SHA-256 verified;
`experiments/bgc_atlas/download.py`):
- every `db/*.tsv.gz` table (37);
- `bgc_uuid_source.parquet` (BGC → genome accession, contig, biome, GTDB lineage);
- both BiG-SLiCE 2.0 runs (`clustering-v3/`: family models for placing new BGCs; the per-gene tables are empty).

Not downloaded: the region sequences (GenBank 180 GB, nucleotide FASTA 178 GB, proteins 59 GB).

**Counts:** 16.1 M BGCs (2.7 M complete), 629,872 listed families (800,091 in all), 750,901 genomes (504k MAGs, 247k
isolates, 168 SAGs) with CheckM2 values and GTDB r232 lineages; genomes below MIMAG medium quality were removed.

## First look: `experiments/bgc_atlas/build_presence.py` → `derived/genome_gcf.parquet`

- 3.89 M BGCs sit in those genomes, and 95.3% of them have a primary family.
- **Fragmentation is the dominant MAG effect:**

  | genome source | BGCs | at a contig edge | median region |
  |---|---|---|---|
  | MAG | 1,986,083 | **80.4%** | 11.0 kb |
  | isolate | 1,900,667 | 34.6% | 21.6 kb |
  | SAG | 595 | 57.5% | 20.9 kb |

  Family absence in MAGs is therefore mostly a detection problem, which is the case the fragmentation-aware
  detection model (P_intact of a span on the MAG's contigs) was built for.
- **Level-1 test scale:** 1,987 GTDB species have ≥ 15 high-quality genomes, and 548 have ≥ 50. Across all
  species, the median number of families per species is 5, and 49,729 species have ≥ 5.

## Uses

- **metapanclust:**
  1. A level-1 held-out test with BGC families as units, on GTDB species with ≥ 15 genomes. Compare species
     prevalence, the gene-content and GTDB trees, and the effect of correcting MAG absences with the detection
     model (contig lengths are in `bgc_features`).
  2. Environment (biome) as a level-2 covariate.
  3. The site's family co-occurrence statistics (Jaccard, observed/expected) as a baseline for context.
     These are only in their database, not in the downloads.
- **danaSeq:** assign antiSMASH 8 regions to Atlas families locally with `bigslice --query` against the two
  published runs. That gives a known-or-novel call, the family's taxonomic range and its biome distribution,
  for the viz.
