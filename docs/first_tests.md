# First tests (2026-11-18)

Data stay outside the repo (workstation: `/data/scratch/metapanclust_data/`).

## Tools

- **Spacedust** 2d68c85 (static build; the bioconda 2.e56c505 build crashes in `clusterhits` with
  "Invalid query lookup record"). The static binary needs glibc >= 2.34 and libatomic, so it runs in a
  small Ubuntu 24.04 container. `scripts/run_spacedust.sh`, MMseqs2 search mode for now.
- **PPanGGOLiN** 2.3.2 (`envs/ppanggolin.yml`): `ppanggolin all` (families, partitions, RGPs, spots,
  modules).

## Operons in E. coli K-12 (curated EcoCyc transcription units)

Adjacent same-strand gene pairs, 2,820 labelled (1,727 in the same transcription unit).
`scripts/operon_pairs.py`, `scripts/eval_pairs.py`.

| Rule | Accuracy | Precision | Recall |
|---|---|---|---|
| gap <= 50 bp | 0.850 | 0.94 | 0.81 |
| Spacedust conserved (6 distant genomes) | 0.593 | 0.90 | 0.38 |
| gap <= 50 OR (conserved AND gap <= 300) | **0.874** | 0.92 | **0.87** |

Conservation across only six distant genomes (E. coli, B. subtilis, P. aeruginosa, M. tuberculosis,
Synechocystis, M. jannaschii) already recovers 117 of the 331 far-apart operon pairs. Closer relatives
should do more.

## Chesterfield co-assembly, 93 DAS Tool bins

Spacedust all-vs-all: 4.5 min, 264,176 pairwise conserved clusters. Median 71% of a bin's genes are
in at least one (range 12-96%; Spacedust paper: 58% over 1,308 reference genomes).

First catalogue (`scripts/sd_families.py`): merging overlapping pairwise clusters into loci and linking
loci by single-linkage collapses into one giant family (14,709 of 15,372 loci, all 93 bins). Next:
represent each locus by its gene-family content and group by content similarity (Panorama/BiG-SCAPE
style), or community detection on the locus graph instead of connected components.

## Species pangenomes and module surprise

Two Chesterfield species with many GTDB r226 genomes (all MAGs), plus our bins from the Chesterfield
and Chesterfield-marine assemblies:

| Species | Genomes | Gene families | Persistent / shell / cloud | Modules |
|---|---|---|---|---|
| MAG-120531 sp000173115 (Flavobacteriaceae) | 67 + 2 | 4,394 | 1,472 / 267 / 2,655 | 41 |
| UBA9145 sp001438145 (Pseudohongiellaceae) | 52 + 2 | 5,878 | 1,818 / 408 / 3,652 | 105 |

Prototype surprise (`scripts/module_surprise.py`): module prevalence p among the public genomes
(Beta(1,1)); present -log2 p, absent -log2(1 - p c), c = bin completeness from persistent families.

- UBA9145: both Chesterfield bins carry 36 modules found in none of the 52 public genomes (5.75 bits
  each), led by two flagellar gene sets (FliA sigma factor, flagellar system and basal-body rod
  proteins); none on geNomad plasmid/virus contigs. Either a motile local population or contamination
  shared by the two overlapping co-assemblies: the next check is per-sample coverage of these modules
  against the bin's core.
- MAG-120531: one 4-family module (hypothetical proteins) in both bins and no public genome.
