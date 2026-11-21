# CAMI II marine: metapanclust v0 against CheckM2 on real binner output (2026-11-21)

Scripts: `experiments/cami_marine/` (prepare, write_bins, map_species, run_grex.sbatch, analyze). grex jobs
8195137 (GTDB-Tk on the 777 source genomes), 8195223 (CheckM2 on all 7,018 bins), 8195229 (scoring). Full
tables: [cami_marine_results.txt](cami_marine_results.txt).

## Setup

- **Bins**: all 15 CAMI II submissions on the short-read pooled gold-standard assembly (MetaBAT, MaxBin,
  CONCOCT, VAMB, MetaBinner, MetaWRAP, Autometa, UltraBinner…). Bins of at least 200 kb whose majority source is
  a genome (not a plasmid or virus) were kept: 7,018. True completeness and contamination are by bp against the
  gold-standard contig-to-genome mapping (AMBER's definitions).
- **Reference per bin**: the PanGBank pangenome of the GTDB species of the bin's majority source genome. This is
  an oracle assignment; in practice the species would come from GTDB-Tk on the bin. Only 205 of 777 source
  genomes (95 species) have a PanGBank pangenome (≥ 15 genomes per species), so 1,748 bins were scorable.
  1,567 were scored; the other 181 belong to species with 1,200-2,000 reference genomes and hit the 1 h limit.
- **CheckM2** 1.x (danaSeq container) on every bin.

Contamination in scored bins: 401 clean (< 1%), 679 mainly from other species, 487 mainly from another genome
of the same species (CAMI marine has 36 named strain groups).

## Results

**1. Flagging contaminated bins: CheckM2 wins, and metapanclust adds almost nothing.** Bins ≥ 50% complete,
n = 987:

| true contamination ≥ | CheckM2 | novel-gene bits | Z+ | CheckM2 + metapanclust (CV) | CheckM2 alone (CV) |
|---|---|---|---|---|---|
| 2% | 0.922 | 0.830 | 0.669 | 0.924 | 0.918 |
| 5% | 0.965 | 0.846 | 0.683 | 0.970 | 0.963 |
| 10% | 0.970 | 0.861 | 0.710 | 0.975 | 0.966 |

Novelty does better on other-species contamination (AUC 0.92) than on same-species (0.83). CheckM2 reaches
0.98-1.00 on both. Spearman with true contamination: CheckM2 0.88, novel bits 0.66, Z+ 0.37.

**2. Completeness: metapanclust is more accurate than CheckM2**, given the right species pangenome. Mean absolute
error is 5.2 points with no bias (Spearman 0.95), against 8.2 points for CheckM2, which overestimates by +6.2
(Spearman 0.91). The estimate is just the detected share of the species' persistent families, weighted by how
often each is carried.

**3. Finding the foreign contigs: this is what CheckM2 cannot do, and metapanclust can.** In bins with ≥ 5%
contamination:

| foreign source | AUC novel fraction | AUC divergence | AUC unanchored | foreign bp on novel-majority contigs | native bp on them |
|---|---|---|---|---|---|
| other species | 0.844 | 0.868 | 0.787 | 82% | 5.3% |
| same species | 0.464 | 0.473 | 0.472 | 7% | 5.3% |

**4. Cleaning by novelty.** Removing contigs whose genes are mostly (> 50%) novel to the species:

| bins | median contamination before → after | mean contamination drop | mean completeness lost |
|---|---|---|---|
| other-species contaminated (679) | 13.0% → 0.5% | 19.6 points | 3.2 points |
| same-species contaminated (487) | 45.3% → 37.8% | 5.6 points | 2.1 points |
| clean (401) | 0.05% → 0% | 0.15 points | 1.0 points |

Requiring the contig to also carry no persistent gene (unanchored) is more conservative: 13.0% → 1.4%, at a cost of
2.2 completeness points. (Contigs without a called gene are left in the bin; they are short and few.)

**5. Unit-level surprise (Z+, rare-present) tracks contamination only weakly** (AUC 0.67-0.71). It rises with
binner contamination (CONCOCT median Z+ +1.6 to +2.9, clean binners −3 to −4.5), but its calibration
offset (scorer_v0.md) makes the absolute values uninterpretable.

## What this means

- **Do not sell surprise as a contamination detector.** CheckM2 already does that job at AUC 0.97.
- **Use the species pangenome for what CheckM2 cannot do:**
  - say *which contigs* are foreign: 82% of other-species contamination sits on contigs that are mostly novel
    to the species;
  - clean them at about 3% completeness cost;
  - estimate completeness without CheckM2's +6-point inflation.
  
  This becomes stage 5 (binning moves): evict contigs novel to the bin's species, keep the rest.
- **Same-species contamination (strain mixtures) is invisible to novelty**, as expected, since strains share the
  species' gene pool. It needs coverage (per-sample depth differs between strains) or per-gene identity (the
  merge_identity test found median per-gene divergence separates strains from 94-95% ANI congeners, not strains
  from each other).
- **Coverage of the method is the bottleneck.** PanGBank covers 95 of 747 classified CAMI species. The
  remaining 78% of bins need a genus-level reference: a pooled pangenome of the GTDB genus, where
  xsp_tree_vs_content says the GTDB tree is the right relatedness. Large species (≥ 1,200 genomes) need
  subsampling to fit the time limit.
