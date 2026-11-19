# Conserved-gene-cluster catalogue, Chesterfield 93 bins (2026-11-18)

Stage 1 of the plan. Code: `experiments/catalogue/` (`run_all.sh` rebuilds everything in ~1 min, plus
~10 min for the MMseqs2 families on 2 threads). Data: `/data/scratch/metapanclust_data/catalogue/`.

## Result

The giant family is gone. 20,849 conserved loci fall into **7,178 cluster families**: 3,057 with two or
more loci (16,728 loci), 2,976 spanning at least two genera. The largest family has 98 loci in 85 bins
(the S10-spc-alpha ribosomal protein superoperon), not 14,709. The top 40 families are textbook operons
(table below). Run on the same locus graph, connected components still give one component of 12,053
loci, so the fix comes from Leiden, not from the edge filter. Two Leiden seeds agree at ARI 0.975 for
families and 0.94 for subfamilies.

| | first attempt (`scripts/sd_families.py`) | this catalogue |
|---|---|---|
| loci | 15,372 | 20,849 |
| genes per locus, median / 90% / 99% / max | 3 / 15 / 58 / 403 | 3 / 9 / 19 / 47 |
| genes in a locus | 109,684 | 93,937 (54% of 173,368) |
| families | 663 (one with 14,709 loci) | 7,178 (largest 98 loci) |

## Method

1. **Gene families** (`01_genes_families.py`): MMseqs2 `easy-cluster` over the 173,368 binned proteins
   only, 30% identity, 80% coverage, `-s 7.5`. Result: 94,987 families, most of them singletons in this
   diverse community.
2. **Loci** (`02_loci.py`): a *link* joins two consecutive genes. A pairwise Spacedust cluster supports
   every link between two of its members that are at most one gene apart. Support counts **distinct
   partner genera other than the bin's own genus** (GTDB-Tk). Two bins of the same genus share synteny
   along whole contigs, and that synteny was what chained loci to 400 genes. A locus is a maximal run
   of links with support ≥ 1. It is split again wherever a link's support falls below 0.34 × the smaller
   of the maxima on either side, so two abutting modules joined by only a few partners become two loci.
   With 2-gene bridging, the dip split cuts the longest locus from 131 genes to 51 and the 99th percentile from 41 to 23 genes.
3. **Families** (`03_families.py`): two loci are compared on **matched genes**. A gene counts as
   matched if the other locus has a gene in the same MMseqs2 family or a Spacedust hit partner.
   Spacedust's search is far more sensitive than 30% families: on gene families alone (Jaccard), 11k of
   19.8k loci ended up as singletons. Similarity is (mA + mB) / (nA + nB). An edge needs similarity
   ≥ 0.4 and ≥ 2 matched genes on each side, so one shared ABC ATPase cannot join two loci. Same-bin
   edges are dropped, and a cross-bin edge must be within 90% of the best match of one of its loci in
   the other bin. Leiden (CPM, weights = similarity) runs at resolution 0.1 for **families**, then
   again at 0.4 inside each family for **subfamilies**.
4. `04_sanity.py` checks known systems by Bakta product. `05_eval_operons.py` scores locus boundaries
   on curated E. coli transcription units.

## Outputs (`/data/scratch/metapanclust_data/catalogue/final/`)

- `chesterfield_catalogue.tsv`: one row per family. Columns: family, n_loci, n_bins, n_genera,
  median_genes, n_subfamilies, loci_per_bin, core gene families (in ≥ 50% of loci) with their commonest
  Bakta products, medoid locus with its **typical order** (gene-family IDs, `(-)` = opposite strand,
  oriented so most genes read forward), medoid products, member_loci.
- `chesterfield_members.tsv`: one row per locus (bin, contig, gene-index range, family, subfamily,
  genus, ordered gene families and products, contig-end flag).
- `chesterfield_loci.tsv`, `chesterfield_links.tsv`: loci and per-link support (genera, bins).
- `../fam30/genes.tsv`: every binned gene with its Spacedust ID, Bakta locus tag, product, gene family.

## Size distribution

| loci per family | 1 | 2 | 3-4 | 5-9 | 10-19 | 20-49 | 50+ |
|---|---|---|---|---|---|---|---|
| families | 4,121 | 1,095 | 897 | 696 | 261 | 94 | 14 |

| bins per family | 1 | 2 | 3-4 | 5-9 | 10-19 | 20-49 | 50+ |
|---|---|---|---|---|---|---|---|
| families | 4,121 | 1,232 | 895 | 640 | 220 | 65 | 5 |

Each family's median locus length is 2.5 genes (median over families; max 30.5). Single-locus families
(4,121) are loci whose partners matched below the similarity cut. 26% of them touch a contig end.

## Top families (medoid products)

| family | loci | bins | genera | contents |
|---|---|---|---|---|
| CF00000 | 98 | 85 | 57 | S10-spc-alpha superoperon: uS10, uL3, uL4, uL23, uL2, uS19, uL22, uS3, uL16, uL29, uS17, uL14, uL24, uL5, uS14, uS8, uL6, uL18, uS5, uL30, uL15, SecY, uS13, uS11, uS4, RpoA, bL17 (median 27 genes) |
| CF00001 | 86 | 46 | 28 | branched-chain amino acid ABC transporters (LivHMGF; 14 subfamilies) |
| CF00002 | 66 | 36 | 29 | peptide/nickel ABC transporters, family 5 (Dpp/Opp/Gsi; 17 subfamilies) |
| CF00003 | 64 | 54 | 43 | trigger factor, ClpP, ClpX |
| CF00004 | 62 | 51 | 40 | uS2, EF-Ts, PyrH, Frr |
| CF00005 | 56 | 41 | 30 | nuoA-N (median 11 genes) |
| CF00009 | 54 | 48 | 37 | SecE, NusG, uL11, uL1, uL10, bL12, RpoB, RpoC |
| CF00010 | 53 | 40 | 36 | 2-oxoglutarate dehydrogenase E1/E2/E3 |
| CF00011 | 53 | 52 | 41 | uL13, uS9 |
| CF00012 | 51 | 43 | 34 | IlvBN, IlvC |
| CF00013 | 51 | 41 | 28 | AmtB, P-II |
| CF00014 | 49 | 42 | 32 | dcw cluster: MraZ, RsmH, FtsI, MraY, MurD, FtsW, MurG, MurC… (median 13 genes) |
| CF00016 | 47 | 39 | 29 | GatCAB |
| CF00017 | 46 | 38 | 30 | GroES, GroEL |
| CF00018 | 46 | 41 | 31 | SucCD |
| CF00021 | 44 | 38 | 33 | fab cluster: FabD, FabG, ACP, FabF |
| CF00022 | 43 | 39 | 31 | his operon: HisB, HisH, HisA, HisF, HisIE |
| CF00024 | 41 | 39 | 32 | DnaA, DnaN, RecF, GyrB |
| CF00030 | 38 | 37 | 24 | alpha operon: IF-1, uS13, uS11, uS4, RpoA, bL17 (where it is not joined to CF00000) |
| CF00031 | 37 | 32 | 27 | SufBCD, cysteine desulfurase, SufU |
| CF00033 | 36 | 32 | 29 | F0F1 ATP synthase, atpB-E-F-H-A-G-D-C |
| CF00034 | 36 | 32 | 27 | DnaK, DnaJ |

## Sanity checks (`final/chesterfield_sanity.tsv`)

| system | genes (in loci) | bins | bins with the system in one locus | families holding 80% of its genes | top family: bins, purity |
|---|---|---|---|---|---|
| ATP synthase F0F1 | 348 (97%) | 54 | 27 | 3 | CF00033: 32, 0.60 |
| nuo | 600 (98%) | 48 | 16 | 3 | CF00005: 41, 0.55 |
| ribosomal proteins | 3,249 (98%) | 91 | 1 | 9 | CF00000: 85, 0.69 |
| flagellum | 231 (79%) | 13 | 1 | 15 | CF00282: 7, 0.45 |
| ABC transporters (any) | 1,653 (85%) | 78 | 1 | 92 | CF00001: 46, 0.38 |
| TRAP transporters | 118 (91%) | 13 | 0 | 8 | CF00065: 15, 0.30 |
| histidine biosynthesis | 388 (91%) | 42 | 5 | 13 | CF00022: 39, 0.55 |

"Bins" here means bins with at least 3 of the system's genes in loci. Purity is the fraction of a
family's genes whose product matches the system; the remainder is mostly "hypothetical protein" from
the Bakta light database, plus conserved neighbours.

- **Known operons come out as single families.** ATP synthase: one main family, CF00033. A second,
  CF00270 (12 loci, 11 bins), is the Na⁺-translocating variant (Bakta: "ATP synthase subunit beta,
  sodium ion specific"), a real biological split. nuo: CF00005 holds 41 bins. Ribosomal proteins fill
  9 families, matching the known operon map (S10-spc-alpha, L11-L10-rpoBC, uL13-uS9, bL21-bL27,
  bS6-bS18-bL9, IF3-bL35-bL20, uS2-tsf-pyrH-frr, …).
- **When an operon is in two loci in one bin, it is almost always a contig break.** That was 24 of 28
  such bins for atp and 27 of 34 for nuo. The dip/gap split caused only 4 (atp) and 5 (nuo).
- **ABC transporter hubs are split.** 1,653 ABC genes fall into 291 families; it takes 92 of them to
  hold 80% of the genes. The two hub families group transporter *classes* (Liv/HAAT; peptide/nickel
  family 5) at about 2 loci per bin. Their subfamilies split paralogues back to about one locus per bin
  (e.g. CF00001.0: 30 loci / 27 bins). The best-match filter removed only 2.5% of edges; the
  ≥ 2-matched-genes rule and CPM do the splitting.
- **Flagella fall into sub-operons** (fliEFGH-ATPase, fliM-PQR, flgFG-L/P-ring, flgK/fliS), each its
  own family. That is expected, since the flagellar regulon is several transcription units; only 13
  bins are motile.
- **Resolution trade-off.** At CPM 0.4 alone, the universal operons split by lineage (ATP synthase
  needs 7 families for 80%, ribosomal proteins 20). At 0.05, ABC hubs reach 94 loci / 45 bins.
  Families at 0.1 with subfamilies at 0.4 keep both levels in one table.

## E. coli transcription units (six-genome Spacedust benchmark)

Adjacent same-strand pairs (2,820; 1,727 same-TU). "Same locus" = both genes in one locus from
`02_loci.py` (support ≥ 1 other genome).

| bridge (max genes skipped) | same locus: acc / prec / recall | gap ≤ 50 OR (same locus AND gap ≤ 300): acc / prec / recall | far same-TU pairs recovered |
|---|---|---|---|
| 0 | 0.612 / 0.87 / 0.38 | 0.853 / 0.87 / 0.87 | 116 / 331 |
| 1 (used) | 0.631 / 0.85 / 0.43 | 0.847 / 0.85 / 0.88 | 134 / 331 |
| 2 | 0.630 / 0.82 / 0.45 | 0.842 / 0.84 / 0.89 | 144 / 331 |

Loci are somewhat larger than transcription units. The precision lost (0.82-0.87, against 0.90 for
plain Spacedust co-membership in `first_tests.md`) comes from conserved clusters that span neighbouring
TUs (superoperons, dcw), which is the intended unit here. The combined rule (0.85) is a little below the
0.874 in `first_tests.md`. That rule counts any shared cluster; here a pair counts only if both genes
are in one segmented locus.

## Caveats and next steps

- Loci need support from **another genus**. A cluster conserved only within one genus (e.g. inside the
  8 *Planktophila* bins) is not a locus. That cost is visible: 54% of genes are in loci, against 71% in
  any Spacedust cluster. To include genus-level accessory cassettes, count close relatives as one
  extra vote (support = other genera + 1 if any same-genus partner) and require ≥ 2.
- Gene families at 30%/80% are too fine across phyla (CF00000 has only 17 core families for 27 genes).
  Core composition and typical order are therefore reported per medoid. Spacedust hits, or Foldseek or
  eggNOG OGs, would give coarser families for the composition columns.
- 26% of singleton families touch a contig end. Partial loci should be matched by containment
  (overlap coefficient) when the shorter locus is truncated by a contig end.
- Scaling to GTDB species representatives (stage 1 proper): locus building is linear in Spacedust hits.
  Pair generation is the sparse product plus hits (210k candidate pairs here) and Leiden is fast. The
  costly steps are Spacedust itself and the MMseqs2 families, both cluster jobs at that scale.
