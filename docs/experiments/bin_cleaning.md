# Surprise-guided bin cleaning of the 93 Chesterfield bins (2026-11-18)

## Result

- **8.77 of 155.84 Mb (5.6%) removed from 45 of 93 bins**, 1,335 of 9,213 contigs. Almost all of it comes from
  the two coverage rules: off-population, 1,025 contigs and 6.56 Mb; no reads, 300 contigs and 1.51 Mb. The
  content rule removed 10 contigs (0.70 Mb). It flagged 11 more (0.24 Mb), which I kept because their coverage and
  composition fit the core.
- **CheckM2 (original → cleaned, 45 changed bins):** summed contamination falls by **75 points (380 → 305)** and
  summed completeness by 126 (2,465 → 2,339). No bin gains more than 0.5 points of contamination. Bins at
  MQ or better (≥ 50% complete, < 10% contamination) go from 14 to 17; HQ bins (≥ 90 / < 5) stay at 3. Of the 45
  bins, 13 lose more than 1 point of contamination and 6 lose more than 5 points of completeness. The losses are
  concentrated:
  - **Clear gains:**
    - TMED14 maxbin_102: contamination 17.6 → 0.7;
    - Pseudothioglobus vamb_711_sub: 21.6 → 8.7 at −0.6 completeness, from the content rule;
    - Flavobacterium lorbin_015: 10.6 → 3.4 at −0.8;
    - Nitrosopumilus maxbin_115: 8.3 → 2.3 at −6.3;
    - UBA9320 maxbin_050: 4.5 → 3.2;
    - UBA9145 maxbin_040: 4.2 → 2.8 at −2.7.
  - **Clear losses:** two clean bins whose core was seen in a single sample, where the off-population rule removed
    the bin's own genome:
    - lorbin_077: completeness 25.1 → 8.4, contamination was 1.2;
    - lorbin_059: 31.7 → 7.5, contamination was 0.1.
    
    Reverting those two leaves −85 completeness against −74 contamination.
  - **Mixed:**
    - CAIWMF01 maxbin_069_sub: 70.9/21.9 → 54.5/10.6;
    - vamb_163_sub: 41.7/8.0 → 30.8/2.8.
- **Surprising modules:** in the 18 main-assembly bins with a species pangenome, rare-present modules fall from
  **176 to 97**. Every drop is in three bins: TMED14 68 → 11, UBA9145 39 → 18, Pseudothioglobus sp913057955 3 → 2.
  The other 15 bins keep all their rare-present modules. These modules sit on contigs that coverage and content
  both accept.
- **UBA9145 (maxbin_040):** the 21 contaminant modules from `coverage_test` disappear, exactly those (17
  off-population, 4 no-reads). 15 of the 18 that are left are on contigs called core- or accessory-like, including the
  three flagellar modules (6, 20, 29). The other 3 (modules 7, 17, 18) are not in the `coverage_test` table. Module surprise drops from 245 to 125 bits. The bin's aligned fraction to its
  species (skani, best public genome) rises from 0.65 to 0.74. This test is partly circular: the same coverage
  classes found the 21 modules in the first place.
- **TMED14 (maxbin_102, 17.6% CheckM2 contamination):** 1.15 Mb (48%) removed, all but 7 kb by the off-population rule.
  Rare-present modules drop 68 → 11 and module surprise 263 → 50 bits. Every module with housekeeping genes
  (ribosomal proteins, gyrase, ATP synthase, IF-2) is gone. The 11 left are mostly hypothetical proteins (two carry an enzyme), 7 of them on
  contigs that carry persistent genes. Two checks do not use coverage at all:
  - skani aligned fraction of the bin to the TMED14 public genomes rises from **0.44 to 0.84**, while the
    reference fraction stays at 0.74, so the species' own content is not lost;
  - the removed 1.15 Mb matches no TMED14 genome and no other Chesterfield bin above skani's ANI floor.
- CheckM2: contamination 17.6 → **0.7**, completeness 84.5 → 61.4. The removed part is itself a coherent
    partial genome: **43.8% complete, 1.1% contaminated**. I read the completeness drop mostly as the original
    figure being inflated by a second genome's markers. Persistent-family completeness (0.70 → 0.68) and the
    reference-side aligned fraction (0.74 → 0.74) barely move. I cannot prove this without ground truth.
- **What cannot be validated:** there is no ground truth for which contigs belong to which population. CheckM2 can
  tell whether removal lowered contamination at the expense of completeness, but it cannot score the removal of a
  contig that carries no marker genes. That is most of the 4.4 Mb of removed short contigs. Bins whose core is
  present in only one sample (n_inf = 1, e.g. lorbin_077, which loses 67%) are where the off-population rule is
  most fragile. See "Risks".

## Rules (per contig, first match decides)

Inputs:

- coverage classes of every binned contig from `coverage_test` (`allbins_contigs.tsv`; 36 samples, contig-level
  depths, core = median profile of the bin's contigs of 10 kb or more);
- per-contig gene-family content for the 18 main-assembly bins whose GTDB-Tk species has a pangenome in
  `pangenome_survey`, plus the two earlier species. For each gene I take its family and k, the number of public
  genomes with that family. Gene surprise is -log2((k+1)/(n+2)), the "present" term of `binning_sim`.

| Rule | Action | Test | Contigs | Mb |
|---|---|---|---|---|
| R1 off-population | remove, any length | > 50% of depth in samples where the bin's core is absent | 1,025 | 6.56 |
| R2 no-reads | remove, any length | no mapped reads in any of the 36 samples | 300 | 1.51 |
| R3 content + coverage | remove | ≥ 10 kb, ≤ 10% of genes in a family that any public genome of the species has, gene bits/Mb above the bin's, **and** coverage does not vouch for it (not core-like, or composition outlier tnf > 1.25) | 10 | 0.70 |
| K1 content only | keep, flagged | the R3 content test passes, but coverage is core-like and composition is core-like | 11 | 0.24 |
| K2 short | keep | < 10 kb and not R1/R2 (ratio_adj, other-bin and surprise are noise there) | 4,349 | 24.39 |
| K3 keep | keep | everything else | 3,518 | 122.44 |

Removed bp by length: < 5 kb 2.73 Mb, 5-10 kb 3.05 Mb, 10-20 kb 1.48 Mb, ≥ 20 kb 1.52 Mb.

**Why R1 and R2 also apply below 10 kb.** The earlier reports warn that *ratios* and *surprise* are noise on short
contigs. The off-population and no-reads tests are not ratios, so I calibrated them separately. I took the known
members: contigs carrying ≥ 2 persistent genes in the 16 pangenome bins outside TMED14 and UBA9145. Of these 763
contigs, **1 is off-population and 0 have no reads**. That includes 286 contigs under 10 kb, and the single
off-population one is 5-10 kb. Across all 18 pangenome bins, 36% of contigs under 5 kb are off-population. The
other-bin class does not pass this test: 5 of 97 short known members are called other-bin, against 4.7% of all
short contigs. So other-bin contigs (all < 10 kb) are kept.

**Why R3 needs coverage to agree.** Among contigs not already removed by R1/R2, the content test flags 21. The 11 it flags with core-like
coverage and composition are mostly in Oceaniferula (vamb_189_sub, 8 contigs). That species has only 5 public
genomes, and its pangenome has no shell partition (3,033 persistent / 0 shell), so "no public family" is weak
evidence there. The other three are single contigs in BACL14 (lorbin_058, an 18%-complete bin), Planktomarina and
UBA10364. The binning simulation's AUC of 1.0 for contigs ≥ 10 kb came from pangenomes of 48-63 genomes. Most
species here have 3-13.

R3 removals:

- **Pseudothioglobus sp913057955 (vamb_711_sub, 3 public genomes): 0.60 Mb.**
  - A 498 kb contig with 548 genes, only 2 of them in a public family. It has GC 0.38 against the bin's
    0.30-0.31 and sits at accessory-like depth (0.53).
  - Three contigs of 24-48 kb at 2.8-5.8× core depth, also with no public family. Removing these removes the
    *nuo* module (module_12), which the survey had already called contamination from a more abundant organism.
  - skani aligned fraction to the public genomes rises from 0.54 to 0.77.
- **UBA10364 (vamb_193):** three contigs of 12-20 kb at ratio_adj 0.34-0.44, r = 1.0, core-like composition.
  This is the one ambiguous R3 call. It could be a co-occurring relative or a genuine island carried by
  about 40% of cells. The bin is 11% complete, so little rides on it.
- **Amylibacter (lorbin_025):** one 23 kb contig with GC 0.32 against the bin's ~0.5 (tnf 1.69).
- **Oceaniferula (vamb_189_sub):** two contigs.

## CheckM2

CheckM2 1.1.0 (danaSeq dana-mag-quality env, uniref100.KO.1 database), one run over the original and cleaned
versions of the 45 changed bins, plus the removed parts of the 13 bins that lost ≥ 0.2 Mb.

Bins that lost ≥ 0.05 Mb, or any bin with a species pangenome (all 45 are in `checkm2_compare.tsv`):

| Bin | Species | Mb | Removed Mb | Completeness orig → clean | Contamination orig → clean | Removed part compl. / contam. |
|---|---|---|---|---|---|---|
| maxbin_102 | tmed14 | 2.42 | 1.16 | 84.48 → 61.41 | 17.58 → 0.72 | 43.76 / 1.08 |
| maxbin_115 |  | 2.29 | 1.11 | 79.07 → 72.79 | 8.25 → 2.3 | 20.73 / 3.55 |
| lorbin_077 |  | 1.59 | 1.06 | 25.07 → 8.43 | 1.16 → 0.05 | 19.09 / 1.19 |
| maxbin_069_sub |  | 2.59 | 0.93 | 70.85 → 54.54 | 21.9 → 10.63 | 26.08 / 3.11 |
| lorbin_015 |  | 4.65 | 0.89 | 73.09 → 72.33 | 10.55 → 3.36 | 11.52 / 0.62 |
| vamb_711_sub | pseudothioglobus2 | 2.03 | 0.6 | 88.14 → 87.55 | 21.62 → 8.74 | 19.61 / 0.0 |
| lorbin_059 |  | 0.71 | 0.42 | 31.72 → 7.51 | 0.1 → 0.06 | 19.27 / 0.09 |
| maxbin_043 |  | 4.63 | 0.31 | 93.76 → 94.22 | 32.2 → 28.86 | 8.03 / 0.3 |
| maxbin_050 | uba9320 | 2.74 | 0.31 | 97.21 → 95.51 | 4.49 → 3.17 | 8.35 / 0.01 |
| semibin_238_sub |  | 3.11 | 0.25 | 44.05 → 41.03 | 10.53 → 9.42 | 6.05 / 0.02 |
| maxbin_040 | uba9145 | 2.81 | 0.23 | 94.79 → 92.08 | 4.23 → 2.75 | 5.39 / 0.0 |
| lorbin_037 |  | 4.05 | 0.22 | 70.73 → 70.25 | 10.95 → 10.58 | 3.69 / 0.0 |
| vamb_163_sub |  | 0.85 | 0.21 | 41.72 → 30.76 | 8.01 → 2.81 | 9.34 / 0.88 |
| semibin_272_sub |  | 3.48 | 0.14 | 69.31 → 68.09 | 45.4 → 45.18 |  |
| semibin_223 |  | 3.52 | 0.1 | 73.64 → 72.54 | 24.91 → 24.54 |  |
| lorbin_073_sub |  | 0.92 | 0.09 | 28.79 → 26.32 | 4.54 → 3.28 |  |
| lorbin_100_sub |  | 0.31 | 0.09 | 12.95 → 9.08 | 0.1 → 0.01 |  |
| semibin_267_sub |  | 4.08 | 0.06 | 81.55 → 80.44 | 27.33 → 27.15 |  |
| semibin_098 |  | 1.25 | 0.05 | 52.03 → 50.53 | 2.65 → 2.62 |  |
| semibin_028 |  | 1.59 | 0.05 | 70.62 → 69.3 | 6.76 → 6.25 |  |
| semibin_035_sub |  | 1.15 | 0.05 | 52.84 → 52.82 | 6.06 → 5.93 |  |
| vamb_193 | uba10364 | 0.33 | 0.05 | 13.17 → 11.22 | 0.02 → 0.01 |  |
| vamb_189_sub | oceaniferula | 2.59 | 0.03 | 51.48 → 51.45 | 4.88 → 4.53 |  |
| lorbin_025 | amylibacter | 1.71 | 0.02 | 93.8 → 94.01 | 0.77 → 0.14 |  |
| semibin_016 | hel1-33-131 | 1.28 | 0.01 | 61.97 → 62.25 | 3.56 → 3.62 |  |

What the removed parts tell us:

- **Coherent, low-contamination partial genomes.** TMED14's (44% complete), maxbin_069_sub's (26%) and
  maxbin_115's (21%) removed parts fit binned-in second organisms.
- **Cut from the bin's own genome.** lorbin_077's removed part (19%) is cut from a bin that had no contamination
  to begin with.
- **Little marker content.** UBA9145, UBA9320 and maxbin_043 lose 5-8% "completeness" worth of markers and almost
  no contamination signal in the removed part. They are mostly short contigs whose origin CheckM2 cannot judge.

## Module surprise after cleaning

Pangenomes were not rebuilt. For each main-assembly bin, the genes on removed contigs were dropped from its
PPanGGOLiN gff. Module completion and persistent completeness were then recomputed from the remaining families,
with the survey's model unchanged: present at completion ≥ 0.5; p = (k+1)/(n+2); rare-present means k/n ≤ 0.1.
Public prevalence is unchanged, and the original bins reproduce the survey's `per_bin.tsv` exactly. Marine-assembly
bins were not cleaned because there are no coverage classes for that assembly, so "carried by both bins" no
longer holds for the cleaned main bins.

| Species | Bin | Rare-present orig → clean | Module bits orig → clean | Persistent completeness | skani AF query orig → clean |
|---|---|---|---|---|---|
| TMED14 | maxbin_102 | 68 → 11 | 263 → 50 | 0.70 → 0.68 | 0.44 → 0.84 |
| UBA9145 | maxbin_040 | 39 → 18 | 245 → 125 | 0.92 → 0.92 | 0.65 → 0.74 |
| Pseudothioglobus sp913057955 | vamb_711_sub | 3 → 2 | 54 → 52 | 0.97 → 0.97 | 0.54 → 0.77 |
| UBA9320 | maxbin_050 | 8 → 8 | 60 → 60 | 0.92 → 0.92 | 0.67 → 0.84 |
| UBA10364 | vamb_193 | 1 → 1 | 4.1 → 4.1 | 0.11 → 0.11 | 0.59 → 0.69 |
| 13 others | | 57 → 57 | unchanged | unchanged | within 1 point |

- **UBA9320** loses 0.31 Mb (65 off-population or no-reads contigs) but no module. Its aligned fraction rises from
  0.67 to 0.84. So the bin carried about 0.3 Mb of foreign sequence that module surprise could not see, because
  foreign genes rarely complete a module. That fits `binning_sim` result 4: module surprise is blind to
  contamination.
- **Persistent completeness** falls only in TMED14 (0.698 → 0.676). The removed contaminant contigs carried genes
  that PPanGGOLiN had clustered into TMED14 persistent families (5 per contig, at 3-12% of their genes), as the
  survey suspected.
- The removed parts of all four bins above (TMED14, UBA9145, UBA9320, vamb_711_sub) match no public genome of the
  species and no other Chesterfield bin above skani's ANI floor (~80%). The contaminants are unbinned organisms,
  so the cleaned contigs cannot simply be moved to another bin.

## Risks and what cannot be validated

- **No ground truth.** Validation here rests on consistency between independent signals: coverage across samples,
  gene-family content against public genomes, CheckM2 single-copy markers, and skani alignment to public genomes.
  Of these, only CheckM2 and skani were not used to make the decisions. Both can only see contigs that carry
  markers or that align to a public genome.
- **UBA9145 is partly circular.** The 21 contaminant modules were found by the same coverage classes (focal mode)
  that R1/R2 apply (all-bins mode). Their disappearance shows that the two modes agree, not that the contigs are
  contamination. The independent evidence comes from earlier work: housekeeping genes in unanchored modules, only
  1 of 237 genes in a public family, and no public-genome alignment.
- **Bins with a one-sample core.** The off-population test compares a contig's depth against the samples where the
  core is present (≥ 0.2×). Where that is a single sample, a member whose population sits just below 0.2× in many
  other samples can be called off-population. lorbin_077 (n_inf 1, unclassified family) loses 1.06 of 1.59 Mb
  this way: 0.74 Mb off-population and 0.32 Mb no reads. Its removed contigs have the same GC as the kept ones.
- CheckM2 confirms that this went wrong: completeness 25.1 → 8.4 with almost no contamination to start with
    (1.2). The same failure hit lorbin_059 (31.7 → 7.5, contamination 0.1).
  - A one-sample-core guard (skip R1 when n_inf = 1) would not have been clean. It would also have blocked
    lorbin_015 (n_inf 1), a good cleaning: contamination 10.6 → 3.4 at −0.8 completeness.
  - A workable guard: do not clean a bin whose removal exceeds ~30% of its length unless an independent
    contamination signal agrees, i.e. CheckM2 contamination ≥ 5% or a rare-present module with housekeeping genes.
    That keeps TMED14, maxbin_115 and maxbin_069_sub and skips lorbin_077 and lorbin_059. It was chosen after
    seeing these results and has not been tested.
- **Short contigs.** 5.8 Mb of the 8.8 Mb removed are contigs under 10 kb. Calibration on 286 short known members
  supports R1/R2 there, but those members come from 16 mostly good bins. Short contigs without markers cannot be
  checked.
- **Content rule coverage.** Only 18 of 93 bins have a species pangenome, so R3 could only act there. The other 75
  bins were cleaned on coverage alone. Species with 3-6 public genomes make "no public family" weak evidence,
  which is why R3 also requires coverage to disagree.
- **Contig-level only.** A chimeric contig is kept or removed whole.

## Files

- Code (`experiments/bin_cleaning/`):
  - `content.py`: per-contig family content;
  - `clean.py`: decisions, cleaned FASTAs and calibration;
  - `rescore.py`: module surprise before and after cleaning;
  - `af_check.sh`: skani aligned fractions;
  - `run_checkm2.sh`;
  - `summarize.py`: CheckM2 comparison.

  Everything runs in under 15 minutes at ≤ 3 threads; CheckM2 takes most of it.
- Data (`/data/scratch/metapanclust_data/bin_cleaning/`):
  - `decisions.tsv`: every contig, with action, rule and evidence;
  - `per_bin.tsv`: Mb removed per rule;
  - `cleaned/*.fna`: all 93 bins, unchanged ones copied as they are;
  - `removed/`: removed parts of the 13 bins that lost ≥ 0.2 Mb;
  - `checkm2/quality_report.tsv` and `checkm2_compare.tsv`;
  - `module_surprise_bins.tsv` and `module_surprise_modules.tsv`;
  - `af/af.tsv`;
  - `clean.log`.

## For grex

Nothing needed.
