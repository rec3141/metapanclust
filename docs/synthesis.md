# Synthesis of the 2026-11-18 experiments

Sources: the reports in [experiments/](experiments/), plus [first_tests.md](first_tests.md),
[background.md](background.md) and two grex results not yet in a report. The first is Spacedust structure
mode on the six benchmark genomes, cited below as "grex, coordinator". The second is the 200-species
merge-rule run (job 8187204), cited as "grex 200-species"; its tables are in
`/data/scratch/metapanclust_data/binning_sim_genus/grex_results/results/summary_grex.md`. Still
running: the Chesterfield structure-mode catalogue run (173k proteins, grex GPU queue, results after
today).

## 1. State of the project

The pieces of the surprise model have each been tested once, on small data. Conserved-cluster detection
works: a Leiden catalogue of Chesterfield bins recovers textbook operons as single families. Conserved
neighbourhoods add a little to operon calling, up to a ceiling of about 0.88 accuracy, but only from
distant genomes. A closed-form detection model for contig breaks matches simulation to 0.5 points.
Within-species level-1 prevalence and a level-2 context model both run. Context beats prevalence by
about 30% of bits out of sample in two species. A surprise-per-Mb merge rule separates foreign from
same-species content. Across 200 species and 1,002 congeneric pairs it does so up to about 94-95% ANI,
provided gene families are cut at 95% identity and the threshold is calibrated per species. The main lesson of the day, though, is that **most large
"surprise" in real Chesterfield bins was contamination, not biology**. Cheap checks catch it: per-sample
coverage against the bin core, anchoring on contigs with persistent genes, and GTDB-Tk alignment
fraction. The model has to put those checks before any biological reading. A coverage-first cleaning pass
removes 5.6% of binned DNA and 57 of 177 rare-present modules; the bins that were not contaminated keep
all of their rare modules. Level-1 tree shrinkage has had one
internal test, leave-one-bin-out over the 93 Chesterfield bins, where it cuts held-out log loss by 16%
and the detection model is worth about as much again. That test is circular, because the catalogue was
built from the same bins. A small external test (2 genera) found that genus shrinkage does not beat a
species-level empirical-Bayes prior once a species has 5-20 genomes. Tree borrowing may therefore matter
mainly for data-poor species, which describes most of ours. Much of the evidence is weak: one or two query genomes per question,
reference sets made only of MAGs and SAGs, two overlapping co-assemblies counted as two bins, and
nanopore depth of about 4x in a single sample. Still untested: level 1 and level 2 at scale on an external reference (the
pangbank_scale grex run, submitted 2026-11-18 19:10: download then a 2-node job; results pending), the λ/ε parts of the detection model, population frequencies, and binning on a gold standard.

## 2. What we learned

### Operons and cluster detection

1. **Close genomes carry no operon evidence. Distant genomes do, and recall grows with their number.**
   For E. coli, conservation in genomes within Enterobacteriaceae has precision 0.62-0.67, which is the
   base rate of 0.61. Adding them drops the combined rule from 0.850 to 0.69-0.74. The best panel, 37
   genomes from Beta and beyond, reaches conserved-pair recall 0.487 at precision 0.897 and combined
   accuracy 0.878. Recall rises roughly with log(n) (1 → 0.08, 30 → 0.40) and has not saturated.
   [phylo_scaling.md](experiments/phylo_scaling.md)
2. **The pattern holds for B. subtilis, but the gain is smaller.** Close genomes again sit at base rate
   (0.74). Breadth gives the same log-like recall growth. Conservation adds only +0.003 to +0.006
   accuracy over a re-tuned gap rule (0.871 → 0.874-0.877 CV), against +0.03 for E. coli.
   B. subtilis spacers are longer: gap ≤ 100 bp is optimal, against 50 bp for E. coli, so gap
   thresholds must be fitted per clade. [phylo_scaling_bsub.md](experiments/phylo_scaling_bsub.md)
3. **Neighbourhood evidence caps operon accuracy at about 0.88.** The remaining errors are adjacent
   different-TU pairs (91 false calls at gap ≤ 50) and lineage-specific operons that no distant genome
   conserves (170 misses). More genomes will not fix these; promoters, terminators or expression would.
   [phylo_scaling.md](experiments/phylo_scaling.md)
4. **Structure mode (ProstT5 3Di) finds many more conserved genes but does not improve operon calls.**
   On the full six-genome set (grex, coordinator) it found 3,762 clusters against 2,144 with MMseqs2.
   Genes covered rose from 18.1% to 26.5%, with the largest gains in distant genomes: mja 4 → 15%,
   syn 6 → 13%, mtb 14 → 24%. Only 170 genes were found by MMseqs2 alone, against 2,173 by 3Di alone.
   E. coli conserved-pair recall rose from 0.378 to 0.453, while precision went from 0.897 to 0.890.
   Combined-rule accuracy was 0.874 in both modes, because the gained pairs are mostly close ones the
   gap rule already calls. The local subsets showed the same pattern: +9-25% genes, combined accuracy
   unchanged, and the extra hits are remote homologs at 15-35% identity.
   [structure_mode.md](experiments/structure_mode.md)
5. **Structure mode is affordable on GPU only.** ProstT5 on CPU takes about 18 core-s per protein. On
   one L40S it runs at 31 proteins/s (23,901 proteins in 767 s), and the Foldseek search took 2 s. That
   is about 1.5 GPU-h for Chesterfield's 173k binned proteins. The workstation GV100 (sm_70) cannot run
   the precompiled build. A rough comparison: 3Di on six distant genomes gives about the recall (0.45)
   that MMseqs2 gets from 37 distant genomes (0.49). Gene calls differ between those runs, so treat this
   as indicative only.
6. **Contig breaks lose clusters in a predictable way, and real assemblies break inside conserved
   clusters less often than random.** Random breaks at Chesterfield contig lengths leave 83% of curated
   operons and 81% of Spacedust clusters intact. That falls to 74/70% at bin N50 10-20 kb, and to 26%
   for clusters spanning 10-20 kb. The closed form P_intact = c·(1 − β(1 − f(s))) matches simulation
   (MAE 0.005, r 0.9995). In real bins, conserved clusters touch contig ends only 25% as often (in
   excess over the null) as random breaks predict (β ≈ 0.25), because breaks fall in variable,
   hypothetical-rich DNA. Without this term, absence surprise in fragmented bins is overstated (example:
   2.0 bits instead of 0.57). [fragmentation.md](experiments/fragmentation.md)

### Catalogue

7. **Leiden on a matched-gene locus graph removes the giant family.** The catalogue has 20,849 loci in
   7,178 families, 3,057 of them multi-locus. The largest family is the S10-spc-alpha superoperon (98 loci
   in 85 bins), where the earlier attempt had a single family of 14,709 loci. Connected components on
   the same graph still give one component of 12,053 loci, so the fix comes from Leiden. Two Leiden
   seeds agree at ARI 0.975. [catalogue.md](experiments/catalogue.md)
8. **Known systems come out right.** ATP synthase (with the Na⁺ variant as a separate family), nuo, the
   9 ribosomal-protein operons, his, dcw, fab and suf each form single families. A split in one bin is
   almost always a contig break (24/28 for atp, 27/34 for nuo). ABC transporters are split into
   subfamilies by class, and flagella into sub-operons. Locus boundaries are somewhat larger than TUs:
   same-locus precision on E. coli TUs is 0.85, against 0.90 for plain Spacedust co-membership.
   [catalogue.md](experiments/catalogue.md)
9. **The catalogue covers 54% of binned genes.** That compares with 71% in any Spacedust cluster and
   58% in the Spacedust paper. The gap comes from requiring support from another genus: clusters
   conserved only within a genus are lost. [catalogue.md](experiments/catalogue.md)

### Population and coverage

10. **Contig-level coverage ratios measure accessory content only on contigs of 10 kb or more.** Below
    10 kb, known members look "accessory" as often as other contigs (46% at 3-5 kb). At 10 kb and above,
    394 contigs (7.5 Mb, 72 bins) are accessory-like, with ratio_adj median 0.44 (IQR 0.32-0.53). That is
    a candidate list, not a measurement. [coverage_test.md](experiments/coverage_test.md)
11. **Per-sample coverage is the strongest contamination signal.** Over all binned DNA, 11.1% of
    contigs (4.2% of bp, 1,025 contigs, 6.6 Mb) are off-population: their reads come from samples where
    the bin's core is absent. Correlation across samples barely separates populations here: 81 of 93 bins
    correlate at r ≥ 0.8 with another bin. Composition (TNF) rejects 92-98% of co-occurring contaminant
    contigs that coverage alone passes. [coverage_test.md](experiments/coverage_test.md)
12. **Genuine rare modules are carried by most cells.** In 8 clean bins, module/core depth ratios run
    0.62-1.17 (10th-90th percentile, 48 modules). Only 4 Pontimonas modules (0.03-0.32×) are candidates
    for within-population variable content. [pangenome_survey.md](experiments/pangenome_survey.md)

### Contamination and binning

13. **Large surprise in real bins is mostly contamination.** Of UBA9145's 36 Chesterfield-only modules,
    21 are contamination: 17 are off-population, almost all from barcode16, and 4 have no reads. Only 1
    of 237 genes on those contigs falls in a public UBA9145 family, against 216 of 290 on contigs that
    follow the core. [coverage_test.md](experiments/coverage_test.md) TMED14's 68 rare modules come from
    a second Flavobacteriales genome: they carry ribosomal proteins, gyrase and ATP synthase, at GTDB AF
    0.68-0.71 and CheckM2 contamination 17.6%. In the other 7 species, 47 of 50 rare modules are anchored
    and none has decoupled coverage. [pangenome_survey.md](experiments/pangenome_survey.md)
13a. **Coverage-first cleaning removes most spurious surprise.** The cleaning pass removed 8.77 of
    155.84 Mb (5.6%) from 45 of 93 bins, 1,335 contigs. Almost all of it came from the two coverage
    rules: off-population (6.56 Mb) and no reads (1.51 Mb). A content rule, which requires coverage to
    disagree as well, removed 10 contigs (0.70 Mb). The two coverage rules hold on contigs under 10 kb:
    of 763 known-member contigs, 1 was off-population and 0 had no reads.

    Rare-present modules in the 18 main-assembly bins with a pangenome fell from 177 to 120, all in three
    bins:
    - TMED14 68 → 11 (263 → 50 bits), with every housekeeping module gone. skani AF to public TMED14 rose
      from 0.44 to 0.84 while the reference fraction stayed at 0.74.
    - UBA9145 39 → 18 (245 → 125 bits); exactly the 21 coverage_test modules went.
    - Pseudothioglobus vamb_711_sub 3 → 2, losing the contaminant nuo module and a 498 kb contig at
      GC 0.38 against the bin's 0.30. AF rose from 0.54 to 0.77.

    UBA9320 lost 0.31 Mb and no module, yet its AF rose from 0.67 to 0.84: contamination that module
    surprise cannot see. The other 15 bins kept all 58 of their rare modules, which supports reading them
    as population content. The removed DNA matches no public genome of the species and no other
    Chesterfield bin above about 80% ANI. It comes from unbinned organisms, so it cannot be moved to
    another bin.

    **CheckM2 before and after (45 changed bins).** Summed contamination fell by 75 points (380 → 305)
    and summed completeness by 126 (2,465 → 2,339). No bin gained more than 0.5 contamination points. MQ
    bins went from 14 to 17 and HQ bins stayed at 3. Clear gains: TMED14 contamination 17.6 → 0.7 (its
    removed part is a coherent 44%-complete second genome), vamb_711_sub 21.6 → 8.7, lorbin_015
    10.6 → 3.4, maxbin_115 8.3 → 2.3, UBA9145 4.2 → 2.8. Two bins whose core sits in a single sample
    were damaged: lorbin_077 (completeness 25 → 8) and lorbin_059 (32 → 8), with almost no contamination
    to begin with. The most contaminated bins barely change (semibin_272_sub 45%, maxbin_043 32%,
    semibin_267_sub 27%, semibin_223 25%). Coverage cannot separate co-abundant genomes; those bins need
    the content rule or re-binning.

    Weak points: the UBA9145 result is partly circular (the same coverage classes found the 21 modules),
    and no ground truth exists. The report proposes a guard: skip cleaning that removes more than 30% of
    a bin unless CheckM2 contamination is ≥ 5% or a housekeeping rare module agrees. It was chosen after
    seeing these results and is untested. [bin_cleaning.md](experiments/bin_cleaning.md)
14. **Cheap flags.** A bin to distrust has GTDB-Tk AF < 0.8 together with housekeeping genes in rare
    modules. A module to distrust sits on contigs without persistent genes. Both flags caught both bad
    bins before any coverage was looked at. A module far above core depth (Pseudothioglobus nuo, 6.7×)
    comes from a more abundant organism. [pangenome_survey.md](experiments/pangenome_survey.md)
15. **Surprise per Mb separates foreign from same-species merges, and marker duplication ranks them the
    wrong way round.** The merge rule is excess = ΔS/Mb_added − S(bin)/Mb. Across two species and 5-40%
    additions it gave AUC 1.0. Marker duplication gave AUC 0.0-0.4, because a strain merge duplicates
    more markers than a contaminant. [binning_sim.md](experiments/binning_sim.md)
16. **With 95%-identity families, the merge rule holds up to the species boundary.** The test covered
    200 focal species in 195 genera: 1,002 focal × contaminant pairs, 796 contaminant species, 0
    failures (grex 200-species). The first wave (8 species, 80% cut) put failure at 93-95% ANI. That
    came from the family cut, not from the rule (C12). Medians per ANI bin, family cuts 80 / 90 / 95%:

    | ANI | AUC (excess) | Calibrated rejection, whole bin | Calibrated rejection, 70% bin | Single 20-50 kb contig AUC |
    |---|---|---|---|---|
    | 88-90% | 0.97 / 1.00 / 1.00 | 0.97 / 1.0 / 1.0 | 0.94 / 1.0 / 1.0 | 0.81 / 0.93 / 0.96 |
    | 90-92% | 0.92 / 0.98 / 1.00 | 0.78 / 1.0 / 1.0 | 0.78 / 0.97 / 1.0 | 0.74 / 0.86 / 0.95 |
    | 92-94% | 0.89 / 0.95 / 0.99 | 0.72 / 0.89 / 1.0 | 0.58 / 0.86 / 1.0 | 0.68 / 0.78 / 0.89 |
    | 94-95% | 0.82 / 0.87 / 0.96 | 0.35 / 0.54 / 0.94 | 0.36 / 0.49 / 0.89 | 0.62 / 0.68 / 0.80 |
    | 95-96% | 0.77 / 0.82 / 0.90 | 0.31 / 0.42 / 0.78 | 0.25 / 0.28 / 0.69 | 0.60 / 0.63 / 0.72 |
    | ≥ 96% | 0.75 / 0.75 / 0.82 | 0.29 / 0.32 / 0.43 | 0.24 / 0.22 / 0.26 | 0.61 / 0.60 / 0.65 |

    Same-species merges are wrongly rejected 4-8% of the time throughout (median calibrated strain
    rejection). Below about 88% ANI, every cut gives AUC ≥ 0.98 and calibrated rejection 1.0. The 10th
    percentile of AUC is lower: at the 95% cut it is 0.91 at 92-94% and 0.83 at 94-95%, so some genera
    fail earlier. Shared-gene fraction predicts failure better than ANI does. At the 95% cut, AUC is
    ≥ 0.99 below a shared fraction of 0.7, 0.92 at 0.7-0.8 and 0.76 at 0.8-0.85. Marker duplication
    again ranks contaminants and strains the wrong way (AUC 0.28-0.52). Congeners at ≥ 95% ANI are
    GTDB species-boundary cases, where gene content cannot be expected to separate them.
    [binning_sim_genus.md](experiments/binning_sim_genus.md), grex 200-species
17. **The zero threshold is miscalibrated for incomplete bins; a per-species calibrated threshold fixes
    it.** A 70%-complete base costs 880-1,290 bits/Mb, so close congener contigs look cheap by comparison.
    At the 95% cut, the zero rule rejects only 32% of 94-95% ANI merges into 70% bins (0.81 at 92-94%).
    The calibrated threshold rejects 89% (1.0 at 92-94%). The calibrated threshold is the 95th percentile
    of strain merges, set on a calibration half of the held-out genomes. Its median cost in strain
    rejections is 4-8%, not the 18-21% seen for the single Bartonella test species. Content-only
    surprise is not the fix (it rejects 23-38% of strain merges).
    [binning_sim_genus.md](experiments/binning_sim_genus.md), grex 200-species
18. **Single contigs under 10 kb cannot be judged on content.** A bin's own short contigs are rejected
    26-73% of the time, because islands sit on short contigs. This agrees with finding 10: 10 kb is the
    floor for both coverage and content evidence. [binning_sim.md](experiments/binning_sim.md)
19. **A reference set can contain contaminated MAGs.** Public MAG-120531 genome GCA_963976715.1 has
    17.8% foreign genes on 51 small contigs, while PPanGGOLiN reports 0.2% contamination for it.
    [binning_sim.md](experiments/binning_sim.md)

### Surprise models

20. **Level 1 as prototyped is within-species Beta(1,1) prevalence, and its bits are not comparable
    across species.** A module never seen among n genomes costs log2(n+2) bits: 5.75 at n = 52, 2.3 at
    n = 3. Counts of rare-present modules travel better. In clean bins of at least 50% completeness the
    median is 37 bits per bin (28 present, 10 absent) and 3 rare-present modules. Absence contributes
    only 2-16 bits per bin. [pangenome_survey.md](experiments/pangenome_survey.md)
20a. **Tree shrinkage helps inside Chesterfield, and detection helps as much.** The test scored the
    catalogue's 7,178 families in 93 bins, each bin left out in turn, with a hierarchical beta-binomial
    along the GTDB-Tk lineage. Held-out log loss fell from 0.1864 to 0.1561 bits per cell (−16%), and
    by 25% on families seen in ≥ 2 bins. Cross-fitted κ gave the same number. Fitted κ is 3 from
    phylum to family and 30 at genus, which is strong shrinkage; a fully pooled prior is 6% worse. In a
    two-state comparison, detection is worth about as much as taxonomy (flat 0.180 → taxonomy 0.155,
    flat + detection 0.167, both 0.137). Calibration is good below P = 0.05 and somewhat
    underconfident above 0.2. This is not an external test: the families come from the same bins (see
    C10). [catalogue_surprise.md](experiments/catalogue_surprise.md)
20b. **Absence surprise flags contamination; presence surprise does not by itself.** Absence bits per
    expected Mb correlate with CheckM2 contamination (Spearman 0.66, partial 0.65 given completeness).
    In bins ≥ 70% complete, absence separates ≥ 10%-contaminated bins with AUC 0.85, against 0.77 for
    presence. The likely mechanism is that contaminant markers inflate CheckM2 completeness. Presence
    bits per Mb do not track contamination over all bins (ρ 0.12). Loci on off-population contigs do
    score higher (median 3.6 vs 1.6 bits, p = 1e-22). Accessory-like contigs score almost as high,
    though (38% > 5 bits vs 42%), so level-1 presence surprise means "not this clade's usual content"
    and needs coverage, anchoring or context to say more. Half of all presence bits (48%) come from
    singleton families. [catalogue_surprise.md](experiments/catalogue_surprise.md)
20c. **Absence of core operons costs at most about 3 bits, and the strongest absences are a mix of
    biology, catalogue splits and binning loss.** No absence exceeds 2.89 bits; detection sets that
    ceiling, since the median P_intact of a missing 11-17 kb operon is 0.07-0.21. Examples of each
    cause: UBA845 lacks nuo and carries Nqr instead (biology); Peredibacter's ATP synthase sits in a
    split 2-bin family (catalogue); Nanopelagicus lacks the S10-spc contig (composition loss λ_k). Per
    family, only CF00000 behaves as universal (p̂ 0.95); for atp and nuo, p̂ is 0.6-0.75. Systems split
    across families need system-level scoring. [catalogue_surprise.md](experiments/catalogue_surprise.md)
20d. **In a first small test, genus shrinkage does not beat a species-level empirical-Bayes prior once
    the species has 5-20 genomes.** This is a local test of 2 genera (Marinisoma, Flavobacterium) on
    PanGBank pangenomes. Against a single species-wide Beta(a, b) fitted on the m training genomes
    (sp_eb), genus shrinkage loses in every stratum and at m = 5, 10 and 20. In all clusters at m = 5 it
    costs 2,516 vs 2,156 bits per test genome, winning for 0/11 focal species. Even the test-set-optimal κ
    loses (2,483). Shrinkage looks 55-60% better against Beta(1,1), but that is an artefact of the
    cluster universe: Beta(1,1) pays for thousands of clusters absent from the species. Fitted κ is small
    (median 0.5), so species prevalences are U-shaped around the genus mean, and most clusters act as
    species markers. The only gain is a 1% improvement from a 50/50 mix on clusters the focal species
    shares with congeners. The 40-genus grex run that would decide this has not been submitted. This is
    2 genera, the MAG-only Marinisoma has no completeness model, and part of the loss is a 50%-identity
    mapping artefact. [pangbank_scale.md](experiments/pangbank_scale.md)
21. **Level 2 (context) beats prevalence out of sample by about 30% of bits.** In leave-one-genome-out
    tests, an L2 logistic on accessory families with C = 0.1 gave 0.395 vs 0.583 bits (UBA9145 modules)
    and 0.374 vs 0.529 (MAG-120531). Hyperparameters transfer between the two species. The gain is still
    rising at full size, about +0.03-0.04 bits per doubling of training genomes. Everything it learns is
    positive co-inheritance of island fragments; no mutually exclusive pair was found. A Jaccard-based
    phylogeny proxy adds only 0.04-0.06 bits. [context_model.md](experiments/context_model.md) In 4
    isolate-rich PanGBank species (20-318 genomes; indicative only), the logistic removes a median
    42-63% of the prevalence loss, against 30% in the two MAG species. A neighbour (clade) baseline takes
    11-32% of that. Much of the extra gain is therefore within-species clade structure, not context. The
    full 291-species grex design separates the two but has not been submitted.
    [pangbank_scale.md](experiments/pangbank_scale.md)
22. **Co-occurring novel modules must be scored as one block.** The 36 UBA9145 modules cost 207 bits
    per bin when scored independently. Given each other, each costs 0.29 bits, so the block is worth
    about 5.8 bits. Context sharpens the tail (LOO maximum 8.7 bits, against 4.7 for prevalence). A 10%
    prevalence mixture caps it at 6.4-7.5 bits. [context_model.md](experiments/context_model.md)

### Detection

23. **Detection model.** The fragmentation closed form is validated (finding 6). The P_lost
    approximation has MAE 0.023, with the error confined to low-completeness bins. model_spec adds a
    three-state observation (intact / partial / absent), composition loss λ_k, an assembly term a_j and
    false presence ε_jk. None of these has been estimated yet.
    [fragmentation.md](experiments/fragmentation.md), [model_spec.md](experiments/model_spec.md)

### Reference data (PanGBank)

24. **PanGBank adds no genomes for our species.** It covers 4 of 19 Chesterfield species in r226 and 2
    in r232; the 2.x filter (representative ≥ 85% complete) drops MAG-120531. Where it does cover a
    species, its genome sets are identical to ours (67/67, 52/52, 23/23, 17/17). Family-level surprise
    agrees within 1-20%. [pangbank.md](experiments/pangbank.md)
25. **PanGBank modules hide novel content.** They are defined on public genomes only, so UBA9145 falls
    from 245 to 42 module bits. The content moves into projection's "specific" genes: 57-68 islands and
    328-391 bits. Spurious surprise is hidden, not removed, unless specific islands are scored as k = 0
    units. [pangbank.md](experiments/pangbank.md)
26. **A sister species' pangenome is a genus prior, not a species prior.** At 80% identity, projection
    onto a sister species drops TMED14 completeness from 0.70 to 0.39. At 50% identity the core comes
    back, but 900-1,200 genes remain specific. Useful pieces: per-genome CheckM2 metadata (for
    correcting public-MAG prevalence), standardized RGPs and spots, and 10-30 s projection per species.
    [pangbank.md](experiments/pangbank.md)

## 3. Contradictions and corrections

| # | Claim | Conflicting result | Trust | Why |
|---|---|---|---|---|
| C1 | first_tests: the 36 UBA9145-only modules (flagella first) mean "either a motile local population or contamination" | coverage_test: 21/36 are contamination; flagellar module_6 follows the core at ratio_adj 0.50 with core composition, so it is "real, carried by about half the population". pangenome_survey: all three flagellar modules are unanchored at 0.38-0.62×, so contamination is "the better explanation" and a motile population is "unlikely". context_model: the block looks foreign (153/208 genes on contigs without persistent genes). | **Contamination is settled for the 21 off-population and no-reads modules. The flagellar modules are unresolved, leaning to population content at φ ≈ 0.5.** | coverage_test is the more specific test: per-sample off-population depth, length-corrected ratio, TNF, and 11/26 module_6 genes in public UBA9145 families. pangenome_survey's case is guilt by association with unanchored housekeeping modules, and coverage_test places those on barcode16 contigs. "Unanchored" means no *persistent* gene on the contig; shell/cloud genes from public UBA9145 are present on it. A co-occurring motile relative at about half abundance also fits, and all of it rests on about 4× depth in one sample. Retire the bare "motile local population" claim in first_tests.md. bin_cleaning keeps all three flagellar modules (6, 20, 29) on core- or accessory-like contigs and removes exactly the 21 contaminant modules, but it uses the same coverage classes, so it adds consistency rather than independent evidence. Decide with reads spanning module_6 into core contigs (roadmap R1). |
| C2 | fragmentation.md: absent score −log2(1 − p·P_intact), β = 1 | model_spec §2: three-state observation, P(A) = p·P(A given z=1) + (1 − p)(1 − ε), with λ, a_j, ε | **Both, in sequence.** Use fragmentation's form now; switch to the three-state model once λ and ε are estimated. | model_spec itself shows fragmentation's form is the conservative bound (fewer absence bits) and a special case of the three-state model. fragmentation's form is validated (MAE 0.005); λ, a_j and ε have no estimates. The two agree that a cluster truncated at a contig end is not an absence. Note also that the prototype used in pangenome_survey, binning_sim and pangbank, −log2(1 − p·c), has no fragmentation term and so overstates absence (example: 2.0 vs 0.57 bits). The effect is small in practice, since absence contributed 2-16 bits per bin. |
| C3 | phylo_scaling (E. coli): count conservation only from genomes outside the query's GTDB order | phylo_scaling_bsub: other GTDB orders of NCBI Bacillales still give precision 0.854, near the 0.846 of same-order families; ≥ 0.90 starts at other Bacilli (class level) | **phylo_scaling_bsub: set the floor by divergence, not rank.** | GTDB splits NCBI Bacillales into about ten orders, so the same rank means different divergence. The E. coli "order" result partly reflects GTDB's broad Enterobacterales; even there, precision reached 0.87 only from other Gamma orders, and the best panel started at Beta. Use an AAI or marker-distance floor calibrated to precision ≥ 0.9. Two queries are not enough to fix the number (R5). |
| C4 | catalogue.md: a locus needs support from ≥ 1 *other genus* | phylo_scaling and phylo_scaling_bsub: same-family and same-order genomes give base-rate precision, so conservation by descent spans TU boundaries | **Open.** Catalogue families are probably partly descent-conserved neighbourhoods rather than co-selected units. | The genus rule was chosen to stop whole-contig synteny chaining, not for precision. Its E. coli precision (0.85) comes from the six distant benchmark genomes, not from Chesterfield's genus-level partners. Re-run locus support counting only distant lineages and compare (R5). |
| C5 | related_work: "Do not rebuild species pangenomes; reuse PanGBank". model_spec O13 leans the same way. | pangbank.md: 4/19 species covered (2/19 in r232), identical genomes, no freshwater species | **pangbank.md.** | It measured coverage directly; related_work reasoned from collection size. Use PanGBank where it exists and keep building our own pangenomes for MAG-only marine and freshwater species. |
| C6 | binning_sim: incompleteness does not affect the merge rule "because it compares a bin with itself" | binning_sim_genus: the zero threshold accepts 83-90% ANI congeners into 70%-complete bins (60-64% rejected, against 87-97% for whole bins) | **binning_sim_genus.** | The first wave's contaminants cost 5,000-5,900 bits/Mb, so any threshold worked. The base term S(bin)/Mb grows with incompleteness, and that matters once the contaminant is close. model_spec §7.3 already warns against τ = 0. grex 200-species confirms this at scale: at the 95% cut, the zero rule rejects 0.32 of 94-95% ANI merges into 70% bins, against 0.89 with the calibrated threshold. |
| C7 | context_model: both UBA9145 bins carry "10-11% CheckM2 contamination" | pangenome_survey table: maxbin_040 at 95% complete, **4.2%** contamination | **pangenome_survey (4.2%).** | bin_cleaning's independent CheckM2 1.1.0 run gives maxbin_040 94.8% complete and 4.23% contaminated (2.75% after cleaning). The "10-11%" in context_model is probably from another bin or run and should be corrected there. It changes no conclusion. |
| C8 | first_tests: combined rule 0.874 | catalogue: 0.847-0.853 with segmented loci | **Both correct.** They measure different things. | Locus segmentation joins neighbouring TUs (superoperons) on purpose. The 0.874 counts any shared cluster. |
| C9 | model_spec §6.5 frames the flagellar question as φ ≈ 1 (biology) vs variable φ (contamination) | coverage_test: φ ≈ 0.5 with constant tracking | **Extend the spec.** | A constant φ well below 1 is a third outcome (accessory content in about half the cells, or a co-abundant relative). With one informative sample, H_const vs H_var cannot be tested at all. |
| C10 | model_spec §3.4 A: EM detection correction of training counts using the hierarchical prior, with leave-one-out by subtracting the scored genome's counts | catalogue_surprise: as written this leaks. It gives 0.110 bits per cell, 30% "better" than the honest 0.156: a held-out bin's presence returns through its clade-mates' posteriors (a TMED14-only family got p̂ 0.14 for TMED14). | **catalogue_surprise; fix the spec.** | This was measured directly with a clean contrast (the flat-prior training posteriors). The spec needs a per-bin refit or leave-two-out posteriors. Any future level-1 number computed with §3.4 A as written is suspect. |
| C11 | binning_sim: absolute surprise per Mb confounds incompleteness with contamination; contamination is a presence (foreign-gene) signal | catalogue_surprise: in real bins, contamination shows in *absence* bits per expected Mb (AUC 0.85), not in presence bits per Mb (ρ 0.12) | **Both, for different settings.** | binning_sim used a true completeness from markers in the species' own pangenome, with foreign genes scored as such. catalogue_surprise uses CheckM2 completeness, which contamination inflates, so the missing content shows as absence. That is useful as a flag, but it rests on a CheckM2 artefact; test it on known chimeras (TMED14, vamb_711_sub) and CAMI before relying on it. |
| C12 | binning_sim_genus: "the rule holds up to about 90% ANI and fails at 93-95%" (8 species, 80% family cut). Also: the honest calibrated threshold rejected 18-21% of strain merges (Bartonella). | grex 200-species: at a 95% cut, calibrated rejection is 1.0 up to 94% ANI and 0.94 at 94-95%, with AUC 0.96 at 94-95%. Strain rejection is 4-8%. | **grex 200-species.** | It covers 25× more species (195 genera, 1,002 pairs) with the same code, and it isolates the cause: the 80% family cut, not the surprise rule. The first wave's 5 pairs above 90% ANI came from 3 genera. The Bartonella strain-rejection figure was one species with 3 calibration genomes. Caveat for both: isolates only, 50-kb pseudo-contigs and 15 training genomes, so real MAG bins are still untested (R2). |
| C13 | catalogue_surprise: taxonomy shrinkage cuts held-out log loss by 16% vs a flat prior. model_spec §3 makes rank shrinkage the core of level 1. | pangbank_scale (local test, 2 genera): genus shrinkage loses to species-only empirical Bayes at m = 5-20 genomes, even with the oracle κ | **Not contradictory yet; both are weak.** Treat shrinkage as a fallback for species with little or no data, not as a default adjustment. | The baselines differ. catalogue_surprise's "flat" is one Beta for all bins, and most Chesterfield bins are the only member of their species or genus, so any clade information helps. pangbank_scale gives the species its own 5-20 genomes and a fitted species-wide Beta, which is the stronger baseline. Taken together: borrowing from the tree pays when the species has no data, and may not pay once it has about 5 genomes. catalogue_surprise is also circular (C10 aside, the families come from the same bins), and pangbank_scale covers 2 genera. The 40-genus run (R4) decides. If it agrees, change model_spec §3 to shrink the species-wide Beta(a, b) per cluster towards the genus, or to use genus information only for shared clusters. |

## 4. What does not work, or was ruled out

- **Single-linkage or connected components for catalogue families.** They give one giant family even
  after edge filtering (12,053 loci). [catalogue.md](experiments/catalogue.md)
- **Close genomes as operon evidence**, in any form: union, k-of-n votes, or per-level weights in a
  logistic model. Same-genus, same-family and same-order genomes sit at base-rate precision.
  [phylo_scaling.md](experiments/phylo_scaling.md), [phylo_scaling_bsub.md](experiments/phylo_scaling_bsub.md)
- **Pushing operon accuracy past about 0.88 with more genomes or structure search.** The cap is set by
  adjacent different-TU pairs and lineage-specific operons.
- **Fixed gap thresholds across clades** (E. coli 50/300 gives 0.77-0.81 on B. subtilis).
- **ProstT5 on CPU or on the workstation GV100.** About 7,400 core-h for the full Chesterfield assembly;
  the GV100 aborts with a CUDA error. [structure_mode.md](experiments/structure_mode.md)
- **Marker duplication**, strict or relaxed, **for telling contamination from a strain merge.** AUC
  0.0-0.6 in both waves. Adding a multicopy term to the merge score also makes it worse (AUC 0.74-1.0).
- **Absolute surprise per Mb as a bin-quality score.** It confounds incompleteness with contamination
  (AUC 0.25-0.34 at 5%). Content surprise per gene works better (0.80-0.91).
  [binning_sim.md](experiments/binning_sim.md)
- **Module-level surprise for contamination** (AUC 0.12-0.37): foreign genes never complete a module.
- **80%-identity gene families for merge scoring.** They push the merge rule's failure point down to
  about 92% ANI (AUC 0.89 at 92-94%). Use 95%. **The zero threshold** (excess > 0) for incomplete bins:
  use a per-species calibrated threshold. **Content-only merge scoring**, which rejects 23-38% of strain
  merges. **Gene content to split congeners at ≥ 95-96% ANI**, even at the 95% cut (calibrated
  rejection 0.43, AUC 0.82). That needs SNV or strain-resolved evidence.
- **Ratio, content or surprise judgements on contigs under 10 kb.** The presence/absence coverage tests
  (off-population, no reads) do work there.
- **The off-population rule alone on bins whose core is present in one sample.** It cut lorbin_077 and
  lorbin_059 to 8% completeness. **Coverage cleaning for co-abundant chimeras**: the 25-45% contaminated
  bins barely change. [bin_cleaning.md](experiments/bin_cleaning.md)
- **Coverage correlation alone to separate Chesterfield populations** (81/93 bins at r ≥ 0.8 with
  another bin).
- **Counting a module twice because both overlapping co-assemblies (main, marine) carry it.** It is the
  same population, at ≥ 98% identity.
- **Summing per-module bits for co-occurring novel modules** (207 bits for one event).
- **Weak shrinkage (C = 1) in the context model** (12-28% of targets worse than prevalence), and a
  Jaccard neighbour kernel as a phylogeny proxy.
- **PanGBank as the foundation of level 1**, a full PanGBank mirror (~200 GB), and sister-species
  pangenomes as species priors. **Ignoring projection's specific genes.**
- **EM-corrected training counts under the hierarchical prior** (model_spec §3.4 A as written). They
  leak the held-out bin back into its own prior (C10).
- **Per-family absence scoring for systems the catalogue splits** (atp, nuo). It understates or
  misplaces surprise; score the union of a system's families instead.
- **Level-1 presence surprise alone to separate contamination from accessory DNA** (accessory-like
  contigs score 38% > 5 bits, off-population 42%).
- **Beta(1,1) as the baseline for level-1 comparisons over a multi-species cluster universe.** It
  inflates the apparent gain from shrinkage by 55-60%. Compare against a fitted species-wide Beta(a, b).
  [pangbank_scale.md](experiments/pangbank_scale.md)
- **Species-level pangenomes for freshwater bins.** No freshwater species has more than 2 GTDB genomes.

## 5. Roadmap

Ordered by value per cost. "Settles" names the decision each item makes.

| # | Experiment | Cost | Settles |
|---|---|---|---|
| R1 | **Finish the cleaned baseline.** Cleaning of the main-assembly bins is done (finding 13a). Remaining steps: apply the proposed > 30%-removal guard and re-check lorbin_077/059; use BAMs to check per-base depth and reads spanning from the UBA9145 flagellar contigs (module_6/20/29) and the 4 Pontimonas low-ratio modules into core contigs; re-score the cleaned bins with block scoring (finding 22), the fragmentation absence term and catalogue-level surprise (absence bits per expected Mb after cleaning). Marine bins need coverage classes before they can be cleaned. | Local, about 1-2 CPU-h; BAMs exist | Whether the 18 remaining UBA9145 modules and the 58 modules in clean bins are biology (C1). Fixes the order of filters in the pipeline: coverage, then anchoring/content, then surprise. Checks whether the absence-contamination signal (20b) disappears once bins are cleaned. |
| R2 | **Merge rule on realistic bins.** The 200-species run (done) fixed the 95% family cut and the calibrated threshold. The next test repeats it with real MAG contig-length distributions (fragmentation model) and with MAG-built training pangenomes. It should calibrate the threshold from training genomes alone, since real use has no held-out isolates. It should also check whether the 95% cut, used for merge scoring, fragments families enough to hurt level-1 prevalence. Reuse the grex scripts. | grex, about 1 node-hour (the 200-species run took 11 min on one node) | Whether the 4-8% strain cost and 94-95% ANI reach survive MAG fragmentation and contamination in the references. Whether to use one family cut or two (95% for assignment, 80% for prevalence). |
| R3 | **Rebuild the Chesterfield catalogue on 3Di hits** (run on grex in progress). Compare family count, gene coverage, cross-phylum family purity and the sanity table against MMseqs2 mode. | About 1.5 GPU-h (running) plus about 15 min local | Whether 3Di becomes the default catalogue search. At GTDB scale it costs on the order of thousands of GPU-h at 31 proteins/s, so this decision needs evidence. |
| R4a | **Submit the pangbank_scale grex run** (scripts tested; 291 species for A, 40 genera for B; 455 PanGBank `.h5`, about 45 GB). | About 100-150 CPU-h, 2 genoa nodes, 30-60 min wall plus a 1-2 h download | Level 1: does genus shrinkage beat species-level empirical Bayes at 5-20 genomes, and if not, should model_spec §3 change (C13)? Level 2: how the context gain scales from N = 15 to 800, net of clade structure. This is the cheapest decisive test of both levels. |
| R4 | **E1: level-1 tree shrinkage against an external reference.** First fix the §3.4 A leave-one-out (C10). Then run the catalogue over GTDB species representatives, or over PanGBank GTDB_refseq pangenomes in 3-5 GTDB families with ≥ 3 genera, using complete isolates only. Hold out genomes, species and genera, and compare against species-only Beta(1,1) and global prevalence. Score Chesterfield bins against this prior. Read the Mattick occupancy-model full text first. Also add system-level scoring (union of a system's families, containment matching for truncated loci). | grex, tens of CPU-h, plus Spacedust at representative scale | Falsification test 1: does taxonomy carry usable signal outside the circular in-sample test (finding 20a)? Rank shrinkage vs the occupancy model (O1). |
| R5 | **Divergence floor for neighbourhood evidence.** Add a third query with curated or RNA-seq operons (e.g. Synechocystis, M. tuberculosis or a Bacteroidota), express every panel genome as AAI to the query, and find the AAI where precision reaches 0.9. Then recount catalogue locus support using only lineages beyond that floor. | Local, about 1 CPU-h | The panel rule (C3), and whether catalogue loci should require distant support (C4). |
| R6 | **Planted fragmentation test.** Cut bins' contigs further, rerun Spacedust and measure recovery of truncated clusters. | Local, under 1 CPU-h | β = 1 vs about 0.25, and how much of the end-depletion is a detection loss. |
| R7 | **Correct public-MAG prevalence** through the detection model, using CheckM2 per genome from PanGBank or GTDB metadata and public contig lengths. Re-score the survey. | Local, a few CPU-h | How much rare-present surprise comes from incomplete MAG/SAG references rather than our bins. |
| R8 | **E3: pooled level-2** with the level-1 offset and genus-grouped CV across many species (needs R4's pangenomes), plus a core-genome phylogeny covariate. | grex, 2-4 CPU-h per species | Whether context beyond island co-inheritance (exclusions, cross-species) is learnable, and whether level 2 is kept per cluster. |
| R9 | **E4/E6 on CAMI II marine and strain-madness.** Bin-level content surprise vs CheckM2 + GUNC; merge rule with a coverage term vs coverage alone. | grex, days of setup; known-membership data | Falsification tests 2 and 3: does content surprise add to quality control and binning at all? |
| R10 | **E5: population frequencies on CAMISIM strain mixtures** (2-5 strains, 5-20 samples), or deeper sequencing of the Chesterfield barcode04 sample. | grex, days | Whether φ_k is measurable at realistic depth, and the detection limit φ_min. This decides whether stage 4 is feasible on nanopore co-assemblies. |

## 6. Open questions for the PI

1. **Which unit?** The catalogue uses Spacedust/Leiden families, while the survey, context model,
   binning simulations and PanGBank test all used PPanGGOLiN modules. The two are not linked yet.
   Should we commit to one (O3), or map modules onto catalogue families?
2. **Is operon prediction still a goal?** Neighbourhood evidence tops out at about 0.88, and only +0.003
   over the gap rule for B. subtilis. Should operon recall stay a validation metric for the catalogue
   only?
3. **Should 3Di become the default at GTDB scale?** It improves coverage of distant genomes (mja
   4 → 15%), but costs thousands of GPU-h. Is that grex GPU budget available?
4. **Reference policy.** Pin GTDB r226 to match GTDB-Tk, or move to r232, which drops MAG-only species
   such as MAG-120531? Is building our own pangenomes for MAG-only species acceptable as the default?
5. **Depth.** Population-frequency claims need more than about 4× in one sample. Is deeper sequencing of
   the Chesterfield samples, or a time series, planned? Otherwise stage 4 should move to simulated or
   CAMI data first.
6. **Overlapping co-assemblies.** Should one assembly (main or marine) be canonical, or should bins be
   dereplicated across assemblies before any counting?
7. **Freshwater.** No species-level references exist. Are genus-level priors (sister PanGBank
   pangenomes, our own genus pangenomes) good enough for the questions you want to ask there?
8. **What should count as "surprising" biology?** The clean rare modules point to algal-polysaccharide
   and sugar-alcohol uptake across unrelated species. Is site-specific accessory content the target
   result, given that public references come from temperate coastal waters (a sampling bias the model
   cannot remove)?
9. **Benchmark.** CAMI II (strain-madness, marine) has gold standards now; CAMI III gut gold standards
   come later. Is CAMI II acceptable for E4-E6?
10. **Two family cuts?** Merge scoring wants 95%-identity families (finding 16). Prevalence and
    cross-species work want coarser ones: 80% for PanGBank and 50% for genus clusters. Is it acceptable
    to carry two family definitions, one for assignment and one for prevalence?
11. **Grex time.** pangbank_scale (about 100-150 CPU-h, 45 GB of downloads) is ready to submit and
    decides both levels (C13). Should it go ahead of the CAMI work?

## Addendum: theory documents (19:55)

Five derivations in docs/theory/ (phylo_prior, merge_test, population_frequency, calibrated_surprise, cluster_model). Corrections they make to the findings above:

- **C14 (calibrated_surprise):** catalogue_surprise's absence-flags-contamination AUC 0.85 is circular through CheckM2 contamination in epsilon (0.54 without it). The usable signal is calibrated presence-direction excess (AUC ~0.96, no CheckM2 input, 10 positives). Split cells by q < 1/2, not by observed state.
- **merge_test:** a gamma per-gene divergence model reproduces the grex failure ANIs (92.3/94.0/95.3-96.4% predicted vs 92-94/94-95/>=96% observed for cuts 0.80/0.90/0.95) and AUC within 0.01-0.03; an identity-based likelihood-ratio test needs ~7-26 genes per decision at the species boundary vs ~145 for the 0.80 cut.
- **phylo_prior:** the rank beta-binomial is a moment match to a gain/loss process on the tree (cov = pi(1-pi)e^{-mu d}); explains kappa ~3 (family) and ~0.5 (genus) and why genus shrinkage loses once a species has >= 5 genomes.
- **population_frequency:** precision is read-count limited (ell_eff ~3 kb); reproduces the >= 10 kb rule; coverage cannot separate f > 0.85 from 1; junction-spanning reads give linkage proof (decisive when abundance profiles are identical) and beat coverage only for elements < ~6 kb.
- **cluster_model:** operon adjacencies are lost at ~0.53 x the neutral rate in both E. coli and B. subtilis (r 0.97); breaks in close relatives are evidence against an operon, which is why close genomes lowered accuracy.

## Addendum: prediction tests (19:55)

- **merge_identity:** median per-gene divergence to the focal pangenome separates 94-95% ANI congeners from strains at AUC 0.994-1.00 (current surprise rule 0.64-0.77); single 50 kb contigs 0.87-0.99 (current rule at chance). merge_test eq. (8) needs a heavy-tail component (0.82-0.98 as written). Divergence is correlated over ~5-10 genes (variance inflation ~4.5).
- **operon_evidence:** cluster_model P2 fails as stated (summed independent evidence 0.682 E. coli, 0.593 B. subtilis vs best panels 0.878/0.877): evidence is strongly correlated across lineages (fitted weight 0.02-0.04); per-pair rate heterogeneity is needed, and the break term needs an independent homology search.
- **calibrated_check:** Z+ drops in the two heavily contaminated bins but across 45 cleaned bins barely beats removing random contigs of equal length (p = 0.074); AUC 0.955 vs CheckM2 not confirmed by coverage labels — CAMI needed. Chow-Liu correction works when a near-duplicate genome is in the reference (224 -> 29-42 bits), not for blocks seen nowhere else.
- **phylo_prior_check:** in Marinisoma, genus shrinkage helps only at m = 1 (crossover m* ~1-2, earlier than predicted ~3); the ancestral-state-first predictor is best at m = 1-3; at m >= 5 a two-state species-only model wins. P1 (distance beats rank) weakly supported on 49 placeable Chesterfield bins.

## Addendum: relatedness for level 1 (2026-11-20)

- **Decision: within a species, relatedness comes from gene content.** In the 200-species grex run
  ([tree_vs_context.md](experiments/tree_vs_context.md)) a UPGMA tree on accessory-family Jaccard (targets masked)
  was the best single relatedness predictor of held-out content for 62-65% of species at N = 25 and 97-98% at
  N ≥ 200, ahead of CDS-ANI and core-protein distance, and trees beat kernels by 9-25 points. Level 1 uses the
  gain/loss chain on a gene-content tree; sequence trees are kept as diagnostics. Context still adds 5.5-11
  points beyond all relatedness (families), growing with N, so level 2 stays.
- **Within species, GC-deviant genomes are gene-content outliers, not ANI outliers,** and look like MAG
  contamination (502 genomes; GC deviation 0.57 vs 0.11). Gene-content relatedness must therefore be computed
  on cleaned bins, or it will put contaminated MAGs on long branches.
- **Between species** (the PI's Marine Group I case: GC-shifted, sequence-tree outlier, close by gene content)
  is the open question. [xsp_tree_vs_content.md](experiments/xsp_tree_vs_content.md) compares the GTDB tree,
  core AAI and orthogroup content across 600 species in 18 clades. Result: between species the GTDB tree is best
  (13/18 clades, 64% of species) but only 0.3 points ahead of gene content (41.8 vs 41.5%); core AAI is worst; GC
  shift does not predict which wins. MGI-like cases exist (Taylorella equigenitalis) but are rare. Use the GTDB tree
  for cross-species priors, the gene-content tree within species.

## Addendum: scorer v0 and CAMI II marine (2026-11-21)

- `metapanclust/` packages the model (reference from PanGBank h5; gene-content tree + context; detection; calibrated
  excess; novelty against the species). [scorer_v0.md](experiments/scorer_v0.md): on 98 held-out genomes in 9
  species, novel genes detect 5% congener contamination at AUC 0.92 and the foreign contigs at 0.92; unit-level
  surprise does not (0.53-0.66); calibration is off (rare units over-predicted about 2×, likely from imprecise
  placement of incomplete queries).
- [cami_marine.md](experiments/cami_marine.md), 1,567 real binner bins with an oracle species reference: CheckM2
  flags contaminated bins far better (AUC 0.965 vs 0.846) and metapanclust adds +0.007 to it. metapanclust
  estimates completeness better (MAE 5.2 vs 8.2 points; CheckM2 +6.2 biased) and finds other-species contigs
  (AUC 0.84-0.87). Removing mostly-novel contigs cuts other-species contamination from a median 13% to 0.5% at
  3.2 points of completeness. Same-species (strain) contamination is invisible to novelty.
- **Re-scoped roadmap:** (1) genus-level references for the 78% of bins without a PanGBank species (GTDB tree
  between species); (2) novelty-based contig eviction as the first binning move, tested on all CAMI bins against
  binner output; (3) coverage/identity for strain mixtures; (4) recalibrate unit surprise on references degraded
  to the query's completeness before quoting surprise values.
