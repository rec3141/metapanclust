# Calibrated surprise: independent checks of Z+ and the Chow–Liu correction

Tests three predictions of `docs/theory/calibrated_surprise.md` against the coverage-based cleaning in
`bin_cleaning.md`. That cleaning labels contaminant contigs without CheckM2.
Code: `experiments/calibrated_check/`. Outputs: `/data/scratch/metapanclust_data/calibrated_check/`.
Run time is about 4 CPU-min in total. Output size is under 1 MB.

## Results

1. **Z+ falls when cleaning removes DNA, but mostly because DNA was removed.** Across the 45 changed
   bins, ΔZ+ correlates with removed Mb at ρ = −0.62 (p = 5e-6). A size-matched null removes random
   contigs with the same total bp, 200 draws per bin. Against that null, the real cleaning lowers Z+ only
   slightly: median −0.07 SD units, 27/45 bins below the null mean, Wilcoxon p = 0.07. Of the 28 bins
   with completeness ≥ 0.5, only 2 have P < 0.05.
2. **The five bins that cleaning changed most:**

   | bin | removed bp | Z+ orig → clean | ΔZ+ | null ΔZ+ (sd) | z vs null | CheckM2 ΔConta |
   |---|---|---|---|---|---|---|
   | TMED14 maxbin_102 | 48% | 5.79 → 1.86 | −3.93 | −2.41 (0.87) | **−1.75** (P 0.05) | −16.9 |
   | Pseudothioglobus vamb_711_sub | 29% | 3.63 → −0.40 | −4.03 | −0.92 (2.33) | −1.33 (P 0.06) | −12.9 |
   | Flavobacterium lorbin_015 | 19% | 21.76 → 19.86 | −1.90 | −2.60 (0.81) | +0.87 | −7.2 |
   | UBA9145 maxbin_040 | 8% | −1.58 → −1.53 | +0.05 | −0.20 (0.28) | +0.88 | −1.5 |
   | Nitrosopumilus maxbin_115 | 48% | −7.28 → −3.73 | +3.55 | +0.86 (0.50) | **+5.35** | −6.0 |

   - The two bins with the largest CheckM2 contamination drops (102, 711) behave as predicted.
   - lorbin_015 keeps Z+ ≈ 20 after cleaning. Its high Z+ is not explained by the removed contigs.
   - maxbin_115 moves the wrong way. Its removed contigs carry less presence excess than its kept ones
     (178 vs 375 bits/Mb).
   - In the 48 unchanged bins, ΔZ+ is exactly 0 by construction.
3. **Original Z+ as a predictor of the coverage label is weak.** Z+ predicts "cleaning removed ≥ 10% of bp"
   with AUC 0.67 (compl ≥ 0.7, only 6 positives). Against CheckM2 contamination ≥ 10% it gives AUC 0.955,
   which reproduces the theory doc's figure.
4. **Per contig, removed contigs carry no more presence excess than kept contigs of similar length.**
   - Comparison: E+ per kb, summed over cells with q(seen) < 1/2 where the family is present. The excess
     of each cell is split evenly over its loci. Each removed contig is compared with kept contigs of the
     same bin that are within 2× its length.
   - Across 130,563 such pairs, P(removed > kept) = 0.36 with 30% ties, giving AUC 0.514.
   - Taking each removed contig minus the median of its matched kept contigs, the median difference is 0.
     The median difference is > 0 in 23 of 45 bins.
   - The effect is limited to long contigs. Among contigs ≥ 20 kb, the median E+/kb is 0.68 for removed
     (n = 35) vs 0.33 for kept (n = 1677). Below 20 kb there is no difference.
   - At bin level, the removed DNA is richer in excess than the kept DNA in maxbin_102 (750 vs 516 bits/Mb),
     vamb_711_sub (686 vs 256), lorbin_015 (616 vs 497) and maxbin_040 (394 vs 338). It is poorer in
     maxbin_115 (178 vs 375).
5. **Chow–Liu chain-rule code on UBA9145 modules (PPanGGOLiN, 54 genomes, leave-one-genome-out).** All
   values are bits, given as independent → Chow–Liu.

   | setting | maxbin_040 modules (K=105) | maxbin_040 module families (K=527) | all genomes, median CL/indep (modules / families) |
   |---|---|---|---|
   | twin genome in reference, +½ per 2×2 cell | 224.5 → 41.5 | 1217.5 → 139.7 | 0.77 / 0.37 |
   | twin in reference, pair prior = product of marginals | 224.5 → 28.7 | 1217.5 → 57.5 | 0.83 / 0.37 |
   | both Chesterfield genomes held out, +½ per cell | 282.4 → 71.7 | 1546.2 → 312.2 | 0.86 / 0.38 |
   | both held out, pair prior = product of marginals | 282.4 → **270.5** | 1546.2 → **1494.5** | 0.86 / 0.38 |

   - "Twin" is chesterfield_marine_maxbin_019, a near-duplicate of maxbin_040. maxbin_040 has 55 modules,
     and 36 of them occur only in these two genomes.
   - When the reference contains the block, the tree code shrinks maxbin_040 by 5–20×. Its rank drops
     from 1st to 14th–35th of 54 genomes, which supports the overcount argument in §2.2.
   - When nothing in the reference has seen the block, a pair prior with no built-in association gives
     back ≈ the independent sum (96–97%). This is the "no data" case of §2.2: correcting it needs catalogue
     adjacency or contig co-location, not within-species co-occurrence.
   - The +½-per-cell smoothing looked like a correction in that case (282 → 72), but it is an artifact:
     the prior says a child is present 50% of the time whenever its parent is.
   - Over typical genomes the correction is modest for modules (median 0.77–0.86) and large for module
     families (0.37), because families within one module are almost redundant.
   - The pmi form (s_ind − Σ pmi) is not a proper code when the single and pair marginals are smoothed
     differently, and it went negative in one setting. The conditional chain is the number to use.

## Methods and conventions

- **Prior.** Not refitted. The original leave-one-bin-out r_j = 0 fit (`fit_flat.npz` `p_r0`, `eps_r0`) is
  used for both the original and the cleaned bins. A cleaned bin's own changes therefore do not feed back
  into its prior. The other bins' priors are also unchanged, apart from a negligible leave-one-out effect.
- **Cleaned observations.** The I/P/A observations are rebuilt from `loci_span.tsv`, keeping only loci on
  contigs that `decisions.tsv` marks `keep`.
- **Detection.** P_intact, P_lost and FS are recomputed on the kept contig lengths. The genome size
  G_j = ΣL/c_j is held at its original value, so the effective completeness falls by the removed fraction.
  No CheckM2 is used anywhere in the cleaned scoring.
  - Removing DNA under fixed G lowers c. That makes the remaining presences more surprising, which partly
    offsets the loss of seen cells. The size-matched null uses the same convention, so the comparison
    against the null is fair.
- **Z+.** Defined as in theory §1.1: the 3-state excess over cells with q(seen) < 1/2, divided by the square
  root of their variance. `per_bin_calibrated.tsv` / check4 is reproduced (AUC 0.955).
- **Chow–Liu.** For each held-out genome, mutual information is computed on smoothed 2×2 tables from the
  other genomes. The tree is the maximum spanning tree rooted at node 0. The code length is
  s_T = −log2 Pr(o_root) − Σ log2 Pr(o_k | o_parent). Module presence means completion ≥ 0.5 in
  `modules_in_genomes.tsv`. Module families are the 527 families in `functional_modules.tsv`, taken from
  `gene_presence_absence.Rtab`.

## Interpretation

- Z+ survives the coverage label only for the two strongly contaminated bins (maxbin_102, vamb_711_sub),
  and only weakly beyond the size effect (z ≈ −1.3 to −1.8).
- The coverage label as a whole is not well predicted by Z+ (AUC 0.67), and removed contigs are not
  enriched in presence excess (AUC 0.51).
- So the AUC 0.955 against CheckM2 is not independently confirmed by coverage-based cleaning. Either Z+
  and CheckM2 share a bias (both respond to foreign marker-like or rare content), or the cleaning labels
  are noisy: many removals are short contigs with no loci, which are 30% ties.
- CAMI ground truth remains the needed test.

## Files

`zplus_cleaned.tsv`, `zplus_null.tsv`, `contig_excess.tsv`, `chowliu_*.tsv` and logs in
`/data/scratch/metapanclust_data/calibrated_check/`.
Scripts in `experiments/calibrated_check/`:
- `zplus_cleaned.py`
- `zplus_null.py`
- `chowliu_uba9145.py`: set the environment variables `LEAVE_BOTH_OUT=1` and `PAIR_PRIOR=indep` to select
  the variants.
