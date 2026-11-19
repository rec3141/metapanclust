# Identity-based merge test on the same-genus bins (2026-11-18)

Tests the prediction in `docs/theory/merge_test.md` §3: an identity likelihood-ratio test should reach
AUC ≥ 0.95 up to 95% ANI for 20%-of-genome merges and beat the surprise-per-Mb rule. It uses the
merges from `binning_sim_genus.md`, regenerated with the same seed and checked to be bit-identical
(992 merges per focal species, same partners, same `excess` values).

**Short answer: the prediction holds, but only with a robust divergence statistic.** At 94-95% ANI:

| Statistic (on the added genes) | S. baltica + S. hafniensis (94.7%) | P. distincta + P. arctica (94.0%) | P. arctica + P. distincta (94.0%) |
|---|---|---|---|
| surprise excess (current rule) | 0.64 | 0.65 | 0.77 |
| hard-cut foreign count, c = 0.80 | 0.78 | 0.69 | 0.83 |
| hard-cut foreign count, c = 0.95 | 0.87 | 0.99 | 0.98 |
| identity LR, eq. (8), raw divergence | 0.82 | 0.98 | 0.93 |
| identity LR, divergence clipped at 0.10 | 0.90 | 1.00 | 1.00 |
| mean divergence x̄ | 0.86 | 0.99 | 0.97 |
| mean of clipped divergence | 0.94 | 1.00 | 1.00 |
| **median divergence** | **0.994** | **1.00** | **1.00** |

The values are AUCs, pooled over 5-40% added and both bases. At 20% added alone, the median
divergence gives 0.997 / 1.00 / 1.00. Rows 2-3 of the table are listed by base genome first
(focal + contaminant).

- **Rejection at 5% false rejects of strains.** For baltica + hafniensis, the share of contaminant
  merges rejected is 0.04 by excess, 0.29 by the c = 0.95 count, 0.09 by the raw LR, 0.49 by the
  clipped LR and **0.97 by the median divergence**. For the Pseudoalteromonas pair it is 0.16-0.26
  by excess and 1.00 by the median.
- **Why the raw LR falls short.** The gamma model in eq. (8) does not fit the tail. Some best hits
  are distant (paralogs, xenologs, partial families at 30-80% identity), and their divergence
  dominates the sum T. The mean strain divergence μ_s is 0.007-0.026 raw but 0.006-0.017 once it is
  clipped at 0.10. The LR also adds the family-less count, which carries little information at
  94-95% ANI and is clustered in islands (see below). Robust location statistics on the shared genes
  (the median, or the clipped mean) recover what the theory predicts.
- **Single 50-kb contigs** (about 44 genes) are harder. Against strain contigs at 94-95% ANI, the
  median divergence gives AUC 0.87-0.99, the clipped LR 0.83-0.95, the c = 0.95 count 0.74-0.89 and
  excess 0.50-0.55, which is chance.

### By ANI bin

The mean over pairs, with the minimum in brackets:

| ANI bin | Pairs | AUC excess | AUC raw LR | AUC x̄ | AUC c = 0.95 count | AUC c = 0.80 count | Rejected at 5% strain FR: excess / raw LR / x̄ | Single-contig AUC: excess / raw LR |
|---|---|---|---|---|---|---|---|---|
| no skani (< ~83%), incl. far | 33 | 1.00 (0.99) | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 / 1.00 / 1.00 | 0.99 / 0.99 |
| 83-85% | 6 | 0.96 (0.84) | 1.00 (0.98) | 1.00 | 1.00 | 0.99 | 0.86 / 0.98 / 1.00 | 0.91 / 0.97 |
| 85-90% | 4 | 0.99 (0.97) | 1.00 (0.99) | 1.00 | 1.00 | 1.00 | 0.96 / 0.99 / 1.00 | 0.87 / 0.93 |
| 90-93% | 2 | 0.90 (0.84) | 0.95 (0.91) | 0.99 (0.98) | 1.00 (0.99) | 0.95 | 0.64 / 0.77 / 0.95 | 0.76 / 0.93 |
| 93-95% | 3 | 0.68 (0.64) | 0.91 (0.82) | 0.94 (0.86) | 0.95 (0.87) | 0.77 | 0.15 / 0.53 / 0.69 | 0.52 / 0.81 |

At ≥ 90% ANI, the median divergence has AUC ≥ 0.99 in every pair. The clipped LR is ≥ 0.95 in
every pair except baltica + hafniensis, where it is 0.90. The full table is in
`pairs.tsv`.

**At 94-95% ANI, the contaminant's mean divergence x̄ is 2-2.5 times a strain's.** The x̄ medians
are 0.035 against 0.014 for baltica + hafniensis, and 0.045-0.048 against 0.021-0.022 for the
Pseudoalteromonas pairs. This fits the theory's κ(1-A) ≈ 0.033-0.039.

## How clustered divergence is along contigs

These figures cover genes in contig order, for each (genome, focal pangenome) combination; the values
are medians. φ is the variance of block means divided by the per-gene variance over B, so it is the
variance inflation for blocks of B consecutive genes. The effective number of independent genes is
then m_eff ≈ m/φ.

| Group | ACF of x at lag 1 / 2 / 5 / 10 / 20 genes | φ_x at B = 10 / 25 / 50 | φ for family-less indicator at B = 10 / 25 / 50 |
|---|---|---|---|
| same species (held-out) | 0.47 / 0.30 / 0.12 / 0.05 / 0.02 | 2.9 / 4.1 / 4.6 | 1.9 / 2.6 / 2.6 |
| congener ≥ 90% ANI | 0.44 / 0.30 / 0.15 / 0.09 / 0.04 | 2.9 / 3.9 / 4.5 | 3.3 / 5.3 / 6.8 |
| congener 83-90% | 0.32 / 0.21 / 0.11 / 0.05 / 0.01 | 2.5 / 3.2 / 3.5 | 3.4 / 5.4 / 7.4 |
| congener < 83% | 0.32 / 0.23 / 0.13 / 0.07 / 0.03 | 2.6 / 3.7 / 4.1 | 3.2 / 4.8 / 6.2 |

- **Divergence is clustered over about 5-10 genes.** The autocorrelation decays to about 0.1 by lag 5,
  and φ_x levels off at about 4-4.5 by 50-gene blocks. That sits inside the theory's guess of
  φ_x ≈ 3-10, at the low end.
- **Gene presence/absence is more clustered in congeners.** φ reaches 6-7 at 50 genes and is still
  rising, consistent with islands. Divergence is the less correlated signal, as the theory assumed.
- **What this means for sample size.** A 20% merge adds about 820 genes, so m_eff ≈ 180, far above
  the 7-26 independent genes the theory needs at 95% ANI. A 50-kb contig adds about 44 genes, so
  m_eff ≈ 10. That is near the theory's limit, and matches the single-contig AUCs of 0.87-0.99.

## Method

Code: `experiments/merge_identity/`:

- `01_simulate.py` rebuilds the binning_sim_genus merges and scores the added genes.
- `02_run.sh` runs it for all 8 focal species.
- `03_analyze.py` computes the AUCs and the clustering tables.

Data: `/data/scratch/metapanclust_data/merge_identity/` (`sim_*.tsv`, `autocorr_*.tsv`, `pairs.tsv`,
`summary.txt`, 2.6 MB). The run took about 2 CPU-minutes at 1 thread under `nice -n 15`.

- **Per-gene input.** Each gene's best mmseqs hit to the focal family representatives (30% identity,
  50% coverage). Divergence is x = 1 − identity. A gene with no hit counts as family-less.
- **Representatives, not nearest members.** Identity is to the family representative, not to the
  nearest member as the theory recommends. That inflates μ_s, so the separation here is a lower
  bound.
- **Raw LR.** Eq. (8), per added gene, with α = 2 and κ = 0.65. H0 (g_s, μ_s) is fitted on the
  focal species' held-out genomes, excluding the base genome and its partner. H1 is a congener at
  95% ANI: μ_A = max(0.0325, 1.5 μ_s) and g_A = g_s − 0.08.
- **Clipped LR.** The same, with x clipped at 0.10 and μ refitted on the clipped values.
- **AUC and rejection.** The AUC is computed as in binning_sim_genus: contaminant against strain
  merges, at equal fraction and base. The 5%-false-reject threshold is the 95th percentile of the
  strain merges in the same cell. That is an oracle calibration, so the rejection rates are
  optimistic for every statistic alike.

## Caveats and next steps

- **The pass rests on the median.** The theory's prediction holds for robust divergence
  statistics (median, clipped mean). The literal gamma-model LR fails on baltica + hafniensis (AUC
  0.82) because the per-gene divergence distribution has a heavy tail that the model leaves out.
  The model needs a tail component, such as a mixture with a broad "distant best hit" term, or a
  restriction to persistent families. Either would probably give a proper LR that matches the
  median.
- **Not tested against within-species tails.** The held-out strains here are 97-100% ANI to each
  other. A focal species with members at 95-96% ANI (the species-concept limit) was not tested.
- **Only 3 hard pairs**, all from 2 genera.
