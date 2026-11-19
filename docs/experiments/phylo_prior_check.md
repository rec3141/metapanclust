# Quick check of phylo_prior.md predictions P1 and P2 (2026-11-18)

Status: P2 tested on Marinisoma only (3 focal species, 49 genomes). The Flavobacterium run hit the 15 CPU-min
cap at its 8th focal species and was killed before writing output. P1 tested on 49 of the 93 Chesterfield
bins, with few pairs per rank. Both results are small-sample indications, not results.

## P2: species-only EB vs genus prior vs ancestral-state-first, Marinisoma

Held-out bits per test genome, summed over the clusters in each stratum. Values are means over 3 focal species
x 5 random splits. Up to 50 test genomes per split; Marinisoma species have only 15-17 genomes, so 5-14 in
practice. Wins are counted over the 3 focal species.

| Stratum | m | sp_eb | shrunk (EB κ=2) | mix | sp_2state | **anc** | shrunk wins vs sp_eb | anc wins vs sp_eb | anc wins vs shrunk |
|---|---|---|---|---|---|---|---|---|---|
| all | 1 | 2,488 | 1,406 | 1,646 | 1,663 | **1,322** | 3/3 | 3/3 | 3/3 |
| all | 2 | 1,298 | 1,309 | 1,219 | 1,297 | **1,196** | 1/3 | 3/3 | 3/3 |
| all | 3 | 1,202 | 1,224 | 1,143 | 1,202 | **1,129** | 1/3 | 3/3 | 3/3 |
| all | 5 | 1,107 | 1,155 | 1,091 | **1,073** | 1,079 | 0/3 | 3/3 | 3/3 |
| all | 10 | 976 | 1,053 | 993 | **941** | 981 | 0/3 | 1/3 | 3/3 |
| shared | 1 | 1,435 | **749** | 901 | 1,095 | 759 | 3/3 | 3/3 | 2/3 |
| shared | 2 | 800 | 717 | **711** | 784 | 721 | 3/3 | 3/3 | 1/3 |
| shared | 5 | 671 | 636 | **633** | 652 | 659 | 3/3 | 2/3 | 1/3 |
| shared | 10 | 561 | 556 | 553 | **547** | 567 | 3/3 | 1/3 | 1/3 |
| focal_only | 1 | 429 | 491 | **370** | 484 | 414 | 0/3 | 2/3 | 3/3 |
| focal_only | 5 | 393 | 455 | 405 | 365 | **360** | 0/3 | 3/3 | 3/3 |
| focal_only | 10 | 389 | 461 | 408 | **349** | 364 | 0/3 | 3/3 | 3/3 |
| congener_only | 1 | 624 | 166 | 375 | **84** | 148 | 3/3 | 3/3 | 2/3 |
| congener_only | 5 | **44** | 64 | 53 | 56 | 59 | 0/3 | 1/3 | 1/3 |
| congener_only | 10 | **26** | 37 | 31 | 44 | 49 | 0/3 | 1/3 | 1/3 |

Predictors:

- `sp_eb`, `shrunk`, `mix`: as in `pangbank_scale.md`, using `worker.py` functions unchanged.
- `sp_2state` and `anc` are new. Both use the same two-state form. For each cluster, the species'
  prevalence is drawn from one of two Betas:
  - "present in the species ancestor": Beta(a1, b1), prevalence near 1;
  - "absent": Beta(a0, b0), prevalence near 0.
- Given the m training genomes, the predictor computes the posterior state weight w from the
  beta-binomial likelihoods. The prediction is w·(x+a1)/(m+a1+b1) + (1−w)·(x+a0)/(m+a0+b0). This
  form can jump: one observation can flip the state.
- `anc`: Pr(present state) = sigmoid(c0 + c1·logit(p_genus) + c2·fraction of congeners carrying the
  cluster).
  - All 7 parameters are fitted by maximum likelihood on the congeners only. Each congener is scored on
    its full data against features computed from the remaining congeners, so the focal species is never
    used.
- `sp_2state`: the same mixture with a constant Pr(present), all 5 parameters fitted on the m training
  genomes. It is the "no genus information" control.

**Reading:**

- **The crossover is real, and it is earlier than predicted.**
  - Genus shrinkage beats `sp_eb` overall only at m = 1 (1,406 vs 2,488 bits, 3/3).
  - It already loses at m = 2-3 (1/3) and at m ≥ 5 (0/3).
  - m* ≈ 1-2, below the predicted ≈ 3. The fitted κ here is 2, not 0.5.
  - In the shared stratum, shrinkage keeps winning up to m = 10 (3/3), but by a margin that shrinks
    from 48% to 1%.
- **The ancestral-state-first predictor is the best valid predictor overall at m = 1-3, and beats
  `shrunk` at every m.**
  - Overall bits: 1,322 / 1,196 / 1,129 vs the best of the others at 1,406 / 1,219 / 1,143.
  - It beats `sp_eb` at m = 1-5 (3/3). At m = 10 it wins 1/3, mean 981 vs 976.
  - It beats `shrunk` 3/3 at every m.
  - In focal_only it matches or beats `sp_eb` at every m, as predicted.
  - In congener_only it loses to `sp_eb` at m ≥ 5 (59 vs 44). That part of P2 ("match or beat in
    congener_only at every m") is not met.
- **Most of the gain at m ≥ 5 comes from the two-state form, not from the genus.**
  - `sp_2state` uses no genus information, yet it is best overall at m = 5 and 10 (1,073 and 941 bits,
    against 1,107 and 976 for `sp_eb`).
  - At m ≥ 5, adding the genus to the two-state prior (`anc`) is neutral to slightly harmful
    (+6 and +40 bits). This is the screening-off prediction.
  - At m = 1-3 the genus helps `anc` by 70-340 bits over `sp_2state`.
- So P2 is supported in shape: genus information pays at small m, the conspecifics screen it off by
  m ≈ 5, and a discrete-state predictor beats the linear pseudo-count. The crossover m* and the
  congener_only claim are not confirmed.
- Caveats:
  - 3 focal species from one MAG-only genus, with 2 congeners each. In the `anc` fit, each congener is
    scored against features from a single remaining species.
  - Completeness is not modelled, and these are MAGs.
  - The `sp_eb` grid fit is weak at m = 1. Its 2,488 bits are partly a fitting artefact; `sp_2state` at
    m = 1 gets 1,663.
- Flavobacterium (8 species, 27,644 genus clusters) did not finish. The fitted `anc` components for its
  first 7 focal species were consistent with each other:
  - a1/b1 ≈ 3.7/0.15 (present-state mean prevalence ≈ 0.96);
  - a0/b0 ≈ 0.03/0.5-2 (absent-state mean ≈ 0.02-0.05);
  - κ = 0.5, as in pangbank_scale.
  - Cost was about 1.3 CPU-min per focal species, almost all of it in L-BFGS with numerical gradients
    over 7 × 27k beta-binomial terms. Analytic gradients, or a smaller mmseqs cluster set, would make the
    full run cheap. The script also needs to write per-focal output incrementally.

## P1: within-rank tree distance vs concordance, Chesterfield

- Distance: patristic distance in GTDB r226 `bac120_r226.tree` between each bin's GTDB-Tk closest
  reference (the ANI-screen reference, otherwise the placement reference).
  - 49 of the 93 bins have one. The other 44 are placed only by topology, and no classify tree was kept.
  - Pairs that share a reference are dropped.
- Concordance: the completeness-corrected Pr(seen_j | seen_i) of `rank_decay.py`.
- Pairs within each deepest shared rank are split at the median distance.

| Rank (pairs) | d range (median) | 2-3 bins near / far | 4-9 | 10-29 | ≥ 30 |
|---|---|---|---|---|---|
| phylum (67) | 2.05-2.69 (2.36) | 0.03 / 0.02 | 0.06 / 0.09 | 0.43 / 0.42 | 0.88 / 0.83 |
| class (131) | 1.24-2.23 (1.65) | 0.07 / 0.04 | 0.30 / 0.21 | 0.70 / 0.65 | 0.93 / 0.94 |
| order (25) | 0.70-1.26 (1.06) | 0.13 / 0.11 | 0.50 / 0.68 | 0.91 / 0.91 | 1.21 / 1.04 |
| family (27) | 0.29-0.99 (0.71) | 0.35 / 0.62 | 0.85 / 0.75 | 0.88 / 0.86 | 0.96 / 0.96 |
| genus (10) | 0.07-0.35 (0.24) | 0.36 / 0.24 | 0.87 / 0.76 | 0.95 / 1.08 | 0.89 / 1.22 |

**Reading:**

- **Within a rank, near beats far in 12 of 20 cells, with 2 ties and 6 losses.**
  - The clearest case is class, the only rank with more than 100 pairs: near is ahead in the 3 lower
    strata (0.07 vs 0.04, 0.30 vs 0.21, 0.70 vs 0.65).
  - Order, family and genus have 10-27 pairs each and are noise-dominated.
  - The within-rank P1 prediction is weakly supported, not established.
- **Across ranks, concordance falls on one decreasing curve in d.** For the 4-9 stratum:

  | d | 0.24 | 0.71 | 1.06 | 1.65 | 2.36 |
  |---|---|---|---|---|---|
  | concordance | 0.83 | 0.80 | 0.64 | 0.26 | 0.08 |

- **Genus pairs are much closer than family pairs** in this subset: median d 0.24 vs 0.71, and the
  ranges barely overlap.
  - This rules out the first explanation in phylo_prior.md §5.1 for "genus no more concordant than
    family", namely that genus pairs are not closer. On these 49 bins genus pairs are as or more
    concordant: 0.33 vs 0.49 (2-3 bins), 0.83 vs 0.80, 1.02 vs 0.87.
  - The remaining explanation (the cross-genus construction filter) is untested.
- Caveats:
  - The reference-to-reference distance is a proxy for the bin-to-bin distance.
  - The 44 dropped bins are the novel ones.
  - Pairs share bins, so they are not independent.

## Files and cost

- Code: `experiments/phylo_prior_check/`
  - `prep.py`: reuses `worker.py` run_B up to the genus-cluster matrices;
  - `p2_eval.py`, `summ.py`, `p1_distance.py`.
- Data: `/data/scratch/metapanclust_data/phylo_prior_check/`
  - `Marinisoma.npz`, `Flavobacterium.npz`, `p2_Marinisoma.tsv`, `p1.tsv`, `prep.log`;
  - 220 KB in total.
- Container: `metapanclust-tools:0.1`, 1 CPU, nice 15. No downloads.
- About 15 CPU-min in total:
  - prep 5 min (of which Flavobacterium 4.3 min);
  - P2 eval: Marinisoma 0.7 min, Flavobacterium 9.3 min (killed);
  - P1: 3 s.
