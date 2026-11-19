# Summed per-genome operon evidence: test of cluster_model.md P2 (2026-11-18)

Question: does scoring each adjacent same-strand pair by summed per-genome log-likelihood ratios
(phi = 0.53, nothing else fitted) over all reference genomes reach >= 0.878 accuracy for E. coli, where
the union rule falls to 0.689? Does a lineage (independence) correction help, and do panels spread over
lineages beat single-clade panels of the same size?

## Result

All accuracies are 5-fold block CV, except the union rules, which have fixed thresholds.

| E. coli (2,820 pairs, 68 genomes) | acc | prec | rec |
|---|---|---|---|
| gap prior only | 0.848 | 0.929 | 0.814 |
| union rule, all 68 | 0.689 | 0.667 | 0.986 |
| union rule, >= level 6 (best panel) | **0.878** | 0.908 | 0.891 |
| P2 as stated: gap prior + sum LLR, all 68, w=1 | **0.682** | 0.874 | 0.562 |
| same, each level collapsed to 1 genome (w = 1/n_level) | 0.858 | 0.903 | 0.860 |
| logistic [gap bins, sum LLR], w=1 / w=1/n_level | 0.866 / 0.866 | 0.93 | 0.85 |
| conserved terms only, levels >= 6, w=1 | 0.875 | 0.899 | 0.896 |
| breaks counted only when both genes have a homologue ('pres'), all 68 | 0.857 | 0.886 | 0.879 |

| B. subtilis (2,621 pairs, 92 genomes) | acc | prec | rec |
|---|---|---|---|
| gap prior only | 0.871 | 0.944 | 0.877 |
| union rule, all 92 (T1=100, T2=250) | 0.781 | 0.780 | 0.980 |
| union rule, >= level 7 (best panel) | **0.877** | 0.934 | 0.896 |
| P2 as stated, all 92, w=1 | **0.593** | 0.901 | 0.502 |
| w = 1/n_level | 0.852 | 0.928 | 0.867 |
| logistic [gap bins, sum LLR], w = 1/n_level | 0.876 | 0.936 | 0.892 |
| 'pres', all 92, w = 1/n_level | 0.874 | 0.944 | 0.881 |

- **P2 is falsified as stated.** With phi = 0.53 and independent genomes, the summed score over all genomes gives
  0.682 (E. coli) and 0.593 (B. subtilis). That is no better than the union rule and far below 0.878. No variant
  over all genomes reaches the best panel. The best is 0.868 for E. coli (logistic, conserved-only,
  level-collapsed) and 0.876 for B. subtilis (best panel 0.877).
- **What does hold:** collapsing each distance level to one effective genome lifts the all-genome score
  from 0.68 to 0.86 (E. coli) and from 0.59 to 0.85 (B. subtilis). That is 0.17-0.26 above the union rule on the same
  genomes. Close genomes then stop hurting, but they still add nothing over the distance-floored panel.
- **The evidence is overcounted 25-50x, not just at close range.** When a logistic is fitted on
  [gap bins, summed LLR], the LLR gets coefficient 0.04 (E. coli) and 0.02 (B. subtilis) with w=1. With
  level-collapsed weights it gets 0.25 and 0.24. A calibrated score would get 1. So even the 5-8 levels left after
  collapsing behave like about 1-2 independent observations. Conservation is strongly correlated across distant
  lineages too: operon age and strength vary per pair, so one genome's conservation predicts the others'.
  A single global phi with independent genomes misses this. The theory's per-pair heterogeneity in
  p_i (its section on recall curves) has to enter the evidence model too, for example as a per-pair latent
  phi or a beta-binomial over lineages, and not only the recall model.
- **The break term is the main failure in P2 as written, but presence data could not fix it here.** Without
  homologue presence, every absent gene counts as a break (ln-ratio about -0.07 per distant genome).
  Summed over about 50 distant genomes, that swamps the positive evidence: E. coli levels >= 4 drop to 0.656. Restricting the
  sum to the conserved terms gives 0.744 for all genomes and 0.875 for levels >= 6. The proxy used for
  presence (the gene appears in any Spacedust cluster of that run) is circular: a gene lands in a cluster almost only
  when its neighbour is conserved. As a result, P(c | non-TU, both 'present') is 0.85-0.9 at distant levels
  (x|present = 0.15-0.2 against 3-5 unconditional). A real test of the "breaks are negative evidence" claim
  needs a homology search (e.g. MMseqs2 query vs each genome) that does not depend on the clusters.
- **Adding levels 0-3 (the stated falsification test):** under the literal score with w=1, adding them changes
  E. coli accuracy from 0.656 (>= level 4) to 0.682 (all), so it does not lower it. But both values are below
  gap-only. Under the conserved-only score, adding them lowers accuracy from 0.875 (>= 6) to 0.744. Under
  'pres', it lowers it from 0.865 to 0.857.
- **Spread vs single-clade panels (n = 3-4, 40 draws each):** accuracy cannot tell them apart. All E. coli panels
  give 0.858-0.866 and all B. subtilis panels give 0.858-0.875, under both the LLR score and the union rule. Conserved-pair recall
  for the spread panels is above the mean of the single-clade panels: E. coli n=4 gives 0.295 against a mean of 0.257
  (clades 0.177-0.314), and B. subtilis distant clades give 0.24 against 0.14-0.29. It is below the nearest single clade
  (Gamma 0.314, other Bacillota 0.293), so distance and spread are confounded at this size. The data give weak support
  for the effective-lineage view and do not support "spread beats single-clade" on accuracy.

## Methods

- `experiments/operon_evidence/evidence.py TAG PANEL MATRIX RUN_DIR BAKTA OUT_DIR`. Inputs are
  `phylo_scaling/results/ext_matrix.tsv` and `phylo_scaling_bsub/results/bsub_matrix.tsv` (conserved = value >= 1),
  with panel levels from `panel_ext.tsv` / `panel_bsub.tsv`.
- Per genome k: x_k = -ln P(conserved | different-TU pair), estimated on the training folds as in
  `theory_cluster/check.py`. Conserved LLR = (1-phi) x_k. Break LLR = ln((1-e^{-phi x_k})/(1-e^{-x_k})).
  phi = 0.53.
- Classifier 'prior + LLR' has no fitted weight. The logit is the empirical gap-bin log-odds (18 bins, from the training folds)
  plus the summed score, and a pair is called when the logit > 0. The 'logistic' variant fits gap-bin indicators and one
  coefficient on the score (L2 1e-2). 5 folds of contiguous pairs.
- Variants: `theory` (every non-conserved genome is a break), `pres` (break only if both genes appear in
  some cluster of that run, and x conditioned on that; circular proxy, see above), `conserved-only` (no break term).
  Weights: w=1, or w = 1/(genomes in level).
- Union rule: `gap <= T1 OR (conserved in any panel genome AND gap <= T2)`, with E. coli 50/300 and B. subtilis 100/250
  (the tuned values from phylo_scaling_bsub).
- Panels: E. coli clades L5 Gamma, L6 Beta, L7 Alpha, Bacillota (Bsu, Sau, Lmo, Cac); spread = one genome from each of n
  distinct levels in {5,6,7,8}. B. subtilis clades L3, L5, L7, L8; spread from {4,...,8}.
- Compute: about 1 CPU-min in total (nice 15, 1 thread).

## Files

- Code: `/data/dev/metapanclust/experiments/operon_evidence/evidence.py`
- Output: `/data/scratch/metapanclust_data/operon_evidence/{eco,bsub}_results.txt` (all rows, per-level x),
  `{eco,bsub}_presence.npy` (pair x genome presence proxy), 0.5 MB.

Reproduce:
```
D=/data/scratch/metapanclust_data; E=/data/dev/metapanclust/experiments
nice -n 15 python3 $E/operon_evidence/evidence.py eco $E/phylo_scaling/panel_ext.tsv $D/phylo_scaling/results/ext_matrix.tsv \
  $D/phylo_scaling/runs/pw /data/dev/metacyc_benchmark/ecoli/mag_out/annotation/bakta/basic/annotation.tsv $D/operon_evidence
nice -n 15 python3 $E/operon_evidence/evidence.py bsub $E/phylo_scaling_bsub/panel_bsub.tsv $D/phylo_scaling_bsub/results/bsub_matrix.tsv \
  $D/phylo_scaling_bsub/runs/pw /data/dev/metacyc_benchmark/bsub/mag_out/annotation/bakta/basic/annotation.tsv $D/operon_evidence
```
