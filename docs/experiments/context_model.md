# Level-2 "context" surprise vs prevalence-only surprise (2026-11-18)

## Result

- **Context beats prevalence by about 30% of the bits, in both species, out of sample.** Leave-one-genome-out
  over the public genomes, an L2 logistic regression on accessory-family presence (the target's own module
  masked, completeness as a covariate, C = 0.1) predicts module presence at:

  | Mean bits per (genome, target), LOO | UBA9145 (52) modules | UBA9145 families | MAG-120531 (67) modules | MAG-120531 families |
  |---|---|---|---|---|
  | Prevalence, Beta(1,1) | 0.583 | 0.777 | 0.529 | 0.753 |
  | Prevalence, completeness-aware (P = q·c) | 0.559 | 0.757 | 0.511 | 0.733 |
  | Jaccard-weighted neighbours (phylogeny proxy, best γ = 2) | 0.526 | 0.716 | 0.490 | 0.697 |
  | Logistic on other modules (C = 1) | 0.443 | 0.670 | 0.444 | 0.648 |
  | **Logistic on accessory families (C = 0.1)** | **0.395** | **0.584** | **0.374** | **0.527** |
  | same, in-sample (overfitting check) | 0.073 | 0.111 | 0.109 | 0.163 |

  Targets: 69 / 40 modules with 2 to n−2 public carriers, plus 200 random accessory families with 5 to n−5
  carriers. The context model is better for 74% / 68% of modules (> 0.05 bits/genome) and worse for 6% / 5%.
  Hyperparameters picked on one species are also the best on the other (all six model-family choices
  transfer, gain lost ≤ 0.016 bits), so the table is not tuned to the test set.
- **Completeness alone is worth little (0.02 bits); the neighbour/phylogeny proxy is weak (0.04-0.06).**
  Accessory-family Jaccard between conspecific MAGs is low, so sharp kernels (γ ≥ 8) degenerate to the
  prior. A real phylogeny (core-SNP tree or ANI) has not been tried.
- **Which modules are predictable: co-inherited pieces of one island, never alternatives.** 41/69 UBA9145
  and 14/40 MAG-120531 modules gain > 0.1 bits/genome (top gains 0.4-0.6 bits/genome, e.g. UBA9145
  module_5 0.97 → 0.38, module_83/84 in spot_13, module_64/65 with φ = 1.0; MAG-120531 module_6, the
  tetrapyrrole/SUF modules 2, 4, 5 in spot_33). 23/69 and 10/40 modules have a partner with |φ| ≥ 0.8:
  PPanGGOLiN often splits one island into several modules. **No module pair is mutually exclusive
  (φ ≤ −0.5: none)** at this sample size. Per-module table: `module_gains.tsv`.
- **The context model sharpens, it does not flatten.** The mean drops, but the tail grows: the largest LOO
  surprise is 8.7 bits (prevalence: 4.7), and 66/105 UBA9145 states above 3 bits stay above 3 bits.
  Mixing 10% prevalence into the context prediction caps the maximum at 6.4-7.5 bits for +0.004-0.009 bits
  mean. Use the mixture for reporting.
- **Chesterfield bins: the 36 Chesterfield-only UBA9145 modules are exactly as surprising under the context
  model (5.75 bits each, 207 bits per bin), because no public genome carries them.** A model trained only
  on public genomes cannot predict content it has never seen; it falls back to the prior. What the context
  model does show is that the 36 modules are **one event, not 36**: given the other bin, a context model
  predicts each of them at 0.29 bits (logistic on modules). The correct total is about one block surprise
  (≈ 5.8 bits) plus small per-module terms, not 207 bits. The second bin is the same contigs from an
  overlapping co-assembly (see coverage_test.md), so it is not independent support. The block also
  looks foreign: 153 of its 208 genes in chesterfield_maxbin_040 sit on 27 contigs with no UBA9145
  persistent gene. Both bins carry 10-11% CheckM2 contamination, and coverage_test.md assigns 21/36 modules
  to contamination.
- **The bins' other module states:** with public training, context explains 25-40% of the remaining surprise
  in the UBA9145 bins (39.5 → 23.7 and 40.1 → 31.0 bits, logistic on modules). It turns the shared
  module_7/17/18 trio (3/52 public carriers, co-occurring) from 3.75 to 2.15 bits each. It makes absence of
  module_61 (39/52) more surprising (1.95 → 3.4-3.8 bits): the rest of the genome says module_61 should be
  there. For MAG-120531, the 30%-complete lorbin_024 drops from 27.5 bits (plain prevalence) to 15.0
  (context), the same as completeness-aware prevalence (16.0). The model learns incompleteness from
  missing context and is not fooled by it.

## Overfitting and data size

- With 51-66 training genomes and 940-1,730 features, every model fits the training set almost perfectly
  (in-sample 0.01-0.16 bits vs LOO 0.37-0.58). Only strong shrinkage helps: C = 0.1 beats C = 1 for
  family features. Weak shrinkage (C = 1) makes 12-28% of targets worse than prevalence.
- Learning curve (module targets, LOO, training sets subsampled, 3 repeats), gain of the family logistic over
  prevalence trained on the same genomes:

  | Training genomes | 10 | 25 | 40 | 51 / 66 (all) |
  |---|---|---|---|---|
  | UBA9145 gain, bits (% of prevalence loss) | 0.098 (16%) | 0.155 (26%) | 0.173 (29%) | 0.189 (32%) |
  | MAG-120531 gain | 0.064 (11%) | 0.106 (19%) | 0.129 (24%) | 0.155 (29%) |

  The gain is still rising at full size, about +0.03-0.04 bits per doubling, with no plateau. At about
  10 genomes the context model already gains something but is unreliable per module. A full context model
  (pairwise exclusions, island structure, clade terms) needs hundreds of genomes per species. With
  ~50 genomes, positive co-occurrence of island fragments is the only signal that is learnable.
  Mutually exclusive alternatives will need either more genomes or pooling across species (shared
  module/family identities, e.g. the same model for all genera in a family).
- Leakage caveat: family targets mask only their own module, so operon neighbours outside any module can
  still predict them. Module-target numbers are clean of this (all of the module's families masked).

## Recommendations

1. Level 2 as an L2 logistic on accessory families plus completeness, C ≈ 0.1, mixed 9:1 with
   level-1 prevalence. Report per-module LOO gain, so that surprise is only "context-explained" where the
   model has earned it.
2. Score co-occurring novel modules as blocks (chain rule over a block ordering, or merge modules with
   identical presence vectors) before summing a bin's surprise. Otherwise one contaminating or HGT island
   counts 36 times.
3. A block of never-seen modules on contigs without core genes is a contamination flag, not biology. Check
   it with coverage (coverage_test.md) before interpreting.

## Files

Code: `experiments/context_model/context_model.py` (LOO predictions, in-sample fit, Chesterfield
predictions), `learning_curve.py` (reduced curve: modules only, sizes 10/25/40/all), `analyse.py`
(tables). Outputs: `/data/scratch/metapanclust_data/context_model/{uba9145,mag120531}/`
(`loo_predictions.tsv`, `insample_predictions.tsv`, `chesterfield_surprise.tsv`, `module_gains.tsv`,
`learning_curve_summary.tsv`, `novel_contigs.txt`) and `analysis.txt`. Runtime about 13 CPU-min per species
for the main run, 3 for the curve, on 1 thread each. Module annotations come from Bakta for
chesterfield_maxbin_040 / chesterfield_semibin_019 only, so modules found only in public genomes are
unannotated.

## For grex (not run here)

A full learning curve (more repeats, family targets, all C/γ) and a version with a core-genome
phylogeny feature are ~2-4 CPU-h per species. Run `context_model.py` with the learning-curve loop restored
(per-target parallel), with no new inputs beyond the PPanGGOLiN outputs. Pooling across many GTDB
species is the real test of data size and would need PPanGGOLiN runs per species first.
