# Scorer v0: the surprise model as one package, and a held-out test (2026-11-21)

Code: `metapanclust/` (`python -m metapanclust ref | score | eval`). Held-out test: `experiments/scorer_eval/`
(grex job 8195129, 9 species, 98 held-out genomes, 1,800 simulated bins). Full tables:
[scorer_eval_results.txt](scorer_eval_results.txt).

## What the package does

- **Reference** (`ref`): one species pangenome from a PanGBank/PPanGGOLiN `.h5`. It holds family presence per
  genome, completeness per genome, and scoring units: each PPanGGOLiN module is one unit scored as a block, and
  every other family is its own unit. It also holds each unit's median span in bp, and every distinct protein
  labelled with its family (an MMseqs2 database).
- **Annotation**: Pyrodigal genes, then MMseqs2 against the reference proteins (≥ 80% identity and coverage, the
  PanGBank clustering cut). Genes without a hit are novel to the species.
- **Level 1**: two-state gain/loss chain on a UPGMA gene-content tree (Jaccard on accessory families) with the
  query as an uncertain tip. Rate and stationary frequency are averaged over a grid. Units are split into 5
  groups; each group's prediction uses a query placement that ignores that group's families. A Platt
  recalibration is fitted on the reference genomes' own leave-one-out posteriors.
  ([tree_vs_context.md](tree_vs_context.md): gene content is the best within-species relatedness.)
- **Level 2**: masked L2 logistic per unit on the other variable accessory families plus completeness. Its
  prediction is averaged 50/50 with level 1.
- **Detection**: completeness is estimated from the detected fraction of persistent families, or CheckM2 is used.
  Units longer than a gene also pay the fragmentation term, relative to a single gene's span, so contig-break
  losses are not counted twice.
- **Output**: per unit p, detection, q, bits and calibrated excess (presence direction where q < 1/2). Per bin:
  Z+, Z−, rare-present and expected-absent counts. Novel genes are scored as a negative-binomial tail
  against the reference's own leave-one-out singleton counts (singletons concentrate in a few genomes, so a
  Poisson is far too tight). Per contig: novel and persistent gene counts, presence excess and median identity.

## Held-out test

Nine genera from nine phyla (Akkermansia, Desulfovibrio, Flavobacterium, Fusobacterium, Helicobacter_D,
Methanosarcina, Sphingomonas, Stutzerimonas, Zag111). For each, the species with the most genomes is the
reference and a congener supplies contamination. Up to 20 complete genomes (CheckM2 ≥ 95%) are held out each.
Each is cut into lognormal pieces (median 15 kb), reduced to 70% or 90% completeness, and either kept clean,
contaminated with 5/10/20% congener pieces, or merged with 10% of another genome of the same species.

**Unit prediction** (clean bins): the model cuts bits per unit by a median 9% against species prevalence (0.4-27%
by species and completeness; more at 90% completeness than at 70%). Context adds 1.2% over the tree alone.
This is smaller than in tree_vs_context (45-56%) because every unit is scored here, including near-universal
and singleton ones, and because the query is fragmented and incomplete.

**Calibration is off in a consistent direction.** Units predicted at 1-5% are seen 0.8% of the time, and those
at 0.1-1% are seen 0.13%. Middle predictions are too low (0.36 predicted against 0.50 seen). Z+ on clean bins
has a median of −4.5, where 0 would mean calibrated. Platt recalibration on the reference genomes barely
changes anything (a = 0.23, b = 1.10), so the tree is roughly calibrated for complete reference genomes. The
likely cause is placement: a fragmented 70-90% query sits less precisely on the tree, which pulls rare units
towards their prevalence. A fix would recalibrate on references degraded to the query's completeness. Until
then, absolute Z values should not be read as standard normal.

**Contamination by a congener shows up as novel genes, not as surprising known units.**

| Score | foreign 5% vs clean | foreign 20% vs clean | strain 10% vs clean |
|---|---|---|---|
| novel-gene bits | 0.92 | 0.99 | 0.49 |
| Z+ | 0.53 | 0.60-0.62 | 0.59-0.63 |
| rare-present units | 0.55-0.57 | 0.64-0.66 | 0.62-0.65 |

(AUC, both completeness levels.) Novel bits also separate a foreign merge from a strain merge (AUC 0.97), and they
should: a congener brings genes the species has never seen, while a strain brings known accessory genes.
Per species the 5% AUC ranges from 0.74 (Stutzerimonas) to 1.0 (Flavobacterium); species with large open
pangenomes have high novelty baselines. Strain merges are only weakly visible (AUC about 0.6 from rare
units).

**Contig level:** novel-gene fraction finds the foreign contigs in contaminated bins at AUC 0.91-0.93. Unit
presence excess does not (0.38-0.46, because foreign genes mostly do not map at all). Identity to the reference
is not simulated here (native genes keep their family labels), so it is tested only on CAMI.

**Completeness estimate** from persistent families: 0.68 for a 70% target and 0.86 for a 90% target. Genes cut by
a contig break are lost in the simulation, so the true gene-level completeness is a little below the target.

## What this says about the model

- The useful signal for bin quality is **novelty against the species pangenome**: genes and contigs with no
  species hit. It works without a tree, a context model or CheckM2. It is what merge_identity found from the
  other side (per-gene divergence to the focal pangenome, AUC 0.99 at 94-95% ANI).
- **Level-1/level-2 unit surprise is not a contamination detector.** It describes the bin's accessory content
  relative to its relatives. That is the biological question we wanted, but it needs the calibration fix before
  "surprising" can be stated with a number.
- The CAMI II marine test ([cami_marine.md](cami_marine.md), running) checks both claims on real binner output
  against CheckM2.
