# Tree vs context: relatedness or gene-gene context? (design, local test, grex submission, 2026-11-19)

Status: scripts pass a local end-to-end test on 3 species (plus a timing run on a 318-genome species). The
grex run has not been submitted. The numbers below come from the small test only. They are not results.

## Question

The context model (pangbank_scale part A) removes 27-70% of the prevalence log loss, but its inputs (accessory
families) also track clade structure, and part A's only relatedness baseline was a Jaccard kernel on the same
families (overnight_grex.md, caveat of 2026-11-19). This run separates the two:

1. How much of the context gain is relatedness (phylogeny), and how much is genuine gene-gene context?
2. Which relatedness measure is right? Sequence trees and gene-content relatedness can disagree in both
   directions (e.g. coral-symbiont Marine Group I archaea: GC-shifted outliers on a sequence tree, close by
   gene content). No single tree is treated as truth: three measures that fail differently are compared, and
   genomes where they disagree are listed.

## Design

**Species.** `select_species.py`: 200 of part A's 291 species (stratified by genome count 25-49, 50-99,
100-199, 200-399, ≥ 400; 40 each, round-robin over phyla; 33 phyla; seed 1). Species with < 25 genomes are
dropped. The .h5 files are already on grex. Per species a pool of at most 400 genomes is drawn (seed 1), and
everything is computed on that pool: 38,287 genomes in total.

**No genome download is needed.** The PanGBank .h5 files store every gene's DNA (`Genes_with_Sequences: true`,
deduplicated table `/annotations/sequences`) plus contig order, and PanGBank metadata includes genome GC
(`/metadata/genomes/pangbank_wf/gc_percentage`). `worker.py prep` reads them directly with PyTables (no
`ppanggolin write_pangenome`).

**Relatedness measures (on the pool).**

| Measure | How | Fails when |
|---|---|---|
| `ani` | skani 0.3.1 `triangle -E --min-af 0` on each genome's CDS concatenated per contig in gene order (≈ 88% of the genome; ANI on coding sequence) | composition shifts (GC/codon usage), near-clonal MAGs (ANI noise ~ signal) |
| `aai` | core-protein distance: Mash-style distance from Jaccard of FracMinHash (1/32) amino-acid 5-mers of persistent families, k-mers keyed by family (estimates 1 - AAI over the core; synonymous and GC-driven changes invisible) | very recent divergence (few amino-acid changes) |
| `jac` | 1 - Jaccard on accessory (shell + cloud) families, with every target family (and its module's families) of the subsample removed | HGT, assembly/annotation incompleteness; shares inputs with the context model |

**Predictors** of held-out presence (10-fold CV; same folds for all; targets as in part A: up to 60 modules and
60 accessory families per subsample; N = 25, 50, 100, 200, 400, 3 subsamples each, 1 for N = pool):

- `prev`: Beta(1,1) prevalence.
- `ker_{ani,aai,jac}`: relatedness kernel w = exp(-d / (q · median d)), bandwidth q ∈ {0.03, 0.1, 0.3, 1, 3}
  chosen per target and fold by leave-one-out on the training genomes.
- `tree_{ani,aai,jac}`: UPGMA tree from the distance; two-state gain/loss chain (phylo_prior.md §1) with
  stationary root; per target a Bayesian average over a grid of turnover rate μ·H ∈ {0.1 … 1000} (H = root
  height; 1000 ≈ no tree signal) and π ∈ {0.05 … 0.95}, weighted by the pruning marginal likelihood of the
  training tips. Test tips are uncertain tips [1, 1]; an up-down pass gives their posterior (and leave-one-out
  posteriors for training tips). Tip model: P(observed | present) = completeness (persistent-family fraction,
  capped 0.99), P(observed | absent) = 0.01.
- `relseq`, `relall`: stacked relatedness, logistic (C = 1) on the logits of ker/tree × (ani, aai) or
  × (ani, aai, jac), trained on the training genomes' leave-one-out relatedness predictions.
- `lr`: context only, part A's model (L2 logistic, C = 0.1, accessory families with the target's module masked,
  plus completeness).
- `lr_relseq`, `lr_relall`: context + relatedness: `lr` plus the top min(10, N/10) principal coordinates (PCoA) of
  each relatedness distance as extra covariates, penalised like the families. They carry no outcome
  information, so training and test rows are treated alike.
- `avg_relseq`, `avg_relall`: untuned 50/50 average of `lr` and `relseq` / `relall`.

Two ways to combine context with relatedness through the outcome failed in testing and were dropped: a fixed
offset from leave-one-out relatedness predictions, and a free coefficient on cross-fitted ones. Leave-one-out
predictions for rare targets are anticorrelated with the outcome at small N (a carrier's own presence is removed
from its prediction), so the training and test offsets come from different distributions. The combined model
then lost up to 2× against `lr` alone (0.71 vs 0.34 bits).

**The answer to question 1** is the decomposition in `aggregate.py` §2, computed per species and then as a median:

- ctx = % of the prevalence loss removed by `lr`;
- rel = the same for `relseq` (or `relall`);
- both = the same for `lr_rel*`;
- unique context = both - rel;
- shared = ctx + rel - both.

Because PCoA captures only coarse structure, `both` and hence unique context are upper bounds on pure context. The
`avg_*` rows give a second, untuned combination. Use `relseq` for "context beyond sequence phylogeny". `relall`
also includes gene-content relatedness, which partly is context.

**Disagreement diagnostics** (pool level, `ID.diag.tsv`, `ID.matcor.tsv`):

- Spearman correlation between the distance matrices.
- For each genome:
  - the percentile of its nearest-neighbour distance and the z-score of its mean distance under each measure;
  - the overlap of its 5 nearest neighbours between measures;
  - its GC content (PanGBank metadata) and the deviation from the species median.
- Flags:
  - `flag_ani_out_jac_in`: outlying under ANI (mean-distance z > 2 or nearest-neighbour distance in the top 5%),
    typical under gene content (nearest neighbour ≤ median, z < 1). This is the Marine Group I pattern.
  - The reverse flag, and the same pair of flags for AAI.
- For flagged genomes, `aggregate.py` compares mean held-out bits under the best sequence-based relatedness and
  the best content-based relatedness. This shows which measure predicts them. The script also reports the
  Spearman correlation, over all genomes, of the content-vs-sequence advantage with |GC deviation|.

## Local test

Docker `metapanclust-tools:0.1`, 1 thread. skani came from the danaSeq conda env `dana-mag-derep`, 2-3
threads. Three species ran through `run_species.sh` (RUNMODE=local):

- Deinococcus wulumuqiensis (3022): 20 genomes, scored at N = 20 (shown in the 25 column);
- Geothrix sp903857495 (3442): 50 MAGs;
- Cutibacterium avidum (2981): 150 genomes; the N = 150 run is shown in the 100 column.

Raw output: `/data/scratch/metapanclust_data/tree_vs_context/test2/` (`aggregate.txt`).

**% of the prevalence loss removed** (median over the species present at each N; 3 / 2 / 1 species):

| Predictor | modules 25 | modules 50 | modules 100-150 | families 25 | families 50 | families 100-150 |
|---|---|---|---|---|---|---|
| ker_ani | 30 | 22 | 50 | 11 | 19 | 53 |
| tree_ani | 50 | 32 | 66 | 14 | 28 | 69 |
| ker_aai | 28 | 33 | 50 | 17 | 29 | 54 |
| tree_aai | 50 | 45 | 63 | 20 | 39 | 67 |
| ker_jac | 31 | 27 | 42 | 21 | 26 | 44 |
| tree_jac | 46 | 44 | 68 | 33 | 44 | 71 |
| relseq (ANI + AAI stack) | 50 | 50 | 67 | 23 | 42 | 70 |
| relall (+ gene content) | 51 | 58 | 69 | 42 | 55 | 72 |
| **lr (context only)** | **57** | **61** | **75** | **44** | **57** | **76** |
| lr + PCoA(seq) | 57 | 62 | 75 | 45 | 58 | 76 |
| lr + PCoA(all) | 57 | 62 | 75 | 46 | 58 | 76 |
| avg(lr, relall) | 57 | 63 | 75 | 48 | 59 | 76 |

Indicative findings (3 species):

- **The tree predictor beats the kernel for every measure.** The gain/loss chain with per-target rate
  averaging adds 10-20 points over the tuned kernel on the same distances. Part A's Jaccard-kNN baseline
  therefore understated relatedness badly: in part A it reached about half of the context gain, here a
  gene-content tree reaches 80-90% of it.
- **Sequence relatedness alone (ANI + AAI) already removes most of what context removes.** Modules: 50 vs 57,
  50 vs 61, 67 vs 75. Families at small N: 23 vs 44 and 42 vs 57. What context adds on top of sequence
  relatedness is unique context = 6-20 points; on top of all relatedness it is 2.5-6.5 points
  (decomposition, `aggregate.txt` §2).
- **Adding relatedness to the context model barely helps** (lr + PCoA is within 1.4 points of lr). The context
  model has already learned the clade structure from gene content, so in this test the overlap ("shared") is
  most of the gain.
- **Which measure wins:**
  - for family targets, the gene-content tree wins in every case (by 6-16% of bits over the best sequence
    tree);
  - for module targets, it splits between tree_aai and tree_jac;
  - ANI never wins outright. On the near-clonal Geothrix MAGs (median ANI on CDS 99.91%, AAI distance 1e-4),
    ANI carries no signal (ker_ani ≈ prevalence). Before the tip-branch floor (see below), the ANI tree did worse
    than prevalence (1.26 vs 0.89 bits).
- **Disagreement:**
  - Spearman correlation between distance matrices (interquartile range over the 3 species): ANI-AAI 0.58-0.80, AAI-Jaccard 0.43-0.85, ANI-Jaccard
    0.29-0.77.
  - |GC deviation| correlates weakly with sequence outlyingness (ρ ≈ 0.24 for both ANI and AAI).
  - 4 genomes were flagged, all in Geothrix. GCA_903857495.1 is outlying under ANI (nearest neighbour in the
    top 5%) but typical by gene content, and it has the largest GC deviation in the species (-0.13 points). It is
    predicted best by gene content (0.81 bits) and AAI (0.89), worst by ANI (1.12 = prevalence). This is the
    direction of the Marine Group I example, though with a tiny GC shift.
  - The three genomes that are outlying by gene content are predicted better by AAI than by gene content.

Bugs found and fixed during the test (the reasons behind the design choices above):

- zero-distance pairs (identical CDS) forced identical content at every rate. A tip-branch floor of 0.01 × H
  lets the rate grid express "relatedness does not inform this target";
- the rate grid now goes up to μH = 1000;
- the tip error of 0.01 bounds a single overconfident miss at about 6.6 bits;
- combining through leave-one-out offsets was dropped (see Design).

Caveats:

- 3 species, two of them small. The table is a smoke test, not a result.
- The AAI estimator is Jaccard-based. A genome that lacks a persistent family looks more distant. Module targets
  that contain a family PPanGGOLiN calls persistent can therefore leak slightly into `aai`; family targets
  (shell or cloud only) cannot. ANI (skani, aligned fraction ignored) and `jac` (targets removed) do not leak.
- ANI is on concatenated CDS, not whole genomes (no NCBI download). This is close to whole-genome ANI for
  bacteria, but it is not identical.
- Family targets keep part A's operon-leakage caveat for `lr`.
- UPGMA assumes a clock. NJ and lineage rate multipliers are not implemented.

## Runtime and memory (grex estimate)

Measured locally, 1 thread (genoa cores are similar or faster):

| Step | 150 genomes (2981) | 318 genomes (3530) | 400-genome pool (estimate) |
|---|---|---|---|
| prep (h5 read, CDS fasta, AAI) | 19 s, 0.86 GB | 80 s, 1.2 GB | 2-5 min, ≤ 4 GB for the 1,500-genome .h5 files (not tested) |
| skani triangle | 45 s (3 threads) | 165 s (2 threads) | 4-6 min (2 threads; slower under a full node) |
| model, per subsample | N = 25: 28 s; 50: 37 s; 100: 50 s; 150: 63 s | not run | N = 200: ~90 s, N = 400: ~3-4 min |
| model, whole species | 7 min (10 subsamples), RSS 0.2 GB | | ~15-20 min (13 subsamples), RSS < 1 GB |

- **Per species:** ≤ 30 min (1 thread plus 2 for skani), about 100 CPU-h for 200 species.
- **Layout:** 2 nodes × 100 species in parallel (PAR = 100), one wave. Expected wall time is 40-75 min
  against the 3-h limit, with a per-species timeout of 9,900 s. The model writes its outputs after every
  subsample, so a timeout loses only the largest N.
- **Memory:** ≤ 100 × 1-4 GB per node, well under `--mem=600G`.
- **Disk:**
  - transient CDS fasta: 0.75 MB per genome, ≤ 30 GB and about 38k files in total; deleted per species right
    after skani;
  - `work/<id>` is removed when the species finishes;
  - final outputs: about 1-3 MB per species.

## Exact submission (coordinator)

The coordinator copies the files over and submits; `grex` below stands for the usual route (via concentration).
All paths are resolved (`/home/rec3141/scratch` is a symlink), so first resolve it on grex:

```bash
# 0. on grex: resolved scratch path (use the printed value for $S below)
readlink -f /home/rec3141/scratch

# 1. from the workstation (S = resolved path from step 0)
S=<resolved scratch>
ssh grex "mkdir -p $S/metapanclust/tree_vs_context/{code,sel,out,logs}"
scp /data/dev/metapanclust/experiments/tree_vs_context/{select_species.py,worker.py,run_species.sh,run_grex.sbatch,aggregate.py} \
    grex:$S/metapanclust/tree_vs_context/code/
scp /data/scratch/metapanclust_data/tree_vs_context/sel/{species_tvc.tsv,tasks_tvc.tsv} \
    grex:$S/metapanclust/tree_vs_context/sel/

# 2. on the grex login node: check inputs (200 .h5 present; both containers exist)
S=$(readlink -f /home/rec3141/scratch); cd $S/metapanclust/tree_vs_context
for i in $(cat sel/tasks_tvc.tsv); do [ -s $S/metapanclust/pangbank_scale/h5/$i.h5 ] || echo missing $i; done
ls -l $S/metapanclust/sif/metapanclust-tools-0.1.sif $(readlink -f /home/rec3141/GENICE/danaSeq/mag_analysis/.danaseq-mag-analysis.sif)

# 3. submit: 2 whole genoa/genlm nodes, 100 species each in one wave
sbatch --array=0-1 code/run_grex.sbatch
# node logs: logs/tvc_<job>_<i>.out; per species: logs/tasks/<id>.log
# finished species: out/res/<id>.done (resubmitting the same command skips them and reuses finished
# prep/skani steps in out/work/<id>/)

# 4. aggregate (login node, < 1 min)
export PATH=/opt/singularity-ce-4.4.1/bin:$PATH
singularity exec -B $S $S/metapanclust/sif/metapanclust-tools-0.1.sif \
    python code/aggregate.py out sel/species_tvc.tsv > out/aggregate.txt
```

Outputs, per species, in `out/res/`:

- `<id>.summary.tsv`: per N, subsample and kind, the mean bits of every predictor;
- `<id>.targets.tsv.gz`: per target;
- `<id>.genomes.tsv.gz`: per genome and subsample;
- `<id>.diag.tsv`: per-genome disagreement diagnostics;
- `<id>.matcor.tsv`: correlations between the distance matrices.

The aggregate in `out/agg/` has `per_species.tsv`, `pct_{module,family}.tsv`, `decomposition.tsv`,
`winner.tsv`, `matcor.tsv`, `genome_diag.tsv.gz`, `discordant_summary.tsv` and `discordant_genomes.tsv`. The
markdown tables go to `out/aggregate.txt`. `out/work/<id>/` is removed when a species finishes.

To regenerate the selection on grex: `singularity exec -B $S $SIF python code/select_species.py
$S/metapanclust/pangbank_scale/sel/species_A.tsv sel`. It is deterministic (seed 1) and gives the same 200 species.

## Files

- Code (`experiments/tree_vs_context/`):
  - `select_species.py`;
  - `worker.py` (prep / model);
  - `run_species.sh` (prep, then skani, then model; RUNMODE=grex uses singularity, RUNMODE=local uses docker
    plus the danaSeq conda env);
  - `run_grex.sbatch`;
  - `aggregate.py`.
- Local data (`/data/scratch/metapanclust_data/tree_vs_context/`):
  - `sel/`: the selection;
  - `test2/`: the end-to-end test (`out/res`, `out/agg`, `logs`, `aggregate.txt`);
  - `test3/`: prep and skani timing for the 318-genome species;
  - `test/`, `probe/`: earlier iterations.

## Results: grex run (job 8190397, 2 nodes, 200 species / 33 phyla, 2026-11-19)

Full tables: `tree_vs_context_results.txt`; per-genome discordance list on grex and in
`/data/scratch/metapanclust_data/tree_vs_context/grex_results/out/`.

**1. Gene content is the best relatedness measure.** The gene-content (Jaccard) tree is the single best relatedness
predictor for 62-65% of species at N = 25, rising to 97-98% at N >= 200; it beats the best sequence-based predictor
(ANI or core-protein distance) by a median 1% (N = 25) to 11% (N = 400) of prevalence loss. Trees beat kernels on
every measure by 9-25 points, so the earlier Jaccard-kNN baseline understated relatedness by about half.

**2. Genuine context gain is real and grows with N** (median % of prevalence loss; family / module targets):

| N | context (lr) | all relatedness | unique to context beyond all relatedness | beyond sequence relatedness only |
|---|---|---|---|---|
| 25 | 26.4 / 33.1 | 18.5 / 23.0 | 5.5 / 8.5 | 7.7 / 10.2 |
| 100 | 51.4 / 53.5 | 45.1 / 42.8 | 6.6 / ~10 | 10.4 / 14.9 |
| 400 | 70.0 / 67.8 | 56.4 / 55.9 | 11.1 / ~12 | 16.3 / 18.5 |

Most of the context model's gain is shared with relatedness, but 83-100% of species gain more than 1 point from
context beyond all relatedness, and that unique part grows with N. Feeding relatedness into the context model
adds nothing (lr_relall = lr): it already learns clade structure from gene content.

**3. The measures disagree substantially.** Median Spearman between pairwise distances: AAI-Jaccard 0.67,
ANI-AAI 0.51, ANI-Jaccard 0.50. A genome's GC deviation from its species median tracks how outlying it is in gene
content (median rho 0.34) more than in ANI (0.12) or AAI (0.25).

**4. Discordant genomes (937 of 38,287).**
- Outlier by gene content but typical by ANI (502 genomes, 111 species): median |GC deviation| 0.57 vs 0.11 for the
  rest. These are mostly MAGs (e.g. gut CAG-873, Hominicoprocola); the likeliest explanation is contamination or
  chimerism adding foreign genes and shifting GC, not biology.
- Outlier by ANI but typical by gene content (250 genomes, 79 species): GC deviation like the rest (0.13 vs 0.11).
- In every discordant class the gene-content tree predicts the genome's content better than the sequence trees
  (e.g. 0.21 vs 0.27 bits for ANI-outlier/gene-content-typical genomes).

**Scope note on the Marine Group I example:** this run compares genomes *within* species, so it cannot test the
cross-genus case (GC-shifted coral symbiont MGI placed as phylogenetic outliers while gene content shows close
relationship). Within species, GC deviation goes with gene-content outlierness, mostly from MAG quality. The MGI
question needs a between-species version (genus- or family-level trees vs gene-content similarity across species).
