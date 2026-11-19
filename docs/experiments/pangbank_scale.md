# PanGBank data-size experiment: design, local test, grex submission (2026-11-18)

Status: scripts tested locally on 4 species (A) and 2 genera (B). The full run on grex has not been
submitted yet. The numbers below come from the small test only. They are not results.

## Questions

- **A (level 2).** How much does the context model (`experiments/context_model/context_model.py`) gain over
  Beta(1,1) prevalence as the number of genomes per species grows? The context model is an L2 logistic
  regression on accessory-family presence plus completeness. The design uses a few hundred GTDB species
  with 15 to 1,500 genomes each, subsampled to N = 15, 25, 50, 100, 200, 400 and 800.
- **B (level 1).** Does shrinking a species' cluster prevalence towards its genus help predict held-out
  genomes, especially when only 5 to 20 training genomes are available?

## Design

Data: PanGBank GTDB_all 1.0.0 (GTDB r226, 4,351 pangenomes, PPanGGOLiN 2.2.4 `.h5`). The `.h5` files
open in the container's PPanGGOLiN 2.2.1.

**Selection (`select_species.py`, seed 1).**

- A uses 291 species from 52 phyla. Species are stratified by genome count into 15-29, 30-59, 60-119,
  120-249, 250-499, 500-999 and 1000-1500 genomes. The selection takes 45 species per stratum (only 21 exist
  in 1000-1500), round-robin over phyla. Species above 1,500 genomes are excluded to keep downloads small.
- B uses 40 genera with at least 3 PanGBank species: 203 species from 23 phyla. Each genus contributes up to
  8 species of at most 300 genomes, smallest first.
- Downloads: 455 `.h5` files, 91,006 genomes, about 40-46 GB. Test files ran at 0.25-0.45 MB per genome,
  and the 15 test files took 47 s.

**A worker (`worker.py species ID`).**

- Inputs: `ppanggolin write_pangenome --Rtab --modules --stats --partitions`. Completeness is
  PPanGGOLiN's own value from `genomes_statistics.tsv`.
- Subsamples: for each N, 3 random genome subsamples (1 at N = 800, and 1 using all genomes if
  n ≤ 400).
- Within a subsample:
  - features are shell and cloud families with between 2 and N−1 carriers;
  - targets are up to 60 modules (completion ≥ 0.5, 2 to N−2 carriers) and up to 60 accessory families
    (max(2, N/20) to N−max(2, N/20) carriers);
  - the target's own module families are masked;
  - each target is scored by 10-fold cross-validation, or leave-one-out when N < 10.
- Models, all trained on the same folds:
  - Beta(1,1) prevalence;
  - Jaccard-weighted neighbours with γ = 2 and 8 (clade proxy, as in context_model.py);
  - logistic regression with C = 0.03, 0.1, 0.3 and 1;
  - a 0.9 / 0.1 mix of logistic (C = 0.1) and prevalence.
- Output: mean bits per (genome, target).
- Families inherit the leakage caveat of context_model.md: operon neighbours outside a module are not
  masked. Module targets are clean.

**B worker (`worker.py genus GENUS ID,...`).**

- Genus clusters:
  - each species is subsampled to at most 200 genomes;
  - family representative proteins come from `ppanggolin fasta --prot_families all`;
  - all representatives are pooled and clustered with `mmseqs easy-cluster --min-seq-id 0.5 -c 0.8`.
    50% identity is the setting that recovered the core across congeners in pangbank.md;
  - a genus cluster is present in a genome if any of the species' families in it is present.
- For each focal species, training size m ∈ {5, 10, 20} and 5 random splits:
  - training uses m genomes, testing uses up to 50 others;
  - the score is held-out bits per test genome, summed over all genus clusters.
- Predictors:
  - `sp`: Beta(1,1);
  - `sp_eb`: one Beta(a, b) for all clusters, fitted on the m training genomes;
  - `genus`: leave-focal-out genus mean, each congener capped at 20 genomes as in model_spec §3.4 A;
  - `shrunk`: (x + κ p_genus) / (m + κ), where κ is the empirical-Bayes beta-binomial maximum (grid
    0.1-2000) fitted on the congeners only, each against its own leave-out genus mean, as in model_spec
    §3.5;
  - `mix`: 0.5 `shrunk` + 0.5 `sp_eb`;
  - `oracle`: the best κ on the test genomes (upper bound for κ).
- Strata:
  - by the focal species' full-data prevalence: core, accessory, rare_or_absent;
  - by where a cluster occurs: shared, focal_only, congener_only.
- The taxonomic family level (several genera) is not implemented yet. Its natural extension is to pool the
  genus clusters of the genera in a family.

## Local test (docker `metapanclust-tools:0.1`, 1 thread per worker)

### A: species tested

| Species (id) | Phylum | Genomes | Wall time | Peak RSS |
|---|---|---|---|---|
| Deinococcus wulumuqiensis (3022) | Deinococcota | 20 | 40 s | 0.4 GB (ppanggolin) |
| Geothrix sp903857495 (3442) | Acidobacteriota | 50 | 125 s | 0.4 GB |
| Cutibacterium avidum (2981) | Actinomycetota | 150 | 704 s | 0.7 GB |
| Helicobacter pylori_CN (3530) | Campylobacterota | 318 | 333 s | 0.9 GB |

- 2981 and 3442 ran the earlier fold scheme, which used leave-one-out up to N = 25 (50-90 s per subsample
  at N = 25). The final code uses 10 folds throughout.
- Each (N, subsample) takes 10-65 s at N ≤ 318, so a full species is about 5-25 min on 1 thread.
- N = 800 has not been timed. A rough extrapolation is 10-25 min for its single subsample. Results are
  written after each subsample, so a timeout loses only the largest sizes.

### A: gain over prevalence (4 species, so indicative only)

Median % of prevalence log loss removed:

| N | Modules: neighbours γ = 2 | Modules: logistic C = 0.1 | Modules: logistic C = 0.3 | Families: logistic C = 0.1 | Families: logistic C = 0.3 |
|---|---|---|---|---|---|
| 15 | 32 | 44 | 51 | 45 | 49 |
| 25 | 29 | 42 | 48 | 42 | 45 |
| 50 | 24 | 50 | 52 | 41 | 40 |
| 100 | 26 | 60 | 62 | 51 | 50 |
| 200 | 11 | 60 | 63 | 40 | 40 |

With all genomes, the logistic (C = 0.1) removes 48-79% of the prevalence loss. The neighbour model removes
8-56%.

These gains are larger than the 30% found for the two Chesterfield MAG species, for two reasons:

- PanGBank species are isolate-rich and clade-structured, so much of the signal is clade identity. The
  neighbour baseline captures part of it.
- Every test species has strongly co-inherited modules.

The full run separates the context model from clade structure: compare `pct_lr` against `pct_knn` per
species, and examine trends by N. The best C drifts from 0.3-1 at small N towards 0.1 for families at
N ≥ 50. The aggregate script reports the best fixed C per N, not a per-species choice.

### B: genera tested

| Genus | Species | Genomes used | Genus clusters | Wall time | Peak RSS |
|---|---|---|---|---|---|
| Marinisoma | 3 | 49 | 3,111 (1,048 in ≥ 2 species) | 24 s | 0.9 GB (MMseqs2, 2 threads) |
| Flavobacterium | 8 | 437 (287 → 200) | see below | ≈ 6 min | 1.2 GB |

Held-out bits per test genome, summed over the genus clusters in each stratum. Values are means over 11
focal species; at m = 20 only 4 focal species are large enough. "Wins" is the fraction of focal species
where `shrunk` beats `sp_eb`.

| Stratum | m | sp Beta(1,1) | sp_eb | genus | shrunk (EB κ) | mix | oracle κ | shrunk wins |
|---|---|---|---|---|---|---|---|---|
| all | 5 | 5,752 | **2,156** | 8,524 | 2,516 | 2,239 | 2,483 | 0/11 |
| all | 10 | 3,877 | **1,918** | 8,487 | 2,163 | 1,976 | 2,140 | 0/11 |
| all | 20 | 3,096 | **1,627** | 9,861 | 1,762 | 1,665 | 1,751 | 0/4 |
| shared (focal + congener) | 5 | 1,038 | 911 | 2,727 | 944 | **901** | 926 | 3/11 |
| shared | 10 | 874 | **820** | 2,724 | 851 | 820 | 836 | 3/11 |
| focal_only | 5 | **1,049** | 1,152 | 3,776 | 1,405 | 1,208 | 1,335 | 0/11 |
| congener_only | 5 | 3,664 | **94** | 2,021 | 167 | 130 | 33 | 1/11 |

Bold marks the best predictor that does not look at the test genomes. The oracle κ is not a valid
predictor; it bounds what a better κ could achieve.

**Test result: genus shrinkage, as specified, loses to species-only empirical Bayes in every stratum.**

- The fair baseline is not Beta(1,1). Most genus clusters are absent from any one species, so Beta(1,1)
  pays about 0.2 bits for each of thousands of absent clusters. Against Beta(1,1), shrinkage "wins" by
  55-60%. That is an artefact of the cluster universe. A single species-wide Beta(a, b) fitted on the m
  training genomes (`sp_eb`) removes it.
- The fitted κ is small (median 0.5, range 0.5-2), meaning species prevalences of a cluster are U-shaped
  around the genus mean. Most clusters behave as species markers rather than genus-wide shell.
- Even the best κ on the test set loses to `sp_eb` overall, so a better κ fit would not rescue it.
- Losses concentrate in two strata:
  - focal_only clusters (`shrunk` 1,405 vs `sp_eb` 1,152 bits at m = 5): the genus prior says these are
    rare;
  - congener_only clusters (167 vs 94): the genus prior says these may be present.
- Shrinkage helps only marginally, and only in the shared stratum: the 50/50 mix is 1% better than
  `sp_eb` at m = 5.

Caveats:

- Only 2 genera are tested. Marinisoma is MAG-only (completeness is not modelled in B), and
  Flavobacterium is diverse.
- Orthologs below 50% identity split into focal_only and congener_only clusters. That split penalises
  `shrunk` alone, so the result is partly a mapping artefact.
- The full run (40 genera) decides it. Two follow-ups are worth adding if the full run agrees:
  - shrink the species-wide Beta(a, b) towards the genus per cluster, instead of replacing it;
  - restrict genus information to shared clusters.

## Runtime and memory plan for grex

- A takes ≤ 25 min per species on 1 thread (up to N = 400), plus about 10-25 min for N = 800 in the 46
  species with more than 800 genomes. Peak RSS is under 1 GB at 318 genomes. A 1,500-genome Rtab export
  is estimated at 2-4 GB (not tested).
- B takes ≤ 10 min per genus, using 2 threads for MMseqs2 only, with RSS ≤ 1.5 GB.
- The whole job is 331 tasks and about 100-150 CPU-h. It runs on 2 genoa nodes, with 180 workers each and
  one wave of tasks per node. Expected wall time is 30-60 min, inside the 2-h limit.
- `--mem=600G` per node leaves at least 3 GB per worker.
- Tasks are sorted most expensive first and spread round-robin over the array elements.
- Each task has a 6,900 s timeout. A tasks write their outputs after every subsample.
- Finished tasks leave a `.done` marker. Resubmitting the same command resumes and skips them.

## Exact submission (coordinator)

Files to copy: code from `/data/dev/metapanclust/experiments/pangbank_scale/`, and the selection already
made locally in `/data/scratch/metapanclust_data/pangbank_scale/sel/`. Both go to grex under
`BASE=/home/rec3141/scratch/metapanclust/pangbank_scale`. Use the usual route to grex
(via concentration); `grex:` below stands for that host.

```bash
# 1. from the workstation
ssh grex 'mkdir -p /home/rec3141/scratch/metapanclust/pangbank_scale/{code,sel,h5,out,logs}'
scp /data/dev/metapanclust/experiments/pangbank_scale/{select_species.py,download.py,worker.py,aggregate.py,run_grex.sbatch} \
    grex:/home/rec3141/scratch/metapanclust/pangbank_scale/code/
scp /data/scratch/metapanclust_data/pangbank_scale/sel/{catalog.tsv,species_A.tsv,genera_B.tsv,tasks.tsv,downloads.tsv} \
    grex:/home/rec3141/scratch/metapanclust/pangbank_scale/sel/

# 2. on the grex login node: sequential download (~455 files, ~45 GB, ~1-2 h at 7 MB/s; resumable, md5-checked)
cd /home/rec3141/scratch/metapanclust/pangbank_scale
export PATH=/opt/singularity-ce-4.4.1/bin:$PATH
SIF=/home/rec3141/scratch/metapanclust/sif/metapanclust-tools-0.1.sif
nohup singularity exec -B /home/rec3141/scratch $SIF python code/download.py sel/downloads.tsv h5 1 > logs/download.log 2>&1 &
# finished when the last line of logs/download.log reads "done: 455/455 present; failed: none"
# (rerun the same command if any failed)

# 3. on the grex login node, after the download: 2 full nodes
cd /home/rec3141/scratch/metapanclust/pangbank_scale
sbatch --array=0-1 code/run_grex.sbatch
# (3 nodes: sbatch --array=0-2 --export=ALL,NCHUNK=3 code/run_grex.sbatch)
# per-task logs: logs/tasks/{A_<id>,B_<genus>}.log; node logs: logs/pbscale_<job>_<i>.out
# resubmitting the same sbatch command skips finished tasks (out/A/<id>.done, out/B/<genus>.done)

# 4. aggregate (login node or any node; seconds)
singularity exec -B /home/rec3141/scratch $SIF python code/aggregate.py out sel/catalog.tsv > out/aggregate.txt
```

Outputs:

- `out/A/<id>.summary.tsv` and `.targets.tsv.gz`;
- `out/B/<genus>.summary.tsv` and `.kappa.tsv`;
- `out/agg/A_per_species.tsv`, `A_by_N.tsv`, `A_by_phylum_N50.tsv`, `B_per_focal.tsv` and `B_by_m.tsv`;
- `out/aggregate.txt`, which holds the markdown tables.

Everything under `out/work/` is temporary and is removed per task.

To regenerate the selection on grex instead (≈ 1 min of API paging, login node):
`singularity exec -B /home/rec3141/scratch $SIF python code/select_species.py sel`.

## Files

- Code: `experiments/pangbank_scale/`:
  - `select_species.py`
  - `download.py`
  - `worker.py`
  - `run_grex.sbatch`
  - `aggregate.py`
- Local data: `/data/scratch/metapanclust_data/pangbank_scale/`:
  - `sel/`: selection and catalogue;
  - `h5/`: 15 test files;
  - `test/out/`: test outputs and `agg/`;
  - `test/*.log`: worker logs.
