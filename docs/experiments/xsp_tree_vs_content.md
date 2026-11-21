# Between species: GTDB tree vs gene content (design, local test, grex submission, 2026-11-20)

Status: done. grex job 8194548 (18 clades, 600 species, about 45 min on one node). Results below; full tables in
[xsp_tree_vs_content_results.txt](xsp_tree_vs_content_results.txt).

## Question

[tree_vs_context.md](tree_vs_context.md) found that within a species, gene-content relatedness predicts held-out
content better than ANI or core-protein distance. The case the PI raised is a different scale: coral-symbiont
Marine Group I archaea have shifted GC, sit as outliers on sequence trees, and look closely related by gene
content. That is a statement about species and genera, so this run uses species as tips:

1. Across species in a family, which relatedness best predicts a held-out species' gene content: the GTDB
   marker tree, core-protein AAI, or gene-content similarity?
2. Are GC-shifted species placed differently by the tree and by gene content, and which one is right for them
   (judged by prediction of their content)?
3. How much does gene-gene context add beyond relatedness at this scale?

## Design

**Clades** (`select_clades.py`): 17 GTDB families with ≥ 20 PanGBank species (spread over Pseudomonadota,
Bacteroidota, Bacillota, Actinomycetota, Campylobacterota, Cyanobacteriota; GC-diverse and marine families
included) plus class Nitrososphaeria (16 species, including Nitrosopumilaceae = Marine Group I). Up to 40 species
per clade, randomly drawn among species with ≤ 250 genomes; 600 species and 31,000 genomes in all. PanGBank needs
≥ 15 genomes per species, so the coral-symbiont MGI themselves are not in it; the class is still the
closest available test.

**Inputs.** Each PanGBank `.h5` stores a representative protein per gene family, and from gene → genome
membership the family's frequency in the species. Per species: families in ≥ 2 genomes and ≥ 5% of genomes (all
persistent families kept), median genome GC and CheckM2 completeness (`worker.py prep`).

**Orthogroups (OGs).** MMseqs2 easy-cluster of all family representatives in the clade (40% identity, 80%
coverage both ways). F[species, OG] = highest frequency among the species' families in the OG; a species carries
the OG if F ≥ 0.5. A target OG's homology neighbours (representatives that hit it at e ≤ 1e-3) are masked with
it, so split families and paralogs cannot leak the target into the predictors.

**Relatedness between species.**

| Measure | How | Expected weakness |
|---|---|---|
| `gtdb` | GTDB r226 bac120/ar53 tree pruned to the clade's species representatives; patristic distance, and the tree itself (not a UPGMA rebuild) for the tree predictor | compositional bias in the marker alignment (long branches for GC-shifted lineages), marker HGT |
| `aai` | 1 − AAI on core proteins (families in ≥ 90% of genomes): best hit per query and target species (≥ 50% coverage), mean identity, averaged over both directions; UPGMA tree | amino-acid composition also shifts with GC |
| `jac` | 1 − Jaccard of OG presence over variable OGs, every target OG and its neighbours removed; UPGMA tree | HGT and convergent gene loss (streamlining) |

**Prediction.** Leave one species out; predict its presence of up to 400 variable OGs (carried by 2 … S − 2
species), using the same models as tree_vs_context: prevalence, kernel and tree (gain/loss chain, Bayesian
average over rate and stationary frequency) for each measure, stacked relatedness (`relseq` = gtdb + aai,
`relall` = + jac), context (`lr`, L2 logistic on the other variable OGs, target and neighbours masked), and
50/50 averages. The tip model uses completeness 0.99 (a species-level majority call is nearly free of detection
loss).

**Diagnostics** (`<clade>.diag.tsv`, `<clade>.matcor.tsv`): per species, GC deviation from the clade median,
GTDB root-to-tip z, nearest neighbour and its GC difference under each measure, 5-NN overlap and row correlation
between measures. Per clade: distance correlations, correlation of |ΔGC| with each distance, and with the GTDB
(or AAI) distance after removing the part explained by gene content (does the sequence measure carry a GC signal
that gene content does not?). `aggregate.py` §5 asks directly whether gene content gains over the GTDB tree on
GC-shifted species; §6 lists discordant species ("tree outlier, content typical" is the MGI pattern).

## Local test: Nitrososphaeria (16 species, 4 CPU, 100 s)

16,008 OGs, 3,592 variable, 400 targets. % of prevalence log loss removed: tree_gtdb 45.9, tree_jac 43.9,
tree_aai 37.3, kernels 16-32, context (`lr`) 61.2.

- The GTDB tree and gene content are tied (GTDB better for 9/16 species, p 0.67); both clearly beat core AAI
  (p < 0.001). Trees beat kernels, as within species.
- The clade's GC outlier is UBA183 (61-63% GC in a class at about 35%). For both UBA183 species the GTDB tree
  predicts content better than gene content (0.41 vs 0.54 bits). Its GTDB root-to-tip length is ordinary
  (z 0.39), and its nearest neighbour is the same under every measure (the other UBA183).
- Only Nitrosopumilus limneticus fits the "tree outlier, content typical" pattern (GC 30.4, root-to-tip
  z 1.6). Gene content predicts it slightly better (0.31 vs 0.33 bits).
- Context beats every relatedness predictor at this scale (61 vs 46%), more than within species. With 16 tips the
  trees carry little information, and co-occurring OGs (operons, split genes not caught by the homology mask)
  are strong predictors. That is a caveat on context, not yet a finding.

One clade of 16 species decides nothing; the grex run has 18 clades and 600 species.

## grex run

`/home/rec3141/scratch/metapanclust/xsp/` (code/, sel/, h5/ = 83 links to pangbank_scale plus 517 new
downloads). `logs/download.log`; the download wrapper submits `code/run_grex.sbatch` (one genoa node: prep for
all species with 150 workers, then all 18 clades at once with 10 threads each). Results land in `out/res/`;
aggregate with `aggregate.py out/res`.

## Results (grex, 18 clades, 600 species)

**1. Between species, the GTDB tree is the best relatedness predictor, but only just ahead of gene content.**
Median % of prevalence log loss removed over clades: tree_gtdb 41.8, tree_jac 41.5, tree_aai 37.6; kernels
17-29. The GTDB tree wins in 13/18 clades and for 64% of species (median 0.007 bits per target better than the
gene-content tree, p 1e-12). Core AAI is worst (GTDB better for 79% of species). This reverses the within-species
result, where gene content beat sequence measures for 62-98% of species. Trees beat kernels at both scales.

**2. Stacking relatedness measures hurts at this scale** (relseq 38.0, relall 38.1, below tree_gtdb alone):
6 logits fitted on 15-39 training species overfit. Pick one tree; do not stack.

**3. Context still adds.** `lr` alone removes 44.4%, more than any relatedness predictor, and beats `relall`
for 71% of species; the 50/50 average with relatedness is best (46.7, +1.4 over the better of the two). The
caveat from the test stands: some of `lr`'s gain may be co-occurring split OGs the homology mask misses.

**4. GC shifts do not make the GTDB tree worse on average.**
- |GC deviation| does not predict whether gene content beats the GTDB tree (Spearman 0.05, p 0.23; 11/18
  clades positive). In the > 10-point GC class gene content wins for 44% of species against 35% below 2 points:
  a weak trend in the MGI direction, not a reversal.
- GC deviation barely tracks GTDB root-to-tip length (Spearman 0.08, p 0.06), so composition-driven long
  branches are not common in these GTDB trees.
- After removing the part explained by gene content, the GTDB distance carries little GC signal (median
  Spearman 0.08 with |ΔGC|).
- The ten "tree outlier, content typical" species are mostly sister pairs that sit together on a long branch
  (Gilliamella apicola, Campylobacter_B ureolyticus, Bifidobacterium piotii/vaginale, Alteromonas). Their
  nearest neighbour is the same under both measures, and the GTDB tree predicts 9 of the 10 better.

**5. MGI-like cases exist but are rare.** The clearest is *Taylorella equigenitalis* (Burkholderiaceae, GC 37.5
in a family at about 64): its GTDB nearest neighbour is Achromobacter (30 GC points away), its gene-content
nearest neighbour is UBA954 (16.5 points away), and gene content predicts its content better (0.31 vs 0.36
bits). UBA629 (Lachnospiraceae, GC −14.5, root-to-tip z 3.4, the longest branch in the run) is similar (0.32 vs
0.34). The opposite also occurs: the GC-shifted UBA183 pair in Nitrososphaeria (61-63% GC) and Gilliamella
(−21 GC) are predicted better by the GTDB tree.

**6. Nitrososphaeria / Marine Group I.** GTDB and gene content tie (46.0 vs 44.0%). Among the Nitrosopumilaceae,
gene content is slightly better for N. limneticus (0.28 vs 0.29 bits) and Nitrosopelagicus brevis (0.33 vs
0.36, the one MGI species whose nearest neighbour differs: Nitrosopumilus limneticus by content,
Nitrosopelagicus sp013390765 by GTDB). Neither is a large effect, and the coral-symbiont MGI are not in PanGBank.

## What this means for the model

- **Use the measure that fits the scale.** Within a species: gene-content tree (tree_vs_context). Across species
  (genus/family priors for species with few or no genomes, level-1 shrinkage): the GTDB tree, which is free
  for any GTDB-Tk-placed bin. The two are within 1% of each other between species, so a gene-content tree is an
  acceptable fallback for bins GTDB-Tk cannot place.
- **Flag discordance rather than pick a winner per species.** Species whose GTDB and gene-content nearest
  neighbours disagree (238/600) are where a single tree is least safe; averaging the two trees' predictions is
  cheap to test there (not done here).
- Your MGI pattern (sequence-tree outlier, close by gene content, gene content right) does occur (Taylorella,
  UBA629), but at family scale it is the exception, and GC shift alone does not predict it.
