# Operon evidence from conserved neighbourhoods vs reference-panel distance (2026-11-18)

Question: how does Spacedust conserved-pair evidence for E. coli K-12 operons depend on how distant,
and how many, the comparison genomes are? What panel should a reference catalogue use?

## Result

- **Close genomes are useless, not just weaker.** Against any genome inside Enterobacteriaceae, "both
  genes in one Spacedust cluster" has precision 0.62-0.67. That is the base rate of the labelled set (1,727 of 2,820
  same-strand neighbours = 0.61), so the call carries almost no information. Conservation by descent
  covers TU boundaries: E. coli O157 Sakai alone "conserves" 666 of the 798 different-TU pairs with a
  50-300 bp gap. Adding such genomes to the panel lowers the combined rule from 0.850 (distance alone) to 0.69-0.74.
- **Precision rises steadily with distance** while per-genome recall falls (table below). Evidence becomes
  useful (precision >= 0.87) only outside GTDB Enterobacterales, i.e. from Pseudomonadales/other
  Gamma orders outward.
- **Best panel: many genomes, all outside the query's (GTDB) order.** 37 genomes from Betaproteobacteria,
  Alphaproteobacteria, 13 other bacterial phyla and 3 archaea (`ge_level6`): conserved-pair recall
  **0.487 at precision 0.897**, against 0.346-0.378 at 0.89-0.90 for the 5-6 distant genomes used before.
  Combined rule `gap <= 50 OR (conserved AND gap <= 300)`: **accuracy 0.878** (precision 0.908, recall
  0.891), against 0.850 for distance alone and 0.871-0.874 for the old panel.
- **Number of genomes:** with distant genomes only, precision stays flat (0.92-0.94) and recall grows
  roughly with log(n): 1 → 0.08, 5 → 0.22, 10 → 0.29, 20 → 0.36, 30 → 0.40 (random draws from the
  Alpha + other-phyla + Archaea pool, mean of 30). It has not saturated at 30. Combined accuracy goes from 0.855 to 0.876 over
  the same range. One extra genome from a nearer level (e.g. Vibrio) buys more recall than five
  far ones, but costs precision.
- **Weighting close genomes does not rescue them.** A cross-validated logistic model (gap bins plus a
  "conserved at level L" indicator per level) reaches 0.880 accuracy with levels 4-9 only and 0.878 with
  all levels. With levels 0-3 only it gets 0.845-0.860, about the same as gap alone (0.850). Fitted log-odds per level: L0 -0.06,
  L1 -0.25, L2 +0.06, L3 +0.19, L4 -0.01, L5 +0.28, L6 +0.51, L7 +0.55, L8 +1.38, L9 +0.46 (archaea:
  3 genomes, few hits).
- **The gain is capped.** With the best panel, the combined rule still makes 91 false calls at gap <= 50
  (adjacent different-TU pairs, which the distance rule makes too) and misses 189 same-TU pairs (170 at
  gap 51-300 that no distant genome conserves, mostly E. coli-specific operons, plus 19 at gap > 300).
  Getting much past ~0.88 needs a different signal (promoters/terminators, expression), not more genomes.

## Per distance level (union over the level's genomes, conserved = co-membership in one cluster)

68 genomes in total (37 first, 31 added at levels 5-8). Levels follow NCBI taxonomy, with GTDB placement noted.

| Level | Genomes | n | Cons. precision | Cons. recall | Combined acc | Comb. prec | Comb. recall |
|---|---|---|---|---|---|---|---|
| 0 | E. coli O157 Sakai, CFT073, Shigella flexneri 2a | 3 | 0.620 | 0.932 | 0.696 | 0.673 | 0.981 |
| 1 | Escherichia fergusonii, albertii, marmotae | 3 | 0.622 | 0.892 | 0.704 | 0.680 | 0.976 |
| 2 | Enterobacteriaceae: Salmonella, Klebsiella, Citrobacter, Enterobacter, Cronobacter | 5 | 0.656 | 0.863 | 0.735 | 0.708 | 0.964 |
| 3 | Enterobacterales other families: Yersinia, Serratia, Pectobacterium, Proteus, Erwinia | 5 | 0.707 | 0.704 | 0.795 | 0.776 | 0.936 |
| 4 | Vibrio, Haemophilus, Shewanella (other NCBI orders, GTDB Enterobacterales) | 3 | 0.799 | 0.436 | 0.852 | 0.873 | 0.888 |
| 5 | Other Gamma orders (Pseudomonas, Acinetobacter, Xanthomonas, Legionella, Methylococcus, Francisella, Coxiella, Alcanivorax, Cellvibrio, Nitrosococcus, Allochromatium, Hahella) | 12 | 0.870 | 0.422 | 0.871 | 0.907 | 0.880 |
| 6 | Betaproteobacteria (GTDB: Gamma order Burkholderiales): Neisseria, Burkholderia, Ralstonia, Bordetella, Nitrosomonas, Cupriavidus, Thiobacillus | 7 | 0.893 | 0.354 | 0.869 | 0.915 | 0.866 |
| 7 | Alphaproteobacteria: Caulobacter, Cereibacter, Agrobacterium, Zymomonas, Rhodopseudomonas, Bradyrhizobium, Paramagnetospirillum | 7 | 0.929 | 0.281 | 0.869 | 0.927 | 0.853 |
| 8 | Other bacterial phyla (Firmicutes, Actinobacteria, Bacteroidota, Cyanobacteria, Thermotogota, Campylobacterota, Desulfobacterota, Myxococcota, Chlorobi, Aquificota, Chloroflexota, Planctomycetota, Deinococcota) | 20 | 0.938 | 0.340 | 0.873 | 0.927 | 0.862 |
| 9 | Archaea: Methanocaldococcus, Halobacterium, Saccharolobus | 3 | 0.969 | 0.036 | 0.853 | 0.939 | 0.813 |

Nested panels (everything at or beyond a level):

| Panel | n | Cons. precision | Cons. recall | Combined acc | Comb. prec | Comb. recall |
|---|---|---|---|---|---|---|
| distance rule only | 0 | - | - | 0.850 | 0.939 | 0.808 |
| old baseline 5 (Bsu, Pae, Mtu, Syn, Mja; NCBI GFF, pairwise) | 5 | 0.886 | 0.346 | 0.871 | 0.917 | 0.867 |
| all 68 | 68 | 0.621 | 0.967 | 0.689 | 0.667 | 0.986 |
| >= level 2 (no Escherichia) | 62 | 0.659 | 0.900 | 0.732 | 0.704 | 0.970 |
| >= level 3 | 57 | 0.719 | 0.773 | 0.796 | 0.772 | 0.946 |
| >= level 4 | 52 | 0.826 | 0.624 | 0.857 | 0.863 | 0.911 |
| >= level 5 (outside GTDB Enterobacterales) | 49 | 0.874 | 0.540 | 0.874 | 0.898 | 0.898 |
| **>= level 6** (Beta and beyond) | 37 | **0.897** | **0.487** | **0.878** | 0.908 | 0.891 |
| >= level 7 (outside Pseudomonadota-Gamma sensu GTDB) | 30 | 0.923 | 0.402 | 0.876 | 0.921 | 0.872 |
| >= level 8 | 23 | 0.937 | 0.342 | 0.874 | 0.927 | 0.862 |

Requiring conservation in k >= 2 genomes of a level raises precision only slightly (level 3: 0.707 → 0.711) and
loses recall. A majority vote over close genomes does not fix them: conservation by descent is shared by
all of them. A stricter "adjacent and same orientation in the target too" test adds 1-2 points of precision
for 2-5 points of recall at every level. It does not change the picture.

## Answers for panel design

1. **Distance floor, not ceiling.** For a query genome, count neighbourhood conservation only in genomes
   outside its own GTDB order (here: outside Enterobacterales incl. Vibrio/Haemophilus/Shewanella).
   Same-genus, same-family and same-order genomes inflate false positives to the base rate.
2. **Then add breadth.** Recall keeps growing with the number of distant genomes (no plateau by 30-37);
   spreading over classes and phyla is better than many genomes from one clade.
3. **For a GTDB-wide catalogue:** score each pair of neighbouring genes by the number of independent distant
   lineages (e.g. distinct orders outside the query's order) that conserve it, not by raw genome counts.
   Within-order genomes belong in the prevalence/pangenome side of the model, not in the operon evidence.
4. Expected ceiling of neighbourhood + distance evidence for E. coli: ~0.88 accuracy on adjacent same-strand
   pairs.

## Methods

- Genomes: 68 RefSeq reference assemblies (complete, NCBI datasets v2, NCBI GFF3 CDS annotation)
  listed in `experiments/phylo_scaling/panel_ext.tsv` (first 37 = `panel.tsv`). Query: E. coli K-12
  MG1655 with the Bakta annotation used in `first_tests.md`.
- Spacedust 2d68c85, MMseqs2 mode, default parameters. The E. coli set DB is the query and each panel
  genome is a separate target (`run_pairwise.sh`, 9-16 s per genome), so any panel is the union of
  pairwise runs. Control: one joint run against the 14 genomes of levels >= 6 (`run_joint.sh`, scored with
  `scripts/eval_pairs.py`) gives 729 conserved calls, precision 0.911, recall 0.384, combined 0.875.
  The union of the pairwise runs gives 736, 0.910, 0.388, 0.876, so the union is equivalent.
- Scoring (`score_panels.py`): same pair definition and rules as `scripts/eval_pairs.py`, on
  `bench/eco_pairs.tsv` (2,820 labelled adjacent same-strand EcoCyc pairs). Random-panel curves:
  30 draws per size. `fit_model.py`: L2 logistic regression, 5-fold CV over contiguous genome blocks
  (neighbouring pairs of one operon stay in one fold).
- Caveats: a single query genome with curated TUs. NCBI and Bakta gene calls differ slightly (the old
  panel reproduces 0.346 vs 0.378 recall). Archaea and very reduced genomes contribute few clusters.
  Results are for E. coli and should be repeated for a second organism with curated operons (B. subtilis,
  SubtiWiki) before fixing a panel rule.

## Files

Code: `experiments/phylo_scaling/` (`run_pairwise.sh`, `run_joint.sh`, `score_panels.py`, `fit_model.py`,
`panel.tsv`, `panel_ext.tsv`, `bin/spacedust` = niced docker wrapper).
Data (workstation): `/data/scratch/metapanclust_data/phylo_scaling/`: `genomes/` (FASTA + GFF3), `runs/pw/<name>/clusters.tsv`,
`runs/joint_ge6/`, `results/ext_panels.tsv` (all panels), `results/ext_curve.tsv` (size curves),
`results/ext_matrix.tsv` (pair x genome conservation matrix), `results/ext_model_cv.txt`; `pw_*` = first 37 only.
Compute: 68 pairwise runs, about 15 CPU-minutes in total at 2-3 threads, under 2 GB disk. No cluster job needed.

Reproduce:
```
cd /data/scratch/metapanclust_data/phylo_scaling
E=/data/dev/metapanclust/experiments/phylo_scaling
$E/run_pairwise.sh $E/panel_ext.tsv genomes runs/pw 2
python3 $E/score_panels.py $E/panel_ext.tsv runs/pw ../bench/eco_pairs.tsv \
  /data/dev/metacyc_benchmark/ecoli/mag_out/annotation/bakta/basic/annotation.tsv results/ext
python3 $E/fit_model.py $E/panel_ext.tsv results/ext_matrix.tsv
```
