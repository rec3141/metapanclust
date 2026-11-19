# Operon evidence vs reference-panel distance: B. subtilis replication (2026-11-18)

Question: does the E. coli rule from `phylo_scaling.md` hold for a Bacillota query? The rule: count
conserved-neighbourhood evidence only from genomes outside the query's GTDB order, and use many lineages.
Query: B. subtilis 168 with SubtiWiki operons, 92 comparison genomes in 10 distance levels.

## Result

- **Close genomes are useless here too.** Against B. subtilis strains and the other subtilis-group species
  (GTDB g__Bacillus), "both genes in one Spacedust cluster" has precision 0.74-0.75. That equals the base rate
  of the labelled set (1,928 of 2,621 = 0.736). All 92 genomes pooled: 0.738.
- **The precision floor sits further out than the GTDB order.** GTDB r226 splits NCBI Bacillales into about
  ten orders (Bacillales, Bacillales_A/B/D/H, Paenibacillales, Staphylococcales, Exiguobacteriales,
  Brevibacillales, Alicyclobacillales). Two groups give almost the same precision:
  - other families inside o__Bacillales (B. cereus/anthracis, Priestia, Geobacillus, Anoxybacillus): 0.846;
  - the other GTDB orders of NCBI Bacillales (Halalkalibacterium, Shouchella, Oceanobacillus, Lysinibacillus,
    Paenibacillus, Brevibacillus, Exiguobacterium, Alicyclobacillus, Staphylococcus, Listeria): 0.854.

  Precision first reaches E. coli's useful level (>= 0.90) with other Bacilli (Lactobacillales,
  Erysipelotrichales: 0.93) and other Bacillota classes (0.91). Other phyla give 0.95-0.96.
  "Outside the query's GTDB order" is therefore not distant enough for B. subtilis. The working boundary is
  roughly the class, or outside the NCBI order. A fixed rank does not mean a fixed divergence across GTDB.
- **Breadth still helps, with the same log-like recall growth.** Random panels from the distant pools keep
  precision flat while conserved-pair recall grows:
  - other phyla + archaea (L6-9): 1 → 0.08, 5 → 0.16, 10 → 0.21, 20 → 0.25, 40 → 0.29, at precision
    0.95-0.96;
  - all outside NCBI Bacillales (L4-9): 1 → 0.09, 5 → 0.22, 10 → 0.27, 20 → 0.34, 40 → 0.41, 60 → 0.43, at
    precision 0.95 → 0.91.

  As with E. coli, recall has not saturated at 40-60 genomes.
- **But conservation barely improves accuracy for B. subtilis.** Accuracy here is for calling an adjacent
  same-strand pair "same operon".
  - **Gap-only baseline: 0.871** (gap <= 100 bp; precision 0.944, recall 0.877). The E. coli thresholds do
    not transfer. B. subtilis same-operon spacers are longer (75th percentile 64 bp, 90th 114 bp), so
    gap <= 50 gives only 0.774. The E. coli combined rule `gap<=50 OR (conserved AND gap<=300)` reaches
    0.80-0.81 at best.
  - **Best panel, with thresholds re-tuned** (rule `gap<=T1 OR (conserved AND gap<=T2)`): all genomes
    outside GTDB class Bacilli plus Pseudomonadota and beyond (`ge_level7`, 54 genomes, T1=100, T2=250):
    **0.877** in-sample, **0.874** with 5-fold block CV. Other distant panels (`ge_level4` to `ge_level8`)
    give 0.874-0.877 in-sample and 0.871-0.874 in CV.
  - Cross-validated logistic model: gap bins alone 0.871; best (gap + per-level indicators or per-genome
    indicators) 0.876-0.877. Close levels 0-2 alone: 0.863-0.872.
  - Net gain from conservation: +0.003 to +0.006, against +0.03 for E. coli (0.848 → 0.878 CV under the same
    re-tuned scoring).
  - Why: past the gap threshold (gap > 100 bp) there are 238 same-operon and 593 different-operon pairs.
    Distant panels mark only 13-18% of those same-operon pairs as conserved, at 57-63% precision. For E. coli
    the comparable figures are 45-49% at 62-67%. B. subtilis long-spacer operons are mostly not conserved
    outside the Bacillales/Bacilli.
- **Verdict on the rule.** Two parts of it replicate: exclude close genomes, and prefer many distant
  lineages. The cut-off "outside the query's GTDB order" does not replicate. For B. subtilis, half of the
  NCBI-order neighbours outside o__Bacillales still inflate false positives. Set the floor by divergence (for
  example AAI or a marker-protein distance, or "outside the class" as a crude proxy), not by GTDB order.
  Expect the operon-accuracy gain from neighbourhood evidence to depend on the organism: about +3 points for
  E. coli, under +1 point for B. subtilis. The intergenic-distance threshold must be fitted per query, or
  per clade (50 bp for E. coli, 100 bp for B. subtilis).

## Per distance level

Union over the level's genomes. Conserved = co-membership in one Spacedust cluster. "E. coli rule" =
`gap<=50 OR (cons AND gap<=300)`. "Tuned" = best `gap<=T1 OR (cons AND gap<=T2)`, in-sample / 5-fold block CV.

| Level | Genomes | n | Cons. prec | Cons. rec | E. coli rule acc | Tuned acc (CV) |
|---|---|---|---|---|---|---|
| - | gap only | 0 | - | - | 0.774 | 0.871 (0.871) |
| 0 | B. subtilis NCIB 3610, BSn5; B. spizizenii W23 | 3 | 0.738 | 0.993 | 0.772 | 0.871 (0.871) |
| 1 | B. velezensis, licheniformis, pumilus, atrophaeus (g__Bacillus) | 4 | 0.754 | 0.888 | 0.778 | 0.872 (0.872) |
| 2 | GTDB o__Bacillales, other families: B. cereus, B. anthracis, Priestia, Geobacillus, Anoxybacillus | 5 | 0.846 | 0.601 | 0.805 | 0.871 (0.865) |
| 3 | NCBI Bacillales in other GTDB orders: Halalkalibacterium, Shouchella, Oceanobacillus, Lysinibacillus (Caryophanales/Bacillales_A), Paenibacillus, Brevibacillus, Exiguobacterium, Alicyclobacillus, Staphylococcus, Listeria | 10 | 0.854 | 0.620 | 0.808 | 0.873 (0.868) |
| 4 | Other Bacilli: Lactococcus, Streptococcus, Enterococcus, Lactiplantibacillus, Erysipelothrix | 5 | 0.927 | 0.251 | 0.794 | 0.876 (0.876) |
| 5 | Other Bacillota: Clostridium, Clostridioides, Desulfitobacterium, Moorella, Veillonella, Caldanaerobacter | 6 | 0.905 | 0.322 | 0.796 | 0.874 (0.872) |
| 6 | Other Terrabacteria: M. tuberculosis, Streptomyces, Synechocystis, Chloroflexus, Deinococcus | 5 | 0.957 | 0.163 | 0.787 | 0.876 (0.874) |
| 7 | Pseudomonadota (Gamma incl. Beta, Alpha; first-wave genomes, one E. coli) | 40 | 0.960 | 0.248 | 0.792 | 0.876 (0.874) |
| 8 | Other phyla (Bacteroidota, Campylobacterota, Desulfobacterota, Myxococcota, Chlorobiota, Aquificota, Planctomycetota, Thermotogota) | 11 | 0.946 | 0.211 | 0.787 | 0.876 (0.873) |
| 9 | Archaea (Methanocaldococcus, Halobacterium, Saccharolobus) | 3 | 0.988 | 0.043 | 0.777 | 0.873 (0.872) |

Single genomes show the same pattern:
- levels 0-1: 0.73-0.80;
- level 2: 0.83-0.86;
- level 3: 0.84-0.91 (Bacillales_D/H 0.84-0.86; Caryophanales, Paenibacillales, Brevibacillales,
  Exiguobacteriales 0.88-0.89; Alicyclobacillus, Staphylococcus, Listeria 0.90-0.91);
- level 4: 0.93-0.97;
- other phyla: 0.95-0.98.

Nested panels (everything at or beyond a level):

| Panel | n | Cons. prec | Cons. rec | E. coli rule acc | Tuned acc (CV) | T1/T2 |
|---|---|---|---|---|---|---|
| all 92 | 92 | 0.738 | 0.997 | 0.773 | 0.871 (0.871) | 0/100 |
| >= level 2 (no g__Bacillus) | 85 | 0.846 | 0.696 | 0.811 | 0.872 (0.872) | 90/100 |
| >= level 3 (outside GTDB o__Bacillales) | 80 | 0.854 | 0.645 | 0.809 | 0.873 (0.868) | 90/100 |
| >= level 4 (outside NCBI Bacillales) | 70 | 0.905 | 0.446 | 0.804 | 0.874 (0.871) | 100/150 |
| >= level 5 (outside class Bacilli) | 65 | 0.911 | 0.413 | 0.801 | 0.875 (0.873) | 100/150 |
| >= level 6 (outside Bacillota) | 59 | 0.943 | 0.309 | 0.795 | 0.876 (0.874) | 100/250 |
| **>= level 7** | 54 | **0.950** | **0.295** | 0.795 | **0.877 (0.874)** | 100/250 |
| baseline5 (E. coli, Pae, Mtu, Syn, Mja) | 5 | 0.962 | 0.169 | 0.785 | - | - |

Requiring k >= 2-5 genomes does not rescue close panels: `ge_level2`, k=5, gives 0.861. Requiring the hit genes
to be adjacent in the target too adds 1-3 precision points (`ge_level3` 0.864, `ge_level6` 0.970) and costs
3-4 recall points, the same as for E. coli.

Cross-validated logistic model (`results/bsub_model_cv.txt`):

| Features | Acc |
|---|---|
| gap bins only | 0.871 |
| gap + any(level >= 3) | 0.866 |
| gap + count(level >= 3) | 0.874 |
| gap + per-level levels 3-9 | 0.876 |
| gap + per-level levels 0-2 | 0.863 |
| gap + per-genome indicators | 0.877 |

Full-data log-odds per level are all small (-0.12 to +0.72), so no level carries much signal beyond gap.

## Methods

- **Labels.** SubtiWiki operons, fetched 2026-11-18:
  - v5 REST API `https://subtiwiki.uni-goettingen.de/v5/api/operon/`: 2,310 operons, 1,038 with more than one
    gene;
  - per-gene `/v5/api/gene/<id>` for strain-168 coordinates (3,925 genes);
  - the v4 CSV export (`/v4/operon/exporter`, 2,295 operons) kept for reference.

  SubtiWiki has no release number; the fetch date stands in for a version. Bakta CDS were mapped to SubtiWiki
  genes by best reciprocal overlap >= 0.8 on NC_000964.3 (3,776 of 4,233 CDS). Same-strand neighbours were
  labelled 1 if they share an operon. They were labelled 0 if both genes are in operons but never the same
  one, and left unlabelled otherwise. Result: 3,097 pairs, 2,621 labelled, 1,928 positive.
  Script: `bsub_pairs.py`, adapted from `scripts/operon_pairs.py`.
- **Genomes.** 92 RefSeq complete assemblies with NCBI GFF3 CDS (`panel_bsub.tsv`):
  - 62 reused from phylo_scaling: E. coli strains/species redundancies dropped; Bsu dropped as the query
    itself; Sau, Lmo, Cac reassigned to levels 3/3/5;
  - 30 downloaded with the NCBI datasets v2 API (`get_genomes.sh`).

  GTDB placements are from `/data/scratch/refdbs/GTDB/release226/taxonomy/gtdb_taxonomy.tsv`.
- **Runs.** Spacedust 2d68c85 (first-wave niced docker wrapper), MMseqs2 mode, defaults. The B. subtilis
  Bakta set DB is the query against each genome separately (`run_pairwise.sh`, 5-10 s per genome, 3
  threads). About 10 min wall time (about 30 CPU-min) in total.
- **Scoring.**
  - `score_panels.py`: first-wave scorer with B. subtilis pools; fixed E. coli 50/300 rule; random panels
    with 30 draws.
  - `score_tuned.py`: per-panel threshold grid (T1 0-200, T2 100-1000), with 5-fold contiguous-block CV;
    distant-pool recall curves. Run on the E. coli matrix too, for the +0.03 comparison.
  - `fit_model.py`: first-wave logistic model, distant = level >= 3.
- **Caveats.**
  - SubtiWiki "operons" include overlapping and alternative transcripts. A pair counts as same-TU if any
    operon holds both genes, as with EcoCyc.
  - The base rate (0.74) is higher than for E. coli (0.61).
  - Level 0 includes NCIB 3610, nearly identical to 168.
  - Only one Bacillota query. The "class-level floor" is an observation for B. subtilis, not yet a rule.

## For grex

Nothing needs grex. All runs fit the workstation limits: about 30 CPU-min, 0.75 GB disk.

## Files

Code: `experiments/phylo_scaling_bsub/`:
- `fetch_subtiwiki.sh`, `bsub_pairs.py`, `get_genomes.sh`, `panel_bsub.tsv`;
- `run_pairwise.sh`, `score_panels.py`, `score_tuned.py`, `fit_model.py`.

Data: `/data/scratch/metapanclust_data/phylo_scaling_bsub/`:
- `dl_subtiwiki/` (operon.json, operon_v4.csv, gene_coords.tsv, genes/*.json);
- `bench/bsub_pairs.tsv`;
- `genomes/` (30 downloaded, 62 symlinks to phylo_scaling);
- `runs/pw/<name>/clusters.tsv`;
- `results/bsub_{panels,curve,matrix}.tsv`, `bsub_tuned.tsv`, `bsub_model_cv.txt`, `eco_tuned.tsv`.

Reproduce:
```
cd /data/scratch/metapanclust_data/phylo_scaling_bsub
E=/data/dev/metapanclust/experiments/phylo_scaling_bsub
B=/data/dev/metacyc_benchmark/bsub/mag_out/annotation/bakta/basic/annotation.tsv
$E/fetch_subtiwiki.sh dl_subtiwiki
python3 $E/bsub_pairs.py $B dl_subtiwiki bench/bsub_pairs.tsv
$E/get_genomes.sh $E/panel_bsub.tsv ../phylo_scaling/genomes genomes dl_ncbi
$E/run_pairwise.sh $E/panel_bsub.tsv genomes runs/pw 3
python3 $E/score_panels.py $E/panel_bsub.tsv runs/pw bench/bsub_pairs.tsv $B results/bsub
python3 $E/score_tuned.py $E/panel_bsub.tsv results/bsub_matrix.tsv results/bsub_tuned.tsv
python3 $E/fit_model.py $E/panel_bsub.tsv results/bsub_matrix.tsv
```
