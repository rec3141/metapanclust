# Spacedust structure mode (Foldseek + ProstT5 3Di) vs MMseqs2 mode (2026-11-18)

Question: does searching with ProstT5-predicted 3Di (Spacedust `clustersearch --search-mode 2`) find more
conserved gene neighbourhoods than MMseqs2 sequence search, and does that help operon prediction?

## Result

- **Structure mode finds more, and loses almost nothing.** On the two subsets that could be run here (E. coli
  against one distant genome, both modes on the identical subset), it kept all but one MMseqs2-covered gene and
  added more:
  - E. coli + M. jannaschii: genes in clusters 65 → 81 (+25%), co-clustered gene pairs 48 → 86, clusters
    23 → 26; no gene lost.
  - E. coli + Synechocystis: genes 254 → 276 (+9%), gene pairs 721 → 821, clusters 71 → 75; 1 gene lost.
- **The extra genes are real remote homologs at 15-35% identity**, almost all extending neighbourhoods MMseqs2
  had already found: ribosomal proteins next to the r-protein block, ABC transporter binding proteins and
  permeases (liv, art, mod, tau, fec), PstS, NADH dehydrogenase subunits, ATP synthase delta, carboxysome /
  microcompartment shell proteins.
- **Operon calls barely change.** Recall of conserved E. coli adjacent pairs inside the subsets rises
  (0.23 → 0.28 against M. jannaschii, 0.51 → 0.59 against Synechocystis) at unchanged precision (1.00, 0.93).
  But the combined rule (gap <= 50 OR conserved AND gap <= 300) keeps the same accuracy (0.904, 0.891), because
  14 of the 16 gained pairs are close (gap <= 50 bp) and the distance rule already calls them. What operon
  calling needs is far same-TU pairs (331 genome-wide). In the subsets structure mode added 1 of 35, plus 1 far
  different-TU pair (gap 290 bp) as a new false positive.
- **Cost is the blocker: ProstT5 on CPU is about 9 s per protein on 2 cores (~34 residues/s,
  ~18 core-seconds per protein).** The Foldseek search itself is cheap (seconds at this size). The six genomes
  (23,901 proteins, 7.45 M residues) would take ~2.5 days on 2 cores (long proteins cost more: the 3,973 aa one dominated the second subset), so only subsets were run. Chesterfield
  (1.5 M proteins) would need ~7,400 core-hours on CPU: GPU only.
- **The workstation GPU cannot be used.** The precompiled Foldseek GPU build targets Turing and newer; on the
  Quadro GV100 (Volta, sm_70) ProstT5 aborts with `ggml-cuda.cu:72: CUDA error`. A sm_70 source build is
  possible (`experiments/structure_mode/build_fsgpu_sm70.sh`, upstream Dockerfile recipe) but was stopped
  (slow on 2 cores, and the GPU is needed by the live nanopore basecalling run).

**Recommendation:** worth it for the cluster catalogue, not for operon calling. Structure mode adds 9-25% more genes to
conserved clusters between distant genomes and loses almost no MMseqs2 hits, which matters for
cross-phylum catalogue building (stage 1) and for divergent MAG lineages. For E. coli operon prediction against
these six genomes the gain is too small to change the combined rule. Run it only where ProstT5 can be on a
Turing-or-newer GPU (grex), compute 3Di once per protein, and keep it alongside the setDB (`setDB_ss`) so every
later clustersearch reuses it. Next check: the full six-genome run on grex (below) for the genome-wide numbers.

## Method

The MMseqs2-mode run (`bench/sd6_mmseqs`, `scripts/run_spacedust.sh`) gives 2,144 clusters covering 4,334 of
23,901 genes (18%; E. coli 31%, M. jannaschii 4%).

Subsets (`subset_regions.py`): take the MMseqs2-mode clusters between E. coli and one distant genome, expand each
cluster's genes by +-FLANK neighbours, write reduced GFF3s, and rerun **both** modes on the same subset setDB, so
the comparison is like for like. The subsets can show clusters being extended or new clusters near known ones,
not wholly new neighbourhoods elsewhere in the genome, and conservation is against one genome instead of five, so
absolute recall is lower than in the six-genome run.

| Subset | Proteins | Residues | E. coli labelled pairs |
|---|---|---|---|
| E. coli + M. jannaschii, flank 3 | 247 | 75,142 | 83 (64 same TU) |
| E. coli + Synechocystis, flank 1 | 453 | 141,985 | 165 (132 same TU) |

3Di was predicted with the CPU Foldseek (same commit as the Spacedust-bundled one, 4f5c00c) from the setDB
sequences (`setdb_to_fasta.py`, keys as headers) and attached to the setDB as `setDB_ss` (`attach_3di.py`, all
lengths checked equal), so gene ids stay identical to the MMseqs2 run. Search: `clustersearch --search-mode 2
--filter-self-match 1`. Everything niced and capped at 2 CPUs / 8 GB (`fscpu.sh`, `run_subset.sh`).

## Numbers

| Subset / mode | Clusters | Genes covered | Gene pairs | Conserved pair recall | Precision | Combined rule acc / prec / recall |
|---|---|---|---|---|---|---|
| mja, MMseqs2 | 23 | 65 / 247 | 48 | 0.234 | 1.000 | 0.904 / 0.983 / 0.891 |
| mja, 3Di | 26 | 81 / 247 | 86 | 0.281 | 1.000 | 0.904 / 0.983 / 0.891 |
| syn, MMseqs2 | 71 | 254 / 453 | 721 | 0.508 | 0.931 | 0.891 / 0.952 / 0.909 |
| syn, 3Di | 75 | 276 / 453 | 821 | 0.591 | 0.929 | 0.891 / 0.945 / 0.917 |

Distance rule alone on the same pairs: 0.880 (mja subset), 0.830 (syn subset). For reference, the full six-genome
MMseqs2 run gives conserved recall 0.719 on the mja-subset pairs: adding more genomes helps far more than
switching search mode against one genome.

Gained E. coli pairs (3Di conserved, MMseqs2 not), mja subset: 3, all gap <= 9 bp (ribosomal proteins, livHMGF,
fec). Syn subset: 13 gained (12 same TU; 11 with gap <= 50 bp: mod, moa, art, nuo, r-proteins, atp, fec, tau, a
sugar ABC transporter; the far ones are a same-TU solute-binding pair at 220 bp and a different-TU pair at 290 bp)
and 1 lost (atpD/atpG, gap 26).
New hits are remote homologs at 15-35% identity with Foldseek E-values 1e-8 to 1e-48 (ABC transporter
binding proteins and permeases, PstS, FecA/TonB receptors, NADH dehydrogenase subunit N vs NdhD, CcmK2 vs PduJ
shell proteins, Cox2 vs ubiquinol oxidase subunit 2, uL24, ATP synthase delta): `new_hits.py` output.

## Runtime

| Step | Measured | Per protein |
|---|---|---|
| ProstT5 3Di, CPU, 2 cores (247 proteins, workstation under load 3-30) | 2,182 s | 8.8 s (34 res/s) |
| ProstT5 3Di, CPU, 2 cores (453 proteins) | 4,272 s | 9.4 s (33 res/s) |
| clustersearch MMseqs2 / Foldseek, 247 proteins | 1 s / 6 s | - |

Scaling, CPU: ~18 core-s per protein, so 23,901 proteins ≈ 120 core-h, 1.5 M (Chesterfield) ≈ 7,400 core-h.
GPU: Foldseek's README claims one to two orders of magnitude faster inference; the grex job below measures it.

## For grex

`experiments/structure_mode/grex/sd6_structure_mode.sbatch` (+ `runtime.def` for the Spacedust static binary,
which needs glibc >= 2.34): full six-genome structure-mode run on one GPU, then `compare_runs.py` and
`eval_pairs.py` against the MMseqs2 run. Inputs listed in the script header (the existing `bench/sd6_mmseqs`
setDB, its clusters.tsv, ProstT5 weights, eco_pairs.tsv, the E. coli Bakta tsv). Set the GPU partition
(`CHECK_GPU_PARTITION`) and paths. Use the precompiled Foldseek GPU build (glibc >= 2.17, Turing or newer).
Expected runtime: ProstT5 probably tens of minutes on an Ampere/Ada GPU (a guess from the README speed-up, not measured: the per-protein rate it reports
sets the Chesterfield estimate), search + clustering < 15 min; CPU fallback ~120 core-hours.

## Files

- Code: `experiments/structure_mode/` (`run_subset.sh`, `subset_regions.py`, `setdb_to_fasta.py`,
  `attach_3di.py`, `compare_runs.py`, `new_hits.py`, `gained_pairs.py`, `fscpu.sh`, `fsgpu.sh`, `build_fsgpu_sm70.sh`, `grex/`).
- Data: `/data/scratch/metapanclust_data/structure_mode/sub_mja/`, `.../sub_syn/` (gff/, run/clusters_{mmseqs,3di}.tsv,
  eco_pairs_sub.tsv).
