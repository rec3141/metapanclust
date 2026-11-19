#!/bin/bash
# Full six-genome Spacedust structure mode (Foldseek + ProstT5 3Di) vs the existing MMseqs2-mode run.
# Inputs (copy from the workstation, ~45 MB + 2.4 GB weights):
#   $IN/setDB*            /data/scratch/metapanclust_data/bench/sd6_mmseqs/setDB* (23,901 proteins, keeps gene ids)
#   $IN/clusters_mmseqs.tsv  the MMseqs2-mode result, same dir, clusters.tsv
#   $IN/prostt5/prostt5-f16.gguf, $IN/eco_pairs.tsv, $IN/ecoli_bakta.tsv
#   $X = this experiments/structure_mode directory (setdb_to_fasta.py, attach_3di.py, compare_runs.py) (run from the repo checkout so ../../scripts/eval_pairs.py resolves)
# Tools: foldseek GPU build (https://mmseqs.com/foldseek/foldseek-linux-gpu.tar.gz, glibc >= 2.17, runs natively;
#   needs Turing or newer: the workstation GV100 is not supported), spacedust static (mmseqs.com/spacedust) in
#   runtime.sif (apptainer build runtime.sif runtime.def).
# Expected: ProstT5 on one Ampere/Ada GPU ~10-30 min for 7.45 M residues; search + clustering < 15 min.
#   CPU fallback (drop --gpu): ~20 res/s per core measured on the workstation -> ~100 core-hours.
#SBATCH --job-name=sd6_3di
#SBATCH --partition=CHECK_GPU_PARTITION
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
set -euo pipefail
IN=${IN:?set IN}; X=${X:?set X}; OUT=${OUT:-$IN/sd6_3di}; FSGPU=${FSGPU:?path to GPU foldseek}; SD=${SD:?path to spacedust}
SIF=${SIF:?runtime.sif}; FSCPU=${FSCPU:?path to static CPU foldseek next to spacedust}
mkdir -p "$OUT/tmp" && cd "$OUT"
cp -a "$IN"/setDB* .
python3 -I "$X/setdb_to_fasta.py" setDB setdb.faa
s=$(date +%s); "$FSGPU" createdb setdb.faa p5 --prostt5-model "$IN/prostt5" --gpu 1 --threads 8 > prostt5.log 2>&1
echo "prostt5 $(( $(date +%s) - s )) s for $(grep -c '>' setdb.faa) proteins" | tee times.txt
python3 -I "$X/attach_3di.py" setDB p5
s=$(date +%s); apptainer exec -B "$OUT" -B "$(dirname "$SD")" "$SIF" "$SD" clustersearch setDB setDB clusters_3di.tsv tmp \
    --filter-self-match 1 --threads 8 --search-mode 2 --foldseek-path "$FSCPU" > cs_3di.log 2>&1
echo "clustersearch $(( $(date +%s) - s )) s" | tee -a times.txt
python3 -I "$X/compare_runs.py" setDB.lookup ecoli,bsub,pae,mtb,syn,mja mmseqs="$IN/clusters_mmseqs.tsv" 3di=clusters_3di.tsv > compare.txt
for f in "$IN/clusters_mmseqs.tsv" clusters_3di.tsv; do python3 -I "$X/../../scripts/eval_pairs.py" "$IN/eco_pairs.tsv" "$f" "$IN/ecoli_bakta.tsv"; done > eval.txt
