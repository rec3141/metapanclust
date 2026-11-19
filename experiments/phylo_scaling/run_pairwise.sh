#!/bin/bash
# Spacedust (MMseqs2 mode) E. coli K-12 query set against each panel genome separately, so any panel is a union of pairwise runs.
# usage: run_pairwise.sh PANEL_TSV GENOME_DIR OUTDIR THREADS
set -euo pipefail
panel=$1 gdir=$(realpath "$2") out=$3 threads=$4
export PATH=$(dirname "$(realpath "$0")")/bin:$PATH
eco_fna=/data/dev/metacyc_benchmark/ecoli/ecoli_k12.fasta
eco_gff=/data/dev/metacyc_benchmark/ecoli/mag_out/annotation/bakta/basic/annotation.gff3
mkdir -p "$out/query/tmp"; out=$(realpath "$out")
if [ ! -s "$out/query/setDB.dbtype" ]; then
    echo "$eco_gff" > "$out/query/gff_list.txt"
    spacedust createsetdb "$eco_fna" "$out/query/setDB" "$out/query/tmp" --gff-dir "$out/query/gff_list.txt" --gff-type CDS --threads "$threads" > "$out/query/createsetdb.log" 2>&1
fi
for n in $(tail -n+2 "$panel" | cut -f1); do
    d=$out/$n; [ -s "$d/clusters.tsv" ] && continue
    mkdir -p "$d/tmp"; echo "$gdir/$n.gff3" > "$d/gff_list.txt"
    t0=$(date +%s)
    spacedust createsetdb "$gdir/$n.fna" "$d/setDB" "$d/tmp" --gff-dir "$d/gff_list.txt" --gff-type CDS --threads "$threads" < /dev/null > "$d/createsetdb.log" 2>&1 || { echo "$n createsetdb failed"; continue; }
    spacedust clustersearch "$out/query/setDB" "$d/setDB" "$d/clusters.tsv.part" "$d/tmp" --threads "$threads" < /dev/null > "$d/clustersearch.log" 2>&1 || { echo "$n clustersearch failed"; continue; }
    mv "$d/clusters.tsv.part" "$d/clusters.tsv"; rm -rf "$d/tmp" "$d"/clusters.tsv.part_* "$d"/clusters.tsv_*
    echo "$n $(( $(date +%s) - t0 ))s clusters $(grep -c '^#' "$d/clusters.tsv") hits $(grep -c '^>' "$d/clusters.tsv")"
done
