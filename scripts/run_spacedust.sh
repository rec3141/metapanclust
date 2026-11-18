#!/bin/bash
# All-against-all Spacedust on genomes given as FASTA + Bakta GFF3 pairs.
# usage: run_spacedust.sh OUTDIR SEARCH_MODE THREADS name1:genome1.fna:genome1.gff3 [name2:...]
#   SEARCH_MODE 0 = MMseqs2 sequence search, 2 = Foldseek with ProstT5 3Di (needs the ProstT5 weights)
set -euo pipefail
out=$1 mode=$2 threads=$3; shift 3
mkdir -p "$out/fna" "$out/tmp"
: > "$out/gff_list.txt"; files=()
for pair in "$@"; do
    IFS=: read -r name fna gff <<< "$pair"
    ln -sfn "$(realpath "$fna")" "$out/fna/$name.fna"; files+=("$out/fna/$name.fna"); realpath "$gff" >> "$out/gff_list.txt"
done
spacedust createsetdb "${files[@]}" "$out/setDB" "$out/tmp" --gff-dir "$out/gff_list.txt" --gff-type CDS --threads "$threads" > "$out/createsetdb.log" 2>&1 || { echo "createsetdb failed: $out/createsetdb.log" >&2; exit 1; }
extra=(); [ "$mode" != 0 ] && extra=(--search-mode "$mode")
spacedust clustersearch "$out/setDB" "$out/setDB" "$out/clusters.tsv" "$out/tmp" --filter-self-match 1 --threads "$threads" "${extra[@]}" > "$out/clustersearch.log" 2>&1 || { echo "clustersearch failed: $out/clustersearch.log" >&2; exit 1; }
echo "clusters: $(grep -c '^#' "$out/clusters.tsv")  hits: $(grep -c '^>' "$out/clusters.tsv")"
