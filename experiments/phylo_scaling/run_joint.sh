#!/bin/bash
# Control: one Spacedust run of the E. coli query set against several panel genomes in a single target set DB,
# to check that a union of pairwise runs matches a joint run.
# usage: run_joint.sh GENOME_DIR QUERY_SETDB OUTDIR THREADS name1 [name2 ...]
set -euo pipefail
gdir=$(realpath "$1") q=$2 out=$3 threads=$4; shift 4
export PATH=$(dirname "$(realpath "$0")")/bin:$PATH
mkdir -p "$out/tmp"; out=$(realpath "$out"); : > "$out/gff_list.txt"; files=()
for n in "$@"; do files+=("$gdir/$n.fna"); echo "$gdir/$n.gff3" >> "$out/gff_list.txt"; done
spacedust createsetdb "${files[@]}" "$out/setDB" "$out/tmp" --gff-dir "$out/gff_list.txt" --gff-type CDS --threads "$threads" < /dev/null > "$out/createsetdb.log" 2>&1
spacedust clustersearch "$q" "$out/setDB" "$out/clusters.tsv" "$out/tmp" --threads "$threads" < /dev/null > "$out/clustersearch.log" 2>&1
rm -rf "$out/tmp" "$out"/clusters.tsv_*
echo "clusters $(grep -c '^#' "$out/clusters.tsv") hits $(grep -c '^>' "$out/clusters.tsv")"
