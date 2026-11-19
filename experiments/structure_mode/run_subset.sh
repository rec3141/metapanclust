#!/bin/bash
# Spacedust MMseqs2 mode vs Foldseek+ProstT5 (3Di) mode on the same subset setDB, CPU, 2 threads, niced.
# usage: run_subset.sh OUTDIR name1:genome1.fna:genome1.gff3 name2:...
set -euo pipefail
X=/data/dev/metapanclust/experiments/structure_mode
W=/data/scratch/metapanclust_data/tools/prostt5/weights
out=$1; shift; mkdir -p "$out/fna" "$out/tmp"; cd "$out"
: > gff_list.txt; files=()
for pair in "$@"; do
    IFS=: read -r name fna gff <<< "$pair"
    ln -sfn "$(realpath "$fna")" "fna/$name.fna"; files+=("fna/$name.fna"); realpath "$gff" >> gff_list.txt
done
sd() { TOOL=spacedust nice -n 15 $X/fscpu.sh "$@"; }
fs() { TOOL=foldseek nice -n 15 $X/fscpu.sh "$@"; }
t() { local s=$(date +%s); "$@"; echo "TIME $(( $(date +%s) - s )) s: ${*:1:3}" >> times.txt; }
[ -f setDB.dbtype ] || t sd createsetdb "${files[@]}" setDB tmp --gff-dir gff_list.txt --gff-type CDS --threads 2 > createsetdb.log 2>&1
python3 -I $X/setdb_to_fasta.py setDB setdb.faa
echo "proteins $(grep -c '>' setdb.faa) residues $(awk '!/^>/{s+=length($0)}END{print s}' setdb.faa)" > size.txt
[ -f p5_ss.dbtype ] || t fs createdb setdb.faa p5 --prostt5-model $W --threads 2 > prostt5.log 2>&1
python3 -I $X/attach_3di.py setDB p5 > attach.log
mkdir -p tmp0 tmp2
t sd clustersearch setDB setDB clusters_mmseqs.tsv tmp0 --filter-self-match 1 --threads 2 > cs_mmseqs.log 2>&1
t sd clustersearch setDB setDB clusters_3di.tsv tmp2 --filter-self-match 1 --threads 2 --search-mode 2 > cs_3di.log 2>&1
echo done > DONE
