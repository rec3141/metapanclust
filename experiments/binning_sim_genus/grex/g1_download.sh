#!/bin/bash
# LOGIN NODE (compute nodes may lack internet): download every genome in plan_full.tsv from NCBI Datasets
# (GENOME_FASTA), sequentially in batches of 100, into BASE/genomes/<accession>.fna.
# Safe to re-run: skips genomes already present; prints the accessions still missing at the end.
# usage: g1_download.sh   (BASE resolved from the scratch symlink)
set -uo pipefail
B=${BASE:-$(readlink -f /home/rec3141/scratch)/metapanclust/binning_sim_genus}
cd $B; mkdir -p genomes dl
cut -f1 plan_full.tsv | sed 's/^[A-Z][A-Z]_//' | sort -u | while read a; do [ -s genomes/$a.fna ] || echo $a; done > dl/todo.txt
echo "$(wc -l < dl/todo.txt) genomes to download"
split -l 100 -d -a 3 dl/todo.txt dl/batch_
for b in dl/batch_[0-9]*; do
  d=$b.d; rm -rf $d; mkdir -p $d
  for try in 1 2 3; do
    curl -sS --fail --retry 3 -o $d/g.zip \
      "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$(paste -sd, $b)/download?include_annotation_type=GENOME_FASTA" \
      && unzip -oq $d/g.zip -d $d && break
    echo "batch $b attempt $try failed" >&2; sleep 30
  done
  for a in $(cat $b); do
    f=$(ls $d/ncbi_dataset/data/$a/*.fna 2>/dev/null | head -1)
    [ -n "$f" ] && mv "$f" genomes/$a.fna
  done
  rm -rf $d $b
  echo "$b done: $(ls genomes | wc -l) genomes"
  sleep 2                                   # NCBI Datasets asks for <= 5 requests/s without an API key
done
cut -f1 plan_full.tsv | sed 's/^[A-Z][A-Z]_//' | sort -u | while read a; do [ -s genomes/$a.fna ] || echo $a; done > dl/missing.txt
echo "missing after download: $(wc -l < dl/missing.txt) (listed in $B/dl/missing.txt)"
# plan.tsv = plan_full.tsv without rows whose genome could not be downloaded (used by the job)
awk -F'\t' 'NR==FNR{m[$1]; next} {a=$1; sub(/^[A-Z][A-Z]_/,"",a); if (!(a in m)) print}' <(echo __none__; cat dl/missing.txt) plan_full.tsv > plan.tsv
echo "plan.tsv: $(wc -l < plan.tsv) of $(wc -l < plan_full.tsv) rows kept"
