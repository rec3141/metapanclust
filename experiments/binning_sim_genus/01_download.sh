#!/bin/bash
# Download the genomes in genomes.tsv from NCBI Datasets (GENOME_FASTA), in batches of 50,
# into W/genomes/<accession>.fna. GTDB accessions carry an RS_/GB_ prefix that NCBI does not use.
set -euo pipefail
W=/data/scratch/metapanclust_data/binning_sim_genus
cd $W; mkdir -p genomes dl
cut -f1 genomes.tsv | sed 's/^[A-Z][A-Z]_//' > accessions.txt
split -l 50 -d accessions.txt dl/batch_
for b in dl/batch_*; do
  d=$b.d; [ -d $d ] && continue
  mkdir -p $d
  curl -sS -o $d/g.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$(paste -sd, $b)/download?include_annotation_type=GENOME_FASTA"
  (cd $d && unzip -oq g.zip)
done
for a in $(cat accessions.txt); do
  f=$(ls dl/batch_*.d/ncbi_dataset/data/$a/*.fna 2>/dev/null | head -1 || true)
  [ -n "$f" ] && cp -n "$f" genomes/$a.fna || echo "missing $a" >&2
done
rm -rf dl
