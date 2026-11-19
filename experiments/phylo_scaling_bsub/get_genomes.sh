#!/bin/bash
# Genome FASTA + NCBI GFF3 for every panel entry: link from the phylo_scaling genome dir when present, else download (NCBI datasets v2 API).
# usage: get_genomes.sh PANEL_TSV OLD_GENOME_DIR OUT_GENOME_DIR DL_DIR
set -euo pipefail
panel=$1 old=$2 out=$3 dl=$4; mkdir -p "$out" "$dl"
tail -n+2 "$panel" | while IFS=$'\t' read -r n acc lvl grp org; do
    [ -s "$out/$n.fna" ] && [ -s "$out/$n.gff3" ] && continue
    if [ -s "$old/$n.fna" ]; then ln -sf "$(realpath "$old/$n.fna")" "$out/$n.fna"; ln -sf "$(realpath "$old/$n.gff3")" "$out/$n.gff3"; continue; fi
    z=$dl/$acc.zip
    [ -s "$z" ] || curl -sfL -m 600 -o "$z" "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF"
    d=$dl/$acc; rm -rf "$d"; mkdir -p "$d"; unzip -q -o "$z" -d "$d"
    cat "$d"/ncbi_dataset/data/$acc/*.fna > "$out/$n.fna"; cp "$d"/ncbi_dataset/data/$acc/genomic.gff "$out/$n.gff3"
    echo "$n $acc $(grep -m1 '>' "$out/$n.fna" | cut -c1-80) CDS=$(grep -cP '\tCDS\t' "$out/$n.gff3")"
    rm -rf "$d"
done
