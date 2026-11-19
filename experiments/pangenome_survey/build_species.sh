#!/usr/bin/env bash
# Build one species pangenome: GTDB r226 genomes of the species (NCBI download) + our Chesterfield
# DAS Tool bins of that species, then `ppanggolin all`.
# usage: build_species.sh SHORTNAME "GTDB species name (s__...)"
set -euo pipefail
short=$1; sp=$2
W=/data/scratch/metapanclust_data/pangenome_survey
C=/data/scratch/metapanclust_data/chesterfield
export PATH=/data/miniforge3/bin:$PATH
d=$W/$short; mkdir -p $d/genomes $d/dl
cd $d
# public genomes
awk -F'\t' -v s="$sp" '$2 ~ (";" s "$") {print substr($1,4)}' $W/gtdb_slim.tsv | sort -u > accessions.txt
if [ ! -s dl/g.zip ]; then
  curl -sS -o dl/g.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$(paste -sd, accessions.txt)/download?include_annotation_type=GENOME_FASTA"
  (cd dl && unzip -oq g.zip)
fi
for a in $(cat accessions.txt); do
  f=$(ls dl/ncbi_dataset/data/$a/*.fna 2>/dev/null | head -1 || true)
  [ -n "$f" ] && cp -n "$f" genomes/$a.fna || echo "missing $a" >&2
done
# our bins
for asm in mag_chesterfield mag_chesterfield_freshwater mag_chesterfield_marine; do
  tag=${asm#mag_}
  awk -F'\t' -v a=$asm -v s="$sp" '$1==a && $3 ~ (";" s "$") {print $2}' $W/our_bins.tsv | while read b; do
    out=genomes/${tag}_${b#dastool-}.fna
    [ -s $out ] && continue
    awk -F'\t' -v b=$b '$2==b{print $1}' $C/$asm/binning/dastool/contig2bin.tsv > /tmp/pgs_$$.ids
    awk 'NR==FNR{k[$1]=1;next} /^>/{p=(substr($1,2) in k)} p' /tmp/pgs_$$.ids $C/$asm/binning/dastool/allbins.fa > $out
    rm -f /tmp/pgs_$$.ids
  done
done
for f in genomes/*.fna; do echo -e "$(basename $f .fna)\t$PWD/$f"; done > genomes.tsv
if [ ! -s ppanggolin/modules_in_genomes.tsv ]; then
  rm -rf ppanggolin
  nice -n 15 mamba run -n ppanggolin ppanggolin all --fasta genomes.tsv -o ppanggolin -c ${THREADS:-2} > ppanggolin.log 2>&1
  nice -n 15 mamba run -n ppanggolin ppanggolin fasta -p ppanggolin/pangenome.h5 --prot_families all -o ppanggolin/fasta -f >> ppanggolin.log 2>&1 || true
fi
echo "$short: $(wc -l < accessions.txt) accessions, $(ls genomes | wc -l) genomes"
