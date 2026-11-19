#!/usr/bin/env bash
# skani: aligned fraction of each pangenome bin (original, cleaned, removed part) against the public genomes of its
# species; reports the best ANI hit per query. Independent of the coverage rules (sequence only).
set -euo pipefail
W=/data/scratch/metapanclust_data/bin_cleaning
S=/data/scratch/metapanclust_data/pangenome_survey
SKANI=/data/dev/metacyc_benchmark/danaSeq/mag_analysis/conda-envs/dana-mag-derep/bin/skani
mkdir -p $W/af/removed; cd $W/af
tail -n +2 $W/per_bin.tsv | awk -F'\t' '$2!=""{print $1"\t"$2}' | while IFS=$'\t' read b sp; do
  d=$(awk -F'\t' -v s=$sp '$1==s{print $3}' $S/survey_species.tsv $S/survey_small.tsv)
  ls $(dirname $d)/genomes/*.fna | grep -v /chesterfield_ > refs_$sp.txt
  awk -F'\t' -v b=$b '$1==b && $4=="remove"{print $2}' $W/decisions.tsv > rm.ids
  awk 'NR==FNR{k[$1]=1;next} /^>/{p=(substr($1,2) in k)} p' rm.ids /data/scratch/metapanclust_data/chesterfield/bins/$b.fna > removed/$b.fna
  for q in /data/scratch/metapanclust_data/chesterfield/bins/$b.fna $W/cleaned/$b.fna removed/$b.fna; do
    tag=$(basename $(dirname $q)); [ -s $q ] || continue
    nice -n 15 $SKANI dist -t 2 -q $q --rl refs_$sp.txt -n 1 2>/dev/null | tail -n +2 | awk -v b=$b -v t=$tag -v s=$sp 'BEGIN{OFS="\t"}{print s,b,t,$3,$4,$5}'
  done
done > af.tsv
