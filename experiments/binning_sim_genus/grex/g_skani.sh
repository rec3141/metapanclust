#!/bin/bash
# Inside the danaSeq container (env dana-mag-derep, skani 0.3.1): ANI of every non-training genome of
# one focal species against its held-out genomes.  Output ani/<tag>.tsv: query  ref  ANI  AF_q  AF_ref
# usage: g_skani.sh BASE TAG THREADS
set -euo pipefail
B=$1; tag=$2; T=$3; cd $B; mkdir -p ani tmp
[ -s ani/$tag.tsv ] && exit 0
awk -F'\t' -v t=$tag -v b=$B '$5==t && $4!="train"{a=$1; sub(/^[A-Z][A-Z]_/,"",a); print b"/genomes/"a".fna"}' plan.tsv > tmp/ani_$tag.q
awk -F'\t' -v t=$tag -v b=$B '$5==t && $4=="held"{a=$1; sub(/^[A-Z][A-Z]_/,"",a); print b"/genomes/"a".fna"}' plan.tsv > tmp/ani_$tag.r
skani dist --ql tmp/ani_$tag.q --rl tmp/ani_$tag.r -t $T 2>/dev/null \
  | awk -F'\t' 'NR>1{n=split($1,a,"/"); m=split($2,b,"/"); sub(/\.fna$/,"",a[n]); sub(/\.fna$/,"",b[m]); print b[m]"\t"a[n]"\t"$3"\t"$5"\t"$4}' > ani/$tag.tsv.part
mv ani/$tag.tsv.part ani/$tag.tsv; rm -f tmp/ani_$tag.q tmp/ani_$tag.r
