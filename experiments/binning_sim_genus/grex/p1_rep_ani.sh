#!/bin/bash
# LOCAL (workstation, run once): skani ANI between GTDB r226 species representatives within each genus,
# for genera with a focal-capable species (>= 21 HQ isolate genomes) and >= 2 species with >= 2 HQ
# isolates. HQ isolate: CheckM2 >= 95% complete, <= 2% contamination, <= 200 contigs, NCBI category none.
# Output OUT/rep_ani.tsv: genus  focal_species  other_species  ANI  AF_focal  AF_other
# usage (under nice -n 15): p1_rep_ani.sh OUTDIR [THREADS]
set -euo pipefail
O=$1; T=${2:-3}
META=/data/scratch/refdbs/GTDB/release226/metadata/bac120_metadata_r226.tsv.gz
DB=/data/scratch/refdbs/GTDB/release226/skani
SKANI=/data/dev/danaSeq/mag_analysis/conda-envs/dana-mag-derep/bin/skani
mkdir -p $O/lists; cd $O
[ -s hq_iso.tsv ] || zcat $META | awk -F'\t' 'NR>1 && $3>=95 && $4<=2 && $14<=200 && $59=="none" {print $1"\t"$18"\t"$20"\t"$14}' > hq_iso.tsv
# species table: genus species n_hq_isolates representative
awk -F'\t' '{split($3,t,";"); n[t[7]]++; g[t[7]]=t[6]; r[t[7]]=$2} END{for(s in n) print g[s]"\t"s"\t"n[s]"\t"r[s]}' hq_iso.tsv | sort > sp.tsv
awk '{print $1"\t"$2}' $DB/genome_paths.tsv > paths.tsv     # file  dir (space-separated in GTDB)
: > rep_ani.tsv
for g in $(awk -F'\t' '$3>=2{c[$1]++} $3>=21{f[$1]++} END{for(g in f) if(c[g]>=2) print g}' sp.tsv | sort); do
  awk -F'\t' -v g=$g 'NR==FNR{p[$1]=$2; next} $1==g && $3>=2 {a=substr($4,4); f=a"_genomic.fna.gz"; if(f in p) print $2"\t"($3>=21)"\t""'$DB'/"p[f]f}' paths.tsv sp.tsv > lists/$g.tsv
  awk -F'\t' '$2==1{print $3}' lists/$g.tsv > lists/$g.q; cut -f3 lists/$g.tsv > lists/$g.r
  [ -s lists/$g.q ] || continue
  $SKANI dist --ql lists/$g.q --rl lists/$g.r -t $T 2>/dev/null | awk -F'\t' -v g=$g 'NR==FNR{s[$3]=$1; next} FNR>1 && $1!=$2 {print g"\t"s[$2]"\t"s[$1]"\t"$3"\t"$5"\t"$4}' lists/$g.tsv - >> rep_ani.tsv
done
