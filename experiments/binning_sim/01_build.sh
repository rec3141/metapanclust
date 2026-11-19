#!/bin/bash
# Rebuild each species pangenome without its held-out genomes, then project genomes onto each reduced
# pangenome (gene -> family; genes in no family become singleton families = foreign):
#   proj_*     held-out genomes of both species, strict (80% id / 80% cov, PPanGGOLiN default)
#   projrel_*  same, relaxed (30% id / 50% cov): stand-in for HMM-style marker detection across species
#   projnear_* / projnearrel_*  Chesterfield bins from the same order/family ("near" contaminants)
# Also skani ANI among the held-out genomes of each species.
set -euo pipefail
# run under `nice -n 15` on the shared workstation
export PATH=/data/miniforge3/bin:$PATH
eval "$(conda shell.bash hook)"; conda activate ppanggolin
W=/data/scratch/metapanclust_data/binning_sim
SKANI=/data/dev/danaSeq/mag_analysis/conda-envs/dana-mag-derep/bin/skani
cd $W
cat held_uba9145.tsv held_mag120531.tsv > held_all.tsv
proj() {  # pangenome input outdir extra-args...
  local pan=$1 in=$2 out=$3; shift 3
  [ -d $out ] || ppanggolin projection -p $pan --fasta $in -o $out -c 1 --table -d --tmpdir $W/tmp "$@" > log_$out.txt 2>&1
}
for s in uba9145 mag120531; do
  [ -s pan_$s/pangenome.h5 ] || ppanggolin all --fasta train_$s.tsv -o pan_$s -c 1 -f -d --tmpdir $W/tmp > log_all_$s.txt 2>&1
  proj pan_$s/pangenome.h5 held_all.tsv proj_$s
  proj pan_$s/pangenome.h5 held_all.tsv projrel_$s --identity 0.3 --coverage 0.5
  proj pan_$s/pangenome.h5 held_near_$s.tsv projnear_$s
  proj pan_$s/pangenome.h5 held_near_$s.tsv projnearrel_$s --identity 0.3 --coverage 0.5
  [ -s ani_$s.tsv ] || { cut -f2 held_$s.tsv > tmp/ani_$s.txt; $SKANI triangle -l tmp/ani_$s.txt -t 1 -E 2>/dev/null \
    | awk -F'\t' 'NR>1{n=split($1,a,"/"); m=split($2,b,"/"); print a[n], b[m], $3, $4, $5}' > ani_$s.tsv; }
done
