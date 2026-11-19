#!/bin/bash
# Simulate and score bins for every focal species, then summarise. A few CPU-seconds per species
# after 02_build.sh.
set -euo pipefail
W=/data/scratch/metapanclust_data/binning_sim_genus
PY=/data/scratch/grid-envs/conda/ppanggolin/bin/python
D=$(dirname "$(readlink -f "$0")")
cd $W; mkdir -p results
for tag in $(awk -F'\t' '$5!=""{print $5}' genomes.tsv | sort -u); do
  [ -s hits/$tag.m8 ] && nice -n 15 $PY -I $D/03_simulate.py $W $tag results/sim_$tag.tsv
done
nice -n 15 $PY -I $D/04_analyze.py $W results > results/summary.md
nice -n 15 $PY -I $D/06_plot.py results $D/fig
nice -n 15 $PY -I $D/07_threshold.py results > results/calibrated_threshold.md
