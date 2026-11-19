#!/bin/bash
# Simulate and score bins for both species (foreign contaminant = the other species, and related
# Chesterfield bins), then summarise. ~1 CPU-minute after 01_build.sh.
set -euo pipefail
W=/data/scratch/metapanclust_data/binning_sim
PY=/data/scratch/grid-envs/conda/ppanggolin/bin/python
D=$(dirname "$(readlink -f "$0")")
cd $W; mkdir -p results
nice -n 15 $PY $D/02_simulate.py $W uba9145 mag120531 results/sim_uba9145_mag120531.tsv
nice -n 15 $PY $D/02_simulate.py $W mag120531 uba9145 results/sim_mag120531_uba9145.tsv
nice -n 15 $PY $D/02_simulate.py $W uba9145 near_uba9145 results/sim_uba9145_near.tsv
nice -n 15 $PY $D/02_simulate.py $W mag120531 near_mag120531 results/sim_mag120531_near.tsv
cp ani_*.tsv results/
nice -n 15 $PY $D/03_analyze.py results results/sim_*.tsv > results/summary_stdout.txt
