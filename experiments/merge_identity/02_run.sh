#!/bin/bash
# Run 01_simulate.py for every focal species, then 03_analyze.py.
set -euo pipefail
W=/data/scratch/metapanclust_data/binning_sim_genus
O=/data/scratch/metapanclust_data/merge_identity
PY=/data/scratch/grid-envs/conda/ppanggolin/bin/python
D=$(dirname "$(readlink -f "$0")")
for tag in ${TAGS:-arctica baltica distincta macleodii mediterranea mediterraneus pontiacus putrefaciens}; do
  [ -s $O/sim_$tag.tsv ] || nice -n 15 $PY -I $D/01_simulate.py $W $tag $O
done
