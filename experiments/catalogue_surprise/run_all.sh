#!/bin/bash
# Level-1 catalogue surprise over the 93 Chesterfield bins. ~6 min at 2 threads, <2 GB RAM.
set -euo pipefail
export PATH=/data/miniforge3/bin:$PATH OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
cd "$(dirname "$0")"
O=/data/scratch/metapanclust_data/catalogue_surprise
mkdir -p $O
py() { nice -n 15 conda run -n metapanclust python "$@"; }
py build.py
for m in none full flat; do EM_MODE=$m py surprise.py > $O/surprise_$m.log 2>&1; done   # flat last: main fit
py analyze.py > $O/analyze.log 2>&1
py extra.py > $O/extra.log 2>&1
py extra2.py > $O/extra2.log 2>&1
py tops.py > $O/tops.log 2>&1
