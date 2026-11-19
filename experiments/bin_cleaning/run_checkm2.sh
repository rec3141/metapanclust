#!/usr/bin/env bash
# CheckM2 on original and cleaned versions of the changed bins (one call, 3 threads).
set -euo pipefail
W=/data/scratch/metapanclust_data/bin_cleaning
ENV=/data/dev/metacyc_benchmark/danaSeq/mag_analysis/conda-envs/dana-mag-quality
export PATH=/data/miniforge3/bin:$PATH
nice -n 15 mamba run -p $ENV checkm2 predict -i $W/ck_in -x fna -o $W/checkm2 -t 3 --force \
  --database_path /data/scratch/refdbs/checkm2/CheckM2_database/uniref100.KO.1.dmnd --tmpdir $W/tmp
