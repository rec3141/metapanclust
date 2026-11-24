#!/bin/bash
# Run worker.py for every selected species once its PanGBank .h5 has been downloaded (6 at a time, local docker).
O=${1:-/data/scratch/metapanclust_data/bgc_atlas/level1}
mkdir -p $O/res $O/logs
one() {
  id=$1; O=$2
  [ -s $O/res/$id.targets.tsv ] && return 0
  for i in $(seq 1 720); do [ -s $O/h5/$id.h5 ] && break; sleep 10; done
  [ -s $O/h5/$id.h5 ] || { echo "$id: no h5"; return 0; }
  docker run --rm --cpus 2 --user $(id -u):$(id -g) -v /data:/data -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
    metapanclust-tools:0.1 python /data/dev/metapanclust/experiments/bgc_atlas/level1/worker.py $id $O/h5/$id.h5 \
    $O/atlas/$id.tsv $O/res > $O/logs/$id.log 2>&1
  echo "$(date +%T) $id exit $?"
}
export -f one
tail -n +2 $O/species.tsv | cut -f4 | xargs -P 6 -I{} bash -c "one {} $O"
echo ALL_SPECIES_DONE
