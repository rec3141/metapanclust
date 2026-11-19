#!/bin/bash
# CPU Foldseek/Spacedust static builds in the Ubuntu 24.04 runtime image, capped at 2 CPUs and 8 GB.
tool=${TOOL:-foldseek}
exec docker run --rm -i --cpus 2 --cpu-shares 64 --memory 8g --user $(id -u):$(id -g) -v /data:/data -w "$PWD" metapanclust-runtime:24.04 \
    /data/scratch/metapanclust_data/tools/spacedust/bin/$tool "$@"
