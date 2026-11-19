#!/bin/bash
# GPU Foldseek build in the Ubuntu 24.04 runtime image; data paths mounted.
exec docker run --rm -i --gpus all --user $(id -u):$(id -g) -v /data:/data -w "$PWD" metapanclust-runtime:24.04 \
    /data/scratch/metapanclust_data/tools/fsgpu/foldseek/bin/foldseek "$@"
