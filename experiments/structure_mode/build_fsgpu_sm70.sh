#!/bin/bash
# Build Foldseek (same commit as the shipped binaries) with CUDA for Volta (sm_70): the precompiled GPU
# build targets Turing and newer only, and its ProstT5 kernels fail on the GV100.
# Mirrors the upstream Dockerfile GPU recipe; build tools live only in the throwaway container.
set -euo pipefail
B=/data/scratch/metapanclust_data/structure_mode/build
exec docker run --rm --cpus 2 --memory 8g -v $B:/b debian:bookworm-slim bash -c '
set -e
apt-get update -qq && apt-get install -y -qq build-essential curl xxd git cmake zlib1g-dev libbz2-dev libatomic1 ca-certificates >/dev/null
curl -sO https://developer.download.nvidia.com/compute/cuda/repos/debian12/x86_64/cuda-keyring_1.1-1_all.deb && dpkg -i cuda-keyring_1.1-1_all.deb >/dev/null
apt-get update -qq && apt-get install -y -qq cuda-nvcc-12-6 cuda-cudart-dev-12-6 libcublas-dev-12-6 ninja-build >/dev/null
curl --proto =https --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- --profile minimal -y >/dev/null 2>&1
export PATH=/root/.cargo/bin:$PATH CUDACXX=/usr/local/cuda/bin/nvcc
mkdir -p /b/build_sm70 && cd /b/build_sm70
LIBGOMP=/usr/lib/gcc/x86_64-linux-gnu/12/
cmake -GNinja -DHAVE_AVX2=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DFORCE_STATIC_DEPS=1 \
  -DENABLE_CUDA=1 -DCMAKE_CUDA_ARCHITECTURES="70-real" \
  -DOpenMP_C_FLAGS="-fopenmp -I${LIBGOMP}" -DOpenMP_C_LIB_NAMES=gomp -DOpenMP_CXX_FLAGS="-fopenmp -I${LIBGOMP}" -DOpenMP_CXX_LIB_NAMES=gomp -DOpenMP_gomp_LIBRARY=${LIBGOMP}/libgomp.a \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. /b/src
cmake --build . -j2
chown -R '"$(id -u):$(id -g)"' /b/build_sm70
'
