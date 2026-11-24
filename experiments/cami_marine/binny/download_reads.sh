#!/bin/bash
# CAMI II marine short reads (10 samples, ~5.5 GB each) -> reads/sample_<i>.fq.gz. Run on a grex LOGIN node
# (compute nodes have no internet). Resumable: finished samples are skipped.
set -euo pipefail
B=${1:?base dir}
URL=https://frl.publisso.de/data/frl:6425521/marine/short_read
mkdir -p $B/reads $B/reads/tmp
for i in $(seq 0 9); do
  out=$B/reads/sample_$i.fq.gz
  [ -s $out ] && { echo "sample $i present"; continue; }
  rm -rf $B/reads/tmp/$i && mkdir -p $B/reads/tmp/$i
  curl -sSfL --retry 5 "$URL/marmgCAMI2_sample_${i}_reads.tar.gz" | tar xz -C $B/reads/tmp/$i --wildcards '*.fq.gz'
  f=$(find $B/reads/tmp/$i -name '*.fq.gz' | head -1)
  [ -n "$f" ] && mv "$f" $out && rm -rf $B/reads/tmp/$i
  echo "$(date +%T) sample $i: $(du -h $out | cut -f1)"
done
echo DOWNLOAD_DONE
