#!/bin/bash
# One species of the tree-vs-context run: prep -> skani -> model. Resumable: finished steps are not redone.
# usage: run_species.sh ID H5_FILE OUTDIR
# env:   RUNMODE=grex  (singularity: SIF metapanclust tools, DSIF danaSeq (resolved path), S resolved scratch)
#        RUNMODE=local (docker image metapanclust-tools:0.1; skani from the danaSeq conda env SKANI_ENV via mamba,
#                       ANI skipped if SKANI_ENV is missing)
#        SKANI_T threads for skani (default 2)
set -euo pipefail
id=$1; h5=$2; out=$3
CODE=$(cd "$(dirname "$0")" && pwd)
SKANI_T=${SKANI_T:-2}
wd=$out/work/$id
[ -e $out/res/$id.done ] && exit 0
if [ "${RUNMODE:-grex}" = grex ]; then
    py() { singularity exec -B $S $SIF env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
           python $CODE/worker.py "$@"; }
    skani_run() { singularity exec -B $S $DSIF bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                  conda activate /opt/conda/envs/dana-mag-derep && $1"; }
else
    py() { docker run --rm --cpus 1 --memory 12g --user $(id -u):$(id -g) -v /data:/data \
           -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 metapanclust-tools:0.1 python $CODE/worker.py "$@"; }
    SKANI_ENV=${SKANI_ENV:-/data/dev/metacyc_benchmark/danaSeq/mag_analysis/conda-envs/dana-mag-derep}
    skani_run() { if [ -d $SKANI_ENV ]; then PATH=/data/miniforge3/bin:$PATH mamba run -p $SKANI_ENV bash -c "$1";
                  else echo "no skani env: ANI skipped"; fi; }
fi
echo "$(date +%T) $id start"
[ -e $wd/pool.npz ] || py prep $id $h5 $out
if [ ! -s $wd/skani.tsv ] && [ -d $wd/cds ]; then
    skani_run "skani triangle -l $wd/cds_list.txt -E --min-af 0 -t $SKANI_T -o $wd/skani.tsv.part"
    [ -s $wd/skani.tsv.part ] && mv $wd/skani.tsv.part $wd/skani.tsv
    echo "$(date +%T) $id skani done"
fi
rm -rf $wd/cds          # the CDS fasta only feeds skani; keeps inode use down
py model $id $out
echo "$(date +%T) $id finished"
