#!/bin/bash
# metapanclust v0 held-out evaluation on grex: build references for 9 focal/foreign species pairs (PanGBank .h5
# from pangbank_scale), then one process per held-out genome (20 per pair; 'eval --shard i/20').
# Resubmitting skips finished shards (bins.tsv present).
#SBATCH --job-name=mpc_eval
#SBATCH --account=def-rec3141
#SBATCH --partition=genoa,genlm
#SBATCH --nodes=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=500G
#SBATCH --time=06:00:00
#SBATCH --output=logs/eval_%j.out
set -u
S=$(readlink -f /home/rec3141/scratch)
BASE=$S/metapanclust/scorer_eval
H5=$S/metapanclust/pangbank_scale/h5
SIF=$S/metapanclust/sif/metapanclust-tools-0.1.sif
export PATH=/opt/singularity-ce-4.4.1/bin:$PATH
cd $BASE; mkdir -p refs out logs/tasks
mpc() { singularity exec -B $S $SIF env PYTHONPATH=$BASE/code OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
        MKL_NUM_THREADS=1 python -m metapanclust "$@"; }
export -f mpc; export S SIF BASE H5
tail -n +2 code/experiments/scorer_eval/pairs.tsv | cut -f2,3 | tr '\t' '\n' | sort -u | xargs -P 18 -I{} bash -c \
  '[ -e refs/{}/mmdb/ref.dbtype ] || { mpc ref $H5/{}.h5 $BASE/refs/{} && singularity exec -B $S $SIF env PYTHONPATH=$BASE/code python -c "from metapanclust.annotate import ref_db; ref_db(\"$BASE/refs/{}\")"; } > logs/tasks/ref_{}.log 2>&1'
echo "$(date) refs: $(ls -d refs/*/mmdb | wc -l)"
tail -n +2 code/experiments/scorer_eval/pairs.tsv | while IFS=$'\t' read g f x; do
  for i in $(seq 0 19); do echo "$g $f $x $i"; done; done | xargs -P 180 -L 1 bash -c \
  '[ -e out/$0/shard_$3/bins.tsv ] && exit 0; mpc eval $BASE/refs/$1 $BASE/refs/$2 $BASE/out/$0/shard_$3 --n 20 --shard $3/20 --threads 1 > logs/tasks/eval_$0_$3.log 2>&1; echo "$(date +%T) $0 $3 exit $?"'
echo "$(date) finished: $(ls out/*/shard_*/bins.tsv 2>/dev/null | wc -l) shards"
