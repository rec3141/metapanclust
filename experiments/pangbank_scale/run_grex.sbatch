#!/bin/bash
# PanGBank data-size experiment (A: context model vs genomes per species; B: genus shrinkage), grex.
# One array element = one full node; it runs its share of tasks.tsv (every NCHUNK-th line, starting at its
# index; tasks.tsv is sorted most expensive first) with up to $PAR single-threaded workers at once.
# Finished tasks (A/<id>.done, B/<genus>.done) are skipped, so resubmitting the same command resumes.
#
# submit (from $BASE):  sbatch --array=0-1 code/run_grex.sbatch
#   NCHUNK must equal the array size (default 2; override with --export=ALL,NCHUNK=3 and --array=0-2).
#SBATCH --job-name=pbscale
#SBATCH --account=def-rec3141
#SBATCH --partition=genoa,genlm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=600G
#SBATCH --time=02:00:00
#SBATCH --output=logs/pbscale_%A_%a.out

set -u
BASE=/home/rec3141/scratch/metapanclust/pangbank_scale
SIF=/home/rec3141/scratch/metapanclust/sif/metapanclust-tools-0.1.sif
NCHUNK=${NCHUNK:-2}
PAR=${PAR:-180}                 # concurrent workers per node (B tasks use 2 threads for MMseqs2 briefly)
TASK_TIMEOUT=${TASK_TIMEOUT:-6900}   # seconds; A workers write results after every subsample size
export PATH=/opt/singularity-ce-4.4.1/bin:$PATH
cd $BASE
mkdir -p out/A out/B out/work logs/tasks
CHUNK=${SLURM_ARRAY_TASK_ID:-0}

awk -v n=$NCHUNK -v c=$CHUNK '(NR-1) % n == c' sel/tasks.tsv > logs/tasks/chunk_${CHUNK}.tsv
echo "$(date) node $(hostname) chunk $CHUNK/$NCHUNK: $(wc -l < logs/tasks/chunk_${CHUNK}.tsv) tasks, $PAR parallel"

run_one() {
    mode=$1; shift
    if [ "$mode" = species ]; then
        id=$1; tag=A_$id
        [ -e out/A/$id.done ] && return 0
        args="species $id h5 out --threads 1"
    else
        genus=$1; ids=$2; tag=B_$genus
        [ -e out/B/$genus.done ] && return 0
        args="genus $genus $ids h5 out --threads 2"
    fi
    timeout $TASK_TIMEOUT singularity exec -B /home/rec3141/scratch $SIF \
        env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        python code/worker.py $args > logs/tasks/$tag.log 2>&1
    echo "$(date +%T) $tag exit $?"
}
export -f run_one
export SIF TASK_TIMEOUT
tr '\t' ' ' < logs/tasks/chunk_${CHUNK}.tsv | xargs -P $PAR -L 1 bash -c 'run_one "$@"' _
echo "$(date) chunk $CHUNK finished; done markers: A $(ls out/A/*.done 2>/dev/null | wc -l), B $(ls out/B/*.done 2>/dev/null | wc -l)"
