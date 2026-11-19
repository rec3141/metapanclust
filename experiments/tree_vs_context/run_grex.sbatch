#!/bin/bash
# Tree vs context (relatedness vs gene-gene context), grex. One array element = one whole node; it runs its share of
# sel/tasks_tvc.tsv (every NCHUNK-th species, starting at its index; the list is sorted most expensive first) with
# up to $PAR species at once. Per species: prep (metapanclust container) -> skani triangle on the pool's CDS
# (danaSeq container, env dana-mag-derep) -> model (metapanclust container); see code/run_species.sh. Finished species (out/res/<id>.done)
# are skipped and finished steps (work/<id>/{pool.npz,skani.tsv}) are not redone, so resubmitting resumes.
#
# submit (from $BASE):  sbatch --array=0-1 code/run_grex.sbatch
#   NCHUNK must equal the array size (default 2; e.g. --array=0-2 --export=ALL,NCHUNK=3).
#SBATCH --job-name=tvc
#SBATCH --account=def-rec3141
#SBATCH --partition=genoa,genlm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=600G
#SBATCH --time=03:00:00
#SBATCH --output=logs/tvc_%A_%a.out

set -u
S=$(readlink -f /home/rec3141/scratch)
BASE=$S/metapanclust/tree_vs_context
H5=$S/metapanclust/pangbank_scale/h5
SIF=$S/metapanclust/sif/metapanclust-tools-0.1.sif
DSIF=$(readlink -f /home/rec3141/GENICE/danaSeq/mag_analysis/.danaseq-mag-analysis.sif)

NCHUNK=${NCHUNK:-2}
PAR=${PAR:-100}                       # species at once per node (100 = one wave for 200 species on 2 nodes)
SKANI_T=${SKANI_T:-2}
TASK_TIMEOUT=${TASK_TIMEOUT:-9900}    # seconds per species; the model step writes results after every subsample
export PATH=/opt/singularity-ce-4.4.1/bin:$PATH
cd $BASE
mkdir -p out/res out/work logs/tasks
CHUNK=${SLURM_ARRAY_TASK_ID:-0}
awk -v n=$NCHUNK -v c=$CHUNK '(NR-1) % n == c' sel/tasks_tvc.tsv > logs/tasks/chunk_${CHUNK}.tsv
echo "$(date) node $(hostname) chunk $CHUNK/$NCHUNK: $(wc -l < logs/tasks/chunk_${CHUNK}.tsv) species, $PAR parallel"
echo "metapanclust sif: $SIF"; echo "danaSeq sif: $DSIF"

run_one() {
    id=$1
    [ -e out/res/$id.done ] && return 0
    timeout $TASK_TIMEOUT bash code/run_species.sh $id $H5/$id.h5 out > logs/tasks/$id.log 2>&1
    echo "$(date +%T) $id exit $?"
}
export -f run_one
export S SIF DSIF H5 SKANI_T TASK_TIMEOUT RUNMODE=grex
xargs -P $PAR -I{} bash -c 'run_one {}' < logs/tasks/chunk_${CHUNK}.tsv
echo "$(date) chunk $CHUNK finished; done: $(ls out/res/*.done 2>/dev/null | wc -l)"
