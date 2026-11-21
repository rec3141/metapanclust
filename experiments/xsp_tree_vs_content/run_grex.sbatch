#!/bin/bash
# Between-species tree vs gene content, grex. One whole node: prep every species (.h5 -> family proteins and
# frequencies), then all clades at once (MMseqs2 orthogroups, core AAI, leave-one-species-out prediction).
# Finished species (out/prep/<id>.json) and clades (out/res/<clade>.done) are skipped, so resubmitting resumes.
#
# layout under $BASE: code/ (this directory), sel/clades.tsv, sel/trees/<clade>.nwk, h5/<id>.h5
# submit (from $BASE):  sbatch code/run_grex.sbatch
#SBATCH --job-name=xsp
#SBATCH --account=def-rec3141
#SBATCH --partition=genoa,genlm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=400G
#SBATCH --time=03:00:00
#SBATCH --output=logs/xsp_%j.out

set -u
S=$(readlink -f /home/rec3141/scratch)
BASE=$S/metapanclust/xsp
SIF=$S/metapanclust/sif/metapanclust-tools-0.1.sif
THREADS=${THREADS:-10}
export PATH=/opt/singularity-ce-4.4.1/bin:$PATH
cd $BASE
mkdir -p out/prep out/res logs/tasks
py() { singularity exec -B $S $SIF env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       python $BASE/code/worker.py "$@"; }
export -f py
export S SIF BASE THREADS

echo "$(date) prep: $(tail -n +2 sel/clades.tsv | cut -f2 | sort -u | wc -l) species"
tail -n +2 sel/clades.tsv | cut -f2 | sort -u | xargs -P 150 -I{} bash -c \
    '[ -e out/prep/{}.json ] || py prep {} $(readlink -f h5/{}.h5) $BASE/out > logs/tasks/prep_{}.log 2>&1'
echo "$(date) prep done: $(ls out/prep/*.json | wc -l)"

run_clade() {
    cl=$1
    [ -e out/res/$cl.done ] && return 0
    ids=$(awk -F'\t' -v c=$cl '$1==c{print $2}' sel/clades.tsv | paste -sd,)
    timeout 10000 singularity exec -B $S $SIF env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
        MKL_NUM_THREADS=1 python $BASE/code/worker.py clade $cl $BASE/sel/trees/$cl.nwk $BASE/out --ids $ids --threads $THREADS \
        --clades $BASE/sel/clades.tsv > logs/tasks/clade_$cl.log 2>&1
    echo "$(date +%T) $cl exit $?"
}
export -f run_clade
tail -n +2 sel/clades.tsv | cut -f1 | sort -u | xargs -P 19 -I{} bash -c 'run_clade {}'
echo "$(date) finished; clades done: $(ls out/res/*.done 2>/dev/null | wc -l)"
