#!/bin/bash
# CAMI II marine, round 2:
#  (a) species references with more than 400 genomes are rebuilt from a random 400 (max_genomes), and the bins that
#      timed out against them are rescored (4 h limit);
#  (b) genus references (genus/refs.tsv: full genus, or genus minus the bin's own species) are built with ref-clade;
#  (c) bins in genus/bin_refs.tsv are scored against their genus reference (tree only, no context model; 4 h limit).
# Resubmitting skips finished references and bins.
#SBATCH --job-name=cami_genus
#SBATCH --account=def-rec3141
#SBATCH --partition=genoa,genlm
#SBATCH --nodes=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=700G
#SBATCH --time=24:00:00
#SBATCH --output=logs/genus_%j.out
set -u
S=$(readlink -f /home/rec3141/scratch)
BASE=$S/metapanclust/cami2_marine
SIF=$S/metapanclust/sif/metapanclust-tools-0.1.sif
export PATH=/opt/singularity-ce-4.4.1/bin:$PATH
cd $BASE; mkdir -p logs/tasks refs_genus out_genus
mpc() { singularity exec -B $S $SIF env PYTHONPATH=$BASE/code OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
        MKL_NUM_THREADS=1 python -m metapanclust "$@"; }
mkdb() { singularity exec -B $S $SIF env PYTHONPATH=$BASE/code python -c \
         "from metapanclust.annotate import ref_db; ref_db('$1')"; }
export -f mpc mkdb; export S SIF BASE

# (a) large species references
for d in refs/*/; do
  id=$(basename $d)
  n=$(python3 -c "import json;print(json.load(open('$d/ref.json')).get('genomes',0))" 2>/dev/null || echo 0)
  if [ "$n" -gt 400 ]; then echo "rebuild species ref $id ($n genomes)"; rm -rf $d; fi
done
tail -n +2 sp/genome_species.tsv | awk -F'\t' '$5!=""{print $5}' | sort -u | xargs -P 40 -I{} bash -c \
  '[ -e refs/{}/mmdb/ref.dbtype ] || { mpc ref $(readlink -f h5/{}.h5) $BASE/refs/{} --max-genomes 400 && mkdb $BASE/refs/{}; } \
   > logs/tasks/ref_{}.log 2>&1'
echo "$(date) species refs: $(ls -d refs/*/mmdb | wc -l)"

# (b) genus references
tail -n +2 genus/refs.tsv | xargs -P 40 -L 1 bash -c \
  '[ -e refs_genus/$0/mmdb/ref.dbtype ] && exit 0; h=""; for i in ${1//,/ }; do h="$h $(readlink -f h5/$i.h5)"; done; \
   { mpc ref-clade $BASE/refs_genus/$0 $h --names $2 --threads 4 && mkdb $BASE/refs_genus/$0; } > logs/tasks/gref_$0.log 2>&1; \
   echo "$(date +%T) genus ref $0 exit $?"'
echo "$(date) genus refs: $(ls -d refs_genus/*/mmdb | wc -l)"

# (c) scoring: species-level bins still missing, then genus bins
{ awk '{print $1, "refs/"$2, "out", "ctx"}' logs/score_tasks.txt
  tail -n +2 genus/bin_refs.tsv | awk '{print $1, "refs_genus/"$2, "out_genus", "noctx"}'; } > logs/genus_tasks.txt
echo "$(date) scoring tasks: $(wc -l < logs/genus_tasks.txt)"
xargs -P 180 -L 1 bash -c '[ -e $2/$0.summary.tsv ] && exit 0; [ -e $1/mmdb/ref.dbtype ] || exit 0; \
  opt=""; [ "$3" = noctx ] && opt="--no-context"; \
  timeout 14400 bash -c "mpc score $BASE/$1 $BASE/bins/$0.fna $BASE/$2/$0 --threads 1 $opt" > logs/tasks/$2_$0.log 2>&1; \
  rm -rf $2/$0.work' < logs/genus_tasks.txt
echo "$(date) scored: species $(ls out/*.summary.tsv | wc -l), genus $(ls out_genus/*.summary.tsv 2>/dev/null | wc -l)"
