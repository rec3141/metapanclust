#!/bin/bash
# CAMI II marine: score every submitted bin with metapanclust v0 against the PanGBank pangenome of its majority
# source genome's GTDB species (oracle species assignment), and run CheckM2 on every bin for comparison.
# Inputs in $BASE: anonymous_gsa_pooled.fasta.gz, members.tsv.gz, bins.tsv, sp/genome_species.tsv, h5/<id>.h5,
# code/metapanclust (package). Finished steps are skipped, so resubmitting resumes.
#SBATCH --job-name=cami_score
#SBATCH --account=def-rec3141
#SBATCH --partition=genoa,genlm
#SBATCH --nodes=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=700G
#SBATCH --time=08:00:00
#SBATCH --output=logs/score_%j.out
set -u
S=$(readlink -f /home/rec3141/scratch)
BASE=$S/metapanclust/cami2_marine
SIF=$S/metapanclust/sif/metapanclust-tools-0.1.sif
DSIF=$(readlink -f /home/rec3141/GENICE/danaSeq/mag_analysis/.danaseq-mag-analysis.sif)
CK2DB=$(ls $S/databases/checkm2_db/*/*.dmnd 2>/dev/null | head -1)
export PATH=/opt/singularity-ce-4.4.1/bin:$PATH
cd $BASE; mkdir -p logs/tasks refs out
mpc() { singularity exec -B $S $SIF env PYTHONPATH=$BASE/code OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
        MKL_NUM_THREADS=1 python -m metapanclust "$@"; }
export -f mpc; export S SIF BASE

[ -e bins/.done ] || { singularity exec -B $S $SIF python write_bins.py anonymous_gsa_pooled.fasta.gz members.tsv.gz bins \
                       && touch bins/.done; }
echo "$(date) bins: $(ls bins | wc -l)"

# CheckM2 on all bins (64 threads, in the background)
if [ ! -s checkm2_out/quality_report.tsv ]; then
  singularity exec -B $S $DSIF bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate /opt/conda/envs/dana-mag-quality && checkm2 predict -i bins -x fna -o checkm2_out -t 64 \
    --database_path $CK2DB --force" > logs/checkm2.log 2>&1 &
  CK=$!
fi

# references for every species with a PanGBank pangenome
tail -n +2 sp/genome_species.tsv | awk -F'\t' '$5!=""{print $5}' | sort -u | xargs -P 60 -I{} bash -c \
  '[ -e refs/{}/mmdb/ref.dbtype ] || { mpc ref $(readlink -f h5/{}.h5) $BASE/refs/{} && \
   singularity exec -B $S $SIF env PYTHONPATH=$BASE/code python -c "from metapanclust.annotate import ref_db; ref_db(\"$BASE/refs/{}\")"; } \
   > logs/tasks/ref_{}.log 2>&1'
echo "$(date) refs: $(ls -d refs/*/mmdb 2>/dev/null | wc -l)"

# score bins whose majority genome has a reference
awk -F'\t' 'NR==FNR{if(FNR>1 && $5!="") sp[$1]=$5; next} FNR>1 && ($5 in sp){print $1, sp[$5]}' \
  sp/genome_species.tsv bins.tsv > logs/score_tasks.txt
echo "$(date) scoring $(wc -l < logs/score_tasks.txt) bins"
xargs -P 120 -L 1 bash -c '[ -e out/$0.summary.tsv ] && exit 0; [ -e refs/$1/mmdb/ref.dbtype ] || exit 0; \
  timeout 3600 mpc score $BASE/refs/$1 $BASE/bins/$0.fna $BASE/out/$0 --threads 1 > logs/tasks/score_$0.log 2>&1; \
  rm -rf out/$0.work' < logs/score_tasks.txt
echo "$(date) scored: $(ls out/*.summary.tsv 2>/dev/null | wc -l)"
[ -n "${CK:-}" ] && wait $CK
echo "$(date) checkm2: $(wc -l < checkm2_out/quality_report.tsv 2>/dev/null)"
