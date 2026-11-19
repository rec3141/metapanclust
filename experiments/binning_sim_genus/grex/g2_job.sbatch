#!/bin/bash
#SBATCH --job-name=binsim_genus
#SBATCH --account=def-rec3141
#SBATCH --partition=genoa,genlm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=700G
#SBATCH --time=2:00:00
#SBATCH --output=%x-%j.out
# One whole node, tasks run in parallel inside it (xargs -P):
#   1. pyrodigal proteins of every non-training genome           (tools container, 1 thread each)
#   2. per focal species: pangenome, mmseqs, simulation x 3 cuts  (tools container, 3 threads each)
#   3. per focal species: skani ANI vs its held-out genomes       (danaSeq container, dana-mag-derep)
#   4. analysis                                                   (tools container)
# Every step skips finished outputs, so the job can simply be resubmitted if it runs out of time.
# Submit from the login node after g1_download.sh: sbatch $BASE/code/grex/g2_job.sbatch
set -uo pipefail
export PATH=/opt/singularity-ce-4.4.1/bin:$PATH
S=$(readlink -f /home/rec3141/scratch)
B=$S/metapanclust/binning_sim_genus
C=$B/code
SIF=$S/metapanclust/sif/metapanclust-tools-0.1.sif
DANA=$(readlink -f /home/rec3141/GENICE/danaSeq/mag_analysis/.danaseq-mag-analysis.sif)
NP=${SLURM_CPUS_PER_TASK:-192}
PAR=$(((NP - 32) / 3))                # concurrent species at 3 threads each (32 cores left for skani)
run()  { singularity exec -B $S $SIF bash -c "$1"; }
runD() { singularity exec -B $S $DANA bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate /opt/conda/envs/dana-mag-derep && $1"; }
cd $B; mkdir -p prot pan hits results logs tmp ani
[ -s plan.tsv ] || { echo "plan.tsv missing: run g1_download.sh first" >&2; exit 1; }
echo "[$(date)] start: $(cut -f5 plan.tsv | sort -u | wc -l) focal species, $(cut -f1 plan.tsv | sort -u | wc -l) genomes"

# 1. proteins
awk -F'\t' '$4!="train"{print $1}' plan.tsv | sort -u > tmp/prot_todo.txt
run "cat $B/tmp/prot_todo.txt | xargs -P $NP -I{} bash $C/grex/g_prot.sh $B {}"
echo "[$(date)] proteins: $(ls prot | wc -l)"

# 2 + 3 in parallel: species pipelines (tools container) and skani (danaSeq container)
cut -f5 plan.tsv | sort -u > tmp/tags.txt
runD "cat $B/tmp/tags.txt | xargs -P 16 -I{} bash $C/grex/g_skani.sh $B {} 2" > logs/skani.log 2>&1 &
SK=$!
run "cat $B/tmp/tags.txt | xargs -P $PAR -I{} bash -c 'bash $C/grex/g_species.sh $B $C {} 3 || echo \"{} FAILED\"'" > logs/species.log 2>&1
wait $SK
echo "[$(date)] species done: $(grep -c ' done$' logs/species.log), failed: $(grep -c FAILED logs/species.log); skani: $(ls ani | wc -l)"

# 4. analysis
run "python3 $C/grex/g3_analyze.py $B" > results/summary_grex.md 2> logs/analyze.err
echo "[$(date)] analysis: $B/results/summary_grex.md"
