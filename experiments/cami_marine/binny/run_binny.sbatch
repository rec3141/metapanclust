#!/bin/bash
# CAMI II marine: map the 10 short-read samples to the pooled gold-standard assembly (danaSeq container:
# minimap2 -ax sr single-end, samtools; MetaBAT jgi_summarize_bam_contig_depths), then run binny (rec3141/binny
# branch modern, in the metapanclust-tools 0.2 image). Finished steps are skipped on resubmission.
#SBATCH --job-name=cami_binny
#SBATCH --account=def-rec3141
#SBATCH --partition=genoa,genlm
#SBATCH --nodes=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=700G
#SBATCH --time=24:00:00
#SBATCH --output=logs/binny_%j.out
set -euo pipefail
S=$(readlink -f /home/rec3141/scratch)
B=$S/metapanclust/cami2_marine
SIF=$S/metapanclust/sif/metapanclust-tools-0.2.sif
DSIF=$(readlink -f /home/rec3141/GENICE/danaSeq/mag_analysis/.danaseq-mag-analysis.sif)
export PATH=/opt/singularity-ce-4.4.1/bin:$PATH
cd $B; mkdir -p map binny logs
dana() { env=$1; shift; singularity exec -B $S $DSIF bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate /opt/conda/envs/$env && $*"; }
ASM=$B/anonymous_gsa_pooled.fasta.gz
[ -s map/asm.mmi ] || dana dana-mag-assembly "minimap2 -x sr -d map/asm.mmi $ASM"
echo "$(date) index done"
for i in $(seq 0 9); do
  [ -s map/sample_$i.bam ] && continue
  dana dana-mag-assembly "minimap2 -ax sr -t 160 map/asm.mmi reads/sample_$i.fq.gz | samtools sort -@ 24 -m 4G -o map/sample_$i.bam.tmp - && mv map/sample_$i.bam.tmp map/sample_$i.bam"
  echo "$(date) sample $i mapped"
done
[ -s map/depth.txt ] || dana dana-mag-binning "jgi_summarize_bam_contig_depths --outputDepth map/depth.txt map/sample_*.bam"
echo "$(date) depth: $(wc -l < map/depth.txt) contigs"
[ -d binny_db/hmms ] || singularity exec -B $S $SIF binny prepare-db $S/databases/checkm_data $B/binny_db
singularity exec -B $S $SIF env OMP_NUM_THREADS=1 binny run --assembly $ASM --depth map/depth.txt --db $B/binny_db \
  --outdir $B/binny --threads 96 --sample cami2_marine
echo "$(date) binny done: $(ls binny/bins | wc -l) bins"
