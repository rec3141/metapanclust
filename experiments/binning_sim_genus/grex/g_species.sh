#!/bin/bash
# Inside the tools container: one focal species end to end.
#   pangenome of its training genomes (ppanggolin all), family representatives, mmseqs search of
#   its held-out and contaminant genes against them (30% id / 50% cov, best hit), then simulated
#   bins at three strict family cuts (80, 90, 95% identity).
# usage: g_species.sh BASE CODE_DIR TAG THREADS
set -euo pipefail
B=$1; C=$2; tag=$3; T=$4
cd $B; P=pan/$tag; mkdir -p pan hits results logs tmp/$tag
fna() { echo $B/genomes/${1#[A-Z][A-Z]_}.fna; }
if [ ! -s $P/fasta.done ]; then
  awk -F'\t' -v t=$tag '$5==t && $4=="train"{print $1}' plan.tsv | while read a; do printf '%s\t%s\n' $a $(fna $a); done > tmp/$tag/train.tsv
  ppanggolin all --fasta tmp/$tag/train.tsv -o $P -c $T -f --tmpdir $B/tmp/$tag > logs/all_$tag.txt 2>&1
  ppanggolin fasta -p $P/pangenome.h5 --prot_families all -o $P/fasta -f > logs/fasta_$tag.txt 2>&1
  touch $P/fasta.done
fi
if [ ! -s hits/$tag.m8 ]; then
  awk -F'\t' -v t=$tag '$5==t && $4!="train"{print "prot/"$1".faa"}' plan.tsv | xargs cat > tmp/$tag/q.faa
  reps=$(ls $P/fasta/*.faa* | head -1)
  mmseqs easy-search tmp/$tag/q.faa $reps hits/$tag.m8.part tmp/$tag/mm --min-seq-id 0.3 -c 0.5 --cov-mode 0 \
    --max-accept 1 --threads $T --split-memory-limit 8G --format-output query,target,pident,qcov,tcov,bits > logs/mm_$tag.txt 2>&1
  mv hits/$tag.m8.part hits/$tag.m8
fi
for id in 0.8 0.9 0.95; do
  out=results/sim_${tag}_id${id}.tsv
  [ -s $out ] || GENOMES=$B/plan.tsv python3 $C/03_simulate.py $B $tag $out $id >> logs/sim_$tag.txt 2>&1
done
rm -rf tmp/$tag
# keep only what the analysis needs from the pangenome (inodes: ppanggolin writes ~100 files)
find $P -mindepth 1 -maxdepth 1 ! -name gene_families.tsv ! -name partitions ! -name pangenome.h5 ! -name fasta.done -exec rm -rf {} +
echo "$tag done"
