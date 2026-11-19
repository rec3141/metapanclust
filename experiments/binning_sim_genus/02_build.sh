#!/bin/bash
# For each focal species: PPanGGOLiN pangenome of its 15 training genomes (held-out genomes excluded),
# representative proteins of every family. Genes of all non-training genomes are called once with
# pyrodigal, then searched against each focal pangenome's representatives (mmseqs, 30% id / 50% cov
# on both sequences); 02_simulate.py derives strict (80/80) and relaxed (30/50) family assignments
# from the best hit. skani ANI of every non-training genome against the held-out genomes.
# Run under `nice -n 15`; 2-3 threads.
set -euo pipefail
export PATH=/data/miniforge3/bin:$PATH
eval "$(conda shell.bash hook)"; conda activate ppanggolin
W=/data/scratch/metapanclust_data/binning_sim_genus
SKANI=/data/dev/danaSeq/mag_analysis/conda-envs/dana-mag-derep/bin/skani
T=${THREADS:-2}
cd $W; mkdir -p pan prot hits tmp logs
fna() { echo $W/genomes/${1#[A-Z][A-Z]_}.fna; }
# 1. proteins of non-training genomes (header: genome|contig|start)
awk -F'\t' '$4!="train"{print $1}' genomes.tsv | while read a; do
  out=prot/$a.faa; [ -s $out ] && continue
  pyrodigal -i $(fna $a) -a tmp/$a.faa > /dev/null 2>&1
  awk -v g=$a '/^>/{split(substr($1,2),x,"_"); n=x[length(x)]; c=substr($1,2,length($1)-length(n)-2);
       print ">" g "|" c "|" $3; next} {gsub(/\*/,""); print}' tmp/$a.faa > $out; rm tmp/$a.faa
done
# 2. pangenomes + 3. searches
for tag in $(awk -F'\t' '$5!=""{print $5}' genomes.tsv | sort -u); do
  P=pan/$tag
  if [ ! -s $P/pangenome.h5 ]; then
    awk -F'\t' -v t=$tag '$5==t && $4=="train"{print $1}' genomes.tsv | while read a; do printf '%s\t%s\n' $a $(fna $a); done > tmp/train_$tag.tsv
    ppanggolin all --fasta tmp/train_$tag.tsv -o $P -c $T -f --tmpdir $W/tmp > logs/all_$tag.txt 2>&1
    ppanggolin fasta -p $P/pangenome.h5 --prot_families all -o $P/fasta -f > logs/fasta_$tag.txt 2>&1
  fi
  genus=$(awk -F'\t' -v t=$tag '$5==t{print $3; exit}' genomes.tsv)
  if [ ! -s hits/$tag.m8 ]; then
    awk -F'\t' -v g=$genus '$4!="train" && ($3==g || $4=="far"){print "prot/"$1".faa"}' genomes.tsv | xargs cat > tmp/q_$tag.faa
    reps=$(ls $P/fasta/*.faa* | head -1)
    mmseqs easy-search tmp/q_$tag.faa $reps hits/$tag.m8 tmp/mm_$tag --min-seq-id 0.3 -c 0.5 --cov-mode 0 \
      --max-accept 1 --threads $T --format-output query,target,pident,qcov,tcov,bits > logs/mm_$tag.txt 2>&1
    rm -rf tmp/mm_$tag tmp/q_$tag.faa
  fi
done
# 4. ANI: every non-training genome vs every held-out genome
if [ ! -s ani.tsv ]; then
  awk -F'\t' '$4!="train"{print $1}' genomes.tsv | while read a; do fna $a; done > tmp/q.txt
  awk -F'\t' '$4=="held"{print $1}' genomes.tsv | while read a; do fna $a; done > tmp/r.txt
  $SKANI dist --ql tmp/q.txt --rl tmp/r.txt -t $T > tmp/ani_raw.tsv 2>/dev/null
  awk -F'\t' 'NR>1{n=split($1,a,"/"); m=split($2,b,"/"); sub(/\.fna$/,"",a[n]); sub(/\.fna$/,"",b[m]); print a[n]"\t"b[m]"\t"$3"\t"$4"\t"$5}' tmp/ani_raw.tsv > ani.tsv
fi
