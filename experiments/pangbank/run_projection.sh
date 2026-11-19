#!/bin/bash
# Export module tables from each PanGBank pangenome and project the Chesterfield bins onto it.
# usage: run_projection.sh   (reads targets.tsv: short<TAB>h5 tag<TAB>bins list tsv)
set -euo pipefail
export PATH=/data/miniforge3/bin:$PATH
eval "$(conda shell.bash hook)"; conda activate ppanggolin
D=/data/scratch/metapanclust_data/pangbank
S=/data/scratch/metapanclust_data/pangenome_survey
P=/data/scratch/metapanclust_data/pangenomes
mkdir -p $D/bins $D/tables $D/proj
# bin FASTAs, as used by the survey builds
ls $S/pseudothioglobus/genomes/chesterfield_*.fna; ls $S/patiriisocius/genomes/chesterfield_*.fna
mk() { for f in "${@:2}"; do b=$(basename $f .fna); printf '%s\t%s\n' $b $f; done > $D/bins/$1.tsv; }
mk pseudothioglobus $S/pseudothioglobus/genomes/chesterfield_*.fna
mk patiriisocius $S/patiriisocius/genomes/chesterfield_*.fna
mk uba9145 $P/uba9145/genomes/chesterfield_*.fna
mk mag120531 $P/mag120531/genomes/chesterfield_*.fna
for t in mag120531_r226 patiriisocius_r226 uba9145_r226 uba9145_r232 pseudothioglobus_r226 pseudothioglobus_r232; do
  sp=${t%_r*}
  [ -s $D/tables/$t/modules_in_genomes.tsv ] || nice -n 15 ppanggolin write_pangenome -p $D/h5/$t.h5 \
      --stats --modules --partitions -o $D/tables/$t -f --log $D/tables/$t.log
  [ -s $D/tables/$t/rtab/gene_presence_absence.Rtab ] || nice -n 15 ppanggolin write_pangenome -p $D/h5/$t.h5 \
      --Rtab -o $D/tables/$t/rtab -f --log $D/tables/$t.rtab.log
  [ -s $D/tables/$t/fam_fasta/all_protein_families.faa ] || nice -n 15 ppanggolin fasta -p $D/h5/$t.h5 \
      --prot_families all -o $D/tables/$t/fam_fasta -f --log $D/tables/$t.fasta.log
  if [ ! -d $D/proj/$t ]; then
    /usr/bin/time -f "$t projection wall %e s, cpu %U+%S s, maxrss %M kB" \
      nice -n 15 ppanggolin projection -p $D/h5/$t.h5 --fasta $D/bins/$sp.tsv -o $D/proj/$t \
      --table -c 3 --log $D/proj/$t.log 2>> $D/proj/time.log
  fi
done
