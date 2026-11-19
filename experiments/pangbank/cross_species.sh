#!/bin/bash
# Project Chesterfield bins of species absent from PanGBank onto a congeneric PanGBank pangenome
# (GTDB_all 2.1.0, r232), at the default 80% identity and at 50%, to see whether a sister species can
# stand in as the species prior.
set -euo pipefail
export PATH=/data/miniforge3/bin:$PATH
eval "$(conda shell.bash hook)"; conda activate ppanggolin
D=/data/scratch/metapanclust_data/pangbank; S=/data/scratch/metapanclust_data/pangenome_survey
mkdir -p $D/xsp
# our species : PanGBank sister pangenome id : tag
for x in tmed14:17084:tmed14_sp003507115 amylibacter:13175:amylibacter_sp000153745 planktomarina:16002:planktomarina_temperata; do
  IFS=: read sp id tag <<< "$x"
  [ -s $D/h5/$tag.h5 ] || nice -n 15 curl -s -m 600 -o $D/h5/$tag.h5 "https://pangbank-api.genoscope.cns.fr/pangenomes/$id/file"
  for f in $S/$sp/genomes/chesterfield_*.fna; do printf '%s\t%s\n' $(basename $f .fna) $f; done > $D/bins/$sp.tsv
  for idt in 0.8 0.5; do
    o=$D/xsp/${sp}_on_${tag}_id$idt
    [ -d $o ] || /usr/bin/time -f "$sp on $tag id $idt: wall %e s" nice -n 15 ppanggolin projection \
        -p $D/h5/$tag.h5 --fasta $D/bins/$sp.tsv -o $o --identity $idt -c 3 --log $o.log 2>> $D/xsp/time.log
  done
done
grep -H -E "Genes:|Completeness|specific families" $D/xsp/*/*/projection_summary.yaml
