#!/bin/bash
# Inside the tools container: pyrodigal proteins of one genome, header genome|contig|start.
# usage: g_prot.sh BASE ACCESSION(with RS_/GB_ prefix)
set -euo pipefail
B=$1; a=$2; out=$B/prot/$a.faa
[ -s $out ] && exit 0
fna=$B/genomes/${a#[A-Z][A-Z]_}.fna
[ -s $fna ] || { echo "missing genome $a" >&2; exit 0; }
tmp=$B/tmp/prot_$a.faa
pyrodigal -i $fna -a $tmp > /dev/null 2>&1
awk -v g=$a '/^>/{split(substr($1,2),x,"_"); n=x[length(x)]; c=substr($1,2,length($1)-length(n)-2);
     print ">" g "|" c "|" $3; next} {gsub(/\*/,""); print}' $tmp > $out.part && mv $out.part $out
rm -f $tmp
