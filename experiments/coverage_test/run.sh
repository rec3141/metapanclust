#!/bin/bash
# Coverage test of surprising modules and all-bin contig classes (Chesterfield main co-assembly).
# ~1 min, 1 thread, < 2 GB RAM.
set -euo pipefail
export PATH=/data/miniforge3/bin:$PATH
H=$(dirname "$(readlink -f "$0")")
S=/data/scratch/metapanclust_data/chesterfield
P=/data/scratch/metapanclust_data/pangenomes
O=/data/scratch/metapanclust_data/coverage_test
A="$S/assembly_chesterfield/mapping/depths.txt $S/mag_chesterfield/binning/dastool/contig2bin.tsv $S/mag_chesterfield/binning/dastool/allbins.fa"
B="$S/mag_chesterfield/annotation/bakta/basic/annotation.tsv $S/mag_chesterfield/mge/genomad"
R="nice -n 15 mamba run -n metapanclust"
$R python $H/coverage_test.py focal $A $B dastool-maxbin_040 $P/uba9145/ppanggolin chesterfield_maxbin_040 $O/uba9145 > $O/uba9145.log
$R python $H/coverage_test.py focal $A $B dastool-semibin_019 $P/mag120531/ppanggolin chesterfield_semibin_019 $O/mag120531 > $O/mag120531.log
$R python $H/coverage_test.py allbins $A $O/allbins $O/uba9145/chesterfield_maxbin_040_contigs.tsv \
    $O/mag120531/chesterfield_semibin_019_contigs.tsv > $O/allbins.log
# same contigs in the marine co-assembly's UBA9145 bin?
M=$O/marine_check; mkdir -p $M
$R python - "$S" "$O" <<'PY'
import sys, pandas as pd
S, O = sys.argv[1:3]
main = set(pd.read_csv(f'{O}/uba9145/chesterfield_maxbin_040_contigs.tsv', sep='\t', index_col=0).index)
c2b = pd.read_csv(f'{S}/mag_chesterfield_marine/binning/dastool/contig2bin.tsv', sep='\t', header=None, index_col=0)[1]
mar = set(c2b[c2b == 'dastool-maxbin_019'].index)
for fa, keep, out in ((f'{S}/mag_chesterfield/binning/dastool/allbins.fa', main, 'main.fa'),
                      (f'{S}/mag_chesterfield_marine/binning/dastool/allbins.fa', mar, 'marine.fa')):
    w = False
    with open(f'{O}/marine_check/{out}', 'w') as fo:
        for l in open(fa):
            if l[0] == '>': w = l[1:].split()[0] in keep
            if w: fo.write(l)
PY
$R mmseqs easy-search $M/main.fa $M/marine.fa $M/hits.m8 $M/tmp --search-type 3 --threads 2 -v 1 > /dev/null
rm -rf $M/tmp
$R python $H/marine_check.py $O/uba9145/chesterfield_maxbin_040_contigs.tsv $M/hits.m8 $M/marine_cov98.tsv > $O/marine_check.log
