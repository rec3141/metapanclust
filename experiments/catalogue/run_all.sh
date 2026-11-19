#!/usr/bin/env bash
# Rebuild the Chesterfield conserved-gene-cluster catalogue and the E. coli operon check.
# Needs: mmseqs (conda env metapanclust), a python with pandas/scipy/igraph/leidenalg (VENV below).
# Runtime on 2 threads: gene families ~10 min (skipped when fam_cluster.tsv exists), the rest ~3 min.
set -euo pipefail
D=/data/scratch/metapanclust_data
OUT=$D/catalogue
E=$(cd "$(dirname "$0")" && pwd)
VENV=$OUT/venv          # python -m venv --system-site-packages venv && venv/bin/pip install igraph leidenalg
PY="nice -n 15 $VENV/bin/python"
export PATH=/data/scratch/grid-envs/conda/metapanclust/bin:$PATH
CH=$D/chesterfield
cd "$OUT"; mkdir -p final bench

# bin -> GTDB classification (genus used for link support)
ls $CH/bins/*.gff3 | xargs -n1 basename | sed 's/\.gff3$//' > bins.txt
grep -wFf bins.txt $CH/mag_chesterfield/taxonomy/gtdbtk/gtdbtk_taxonomy.tsv | cut -f1,2 > bin_taxonomy.tsv

# 1. gene table + MMseqs2 gene families (30% id, 80% coverage) over binned proteins only
$PY $E/01_genes_families.py $CH/bins $CH/mag_chesterfield/annotation/bakta/basic/annotation.faa fam30 0.3 0.8 2
# 2. loci: links supported by >= 1 other genus, split at support dips < 0.34 x flanking max, bridge <= 1 gene
$PY $E/02_loci.py $CH/sd_bins_mmseqs/clusters.tsv fam30/genes.tsv bin_taxonomy.tsv final/chesterfield 1 0.34 1
# 3. families: matched-gene similarity >= 0.4, >= 2 matched genes, best-match filter 0.9, Leiden CPM 0.1 / 0.4
$PY $E/03_families.py fam30/genes.tsv final/chesterfield_loci.tsv bin_taxonomy.tsv $CH/sd_bins_mmseqs/clusters.tsv \
    final/chesterfield 0.4 0.1 2 1 0.9 0.4
# 4. sanity checks on known systems
$PY $E/04_sanity.py final/chesterfield_members.tsv final/chesterfield_catalogue.tsv fam30/genes.tsv > final/chesterfield_sanity.tsv

# 5. E. coli transcription units: locus boundaries on the six-genome benchmark
mkdir -p bench/gff
for g in $(cat $D/bench/sd6_mmseqs/gff_list.txt); do n=$(echo $g | cut -d/ -f5); ln -sf $g bench/gff/$n.gff3; done
[ -s bench/all.faa ] || for g in $(cat $D/bench/sd6_mmseqs/gff_list.txt); do cat $(dirname $g)/annotation.faa; done > bench/all.faa
for g in bsub ecoli mja mtb pae syn; do printf '%s\td__X;p__%s;c__;o__;f__;g__%s;s__\n' $g $g $g; done > bench/tax.tsv
$PY $E/01_genes_families.py bench/gff bench/all.faa bench/fam30 0.3 0.8 2
for gap in 0 1 2; do
  $PY $E/02_loci.py $D/bench/sd6_mmseqs/clusters.tsv bench/fam30/genes.tsv bench/tax.tsv bench/l_g$gap 1 0.34 $gap
  echo "max_gap=$gap"; $PY $E/05_eval_operons.py $D/bench/eco_pairs.tsv bench/fam30/genes.tsv bench/l_g${gap}_loci.tsv
done > bench/eval.txt
