#!/usr/bin/env bash
# Reproduce the pangenome survey (docs/experiments/pangenome_survey.md). 2 threads, ~1 CPU-h in total.
set -euo pipefail
X=/data/dev/metapanclust/experiments/pangenome_survey
W=/data/scratch/metapanclust_data/pangenome_survey
P=/data/scratch/metapanclust_data/pangenomes      # the two earlier species (read-only)
C=/data/scratch/metapanclust_data/chesterfield
export PATH=/data/miniforge3/bin:$PATH
mkdir -p $W && cd $W

# GTDB r226 slim metadata and our DAS Tool bins with their GTDB-Tk species
zcat /data/scratch/refdbs/GTDB/release226/metadata/bac120_metadata_r226.tsv.gz | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
  {print $h["accession"]"\t"$h["gtdb_taxonomy"]"\t"$h["checkm2_completeness"]"\t"$h["checkm2_contamination"]"\t"$h["ncbi_genome_category"]"\t"$h["genome_size"]}' > gtdb_slim.tsv
for a in mag_chesterfield mag_chesterfield_freshwater mag_chesterfield_marine; do
  awk -F'\t' -v a=$a '$1~/^dastool-/{print a"\t"$1"\t"$2}' $C/$a/taxonomy/gtdbtk/gtdbtk_taxonomy.tsv; done > our_bins.tsv
# species.tsv (>= 8 GTDB genomes, minus the two earlier species) and small_species.tsv (3-7) were picked from:
#   awk -F'\t' 'NR==FNR{n[$2]++; next}{print n[$3]+0"\t"$1"\t"$3}' gtdb_slim.tsv our_bins.tsv | sort -k1,1nr
cp $X/species.tsv $X/small_species.tsv .

for list in species.tsv small_species.tsv; do
  while IFS=$'\t' read s n; do nice -n 15 $X/build_species.sh $s "s__$n"; done < $list
done

# earlier two species: family representative proteins, from a copy of their pangenome.h5
for s in mag120531 uba9145; do
  mkdir -p prev_$s; cp $P/$s/ppanggolin/pangenome.h5 prev_$s/
  mamba run -n ppanggolin ppanggolin fasta -p prev_$s/pangenome.h5 --prot_families all -o prev_$s/fasta -f
  rm prev_$s/pangenome.h5
done
{ awk -F'\t' -v W=$W '{print $1"\t"$2"\t"W"/"$1"/ppanggolin"}' species.tsv
  printf "mag120531\tMAG-120531 sp000173115\t$P/mag120531/ppanggolin\t$W/prev_mag120531/fasta/all_protein_families.faa\n"
  printf "uba9145\tUBA9145 sp001438145\t$P/uba9145/ppanggolin\t$W/prev_uba9145/fasta/all_protein_families.faa\n"; } > survey_species.tsv
awk -F'\t' -v W=$W '{print $1"\t"$2"\t"W"/"$1"/ppanggolin"}' small_species.tsv > survey_small.tsv

for t in all:survey_species.tsv small:survey_small.tsv; do
  o=${t%%:*}; l=${t#*:}
  nice -n 15 mamba run -n ppanggolin python -I $X/survey.py $W/$o $W/$l
  nice -n 15 mamba run -n metapanclust python -I $X/coverage.py $W/$o $W/$l > $W/$o/coverage.md
  python3 -I $X/summarize.py $W/$o > $W/$o/summary.md
done

# Pfam domains for representatives of rare-present module families (TMED14 left out: contamination)
mkdir -p pfam
python3 -I $X/pfam_query.py $W > pfam/query.faa
nice -n 15 /data/dev/danaSeq/mag_analysis/conda-envs/dana-mag-annotate/bin/hmmscan --cpu 2 --cut_ga --noali \
  --domtblout pfam/query.domtbl -o /dev/null /data/scratch/refdbs/antismash_db/pfam/35.0/Pfam-A.hmm pfam/query.faa
python3 -I $X/pfam_summary.py $W > pfam/summary.md
