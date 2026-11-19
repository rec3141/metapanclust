# Coverage test of surprising accessory modules (2026-11-18)

## Result

- **Flagellar modules: real UBA9145 content, carried by about half of the population.** The FliA / flagellar
  locus (module_6, 26 kb contig) follows the core across samples. Its depth relative to the core, after
  correcting for contig length (ratio_adj), is 0.50. All 26 persistent-family contigs of 20 kb or more sit at
  0.82 or above. Its tetranucleotide composition is typical of the core (0.34 of the core's 95th-percentile
  distance). 11 of its 26 genes, glycine cleavage gcvP/H/T among them, belong to families found in public
  UBA9145 genomes. The FlgC module (module_20, 9.6 kb) also follows the core (r 0.98, no reads outside the
  population's samples) at ratio_adj 0.64, and its composition fits the core. But 13% of core contigs of that
  length fall that low, and none of its genes are in public families. So it fits the population, but the
  coverage alone cannot tell it apart from core.
- **Most of the 36 "surprising" UBA9145 modules are contamination: 21 of 36.** 17 modules sit on contigs whose
  reads come only from samples where the UBA9145 core is absent, almost all from FAZ82916_barcode16 (the core
  has about 4x in barcode04 and 0.2x in barcode13). 4 more modules sit on contigs with no mapped reads at all.
  Only 1 of the 237 genes on the 17 off-population contigs is in a family that any public UBA9145 genome has;
  on contigs that follow the core it is 216 of 290. The marine co-assembly's UBA9145 bin contains these same
  contigs at ≥98% identity (median 99% of each contig covered). The second "surprise", that both bins carry
  the modules, is therefore one artefact repeated by two overlapping co-assemblies, not independent support.
  - Of the remaining 15 modules: 8 are core-like and 5 of them share a contig with persistent genes; 7 are
    accessory-like (ratio_adj 0.09-0.57; 4 of those have noisy profiles).
- **MAG-120531 (chesterfield_semibin_019):** its single surprising module sits on a 20.5 kb contig that also
  carries persistent genes. It is core-like (ratio_adj 0.94, r 0.81) with core composition, so it belongs to
  the population.
- **All 93 bins:** contig-level ratios can detect population accessory content only on contigs ≥10 kb.
  Below that, known members (persistent-family contigs) look "accessory" just as often as other contigs (see
  the table). Per class, as a share of the 155.8 Mb binned:

| Class (9,213 binned contigs) | Contigs | Share of bp | Meaning |
|---|---|---|---|
| core-like | 65% | 83% | ratio_adj 0.6-2 |
| accessory-like | 14.5% | 8.1% | ratio_adj < 0.6, coverage follows the bin |
| high-copy | 4.9% | 3.9% | ratio_adj > 2 |
| off-population | 11.1% | 4.2% | > 50% of depth in samples where the bin's core is absent (likely misbinned) |
| no-reads | 3.3% | 1.0% | no mapped reads in any sample |
| other-bin | 0.5% | 0.1% | another bin's core fits clearly better and composition is off |

| Contig length | Accessory-like, all bins | Accessory-like, known members |
|---|---|---|
| 3-5 kb | 22% | 46% |
| 5-10 kb | 21-23% | 21-27% |
| 10-15 kb | 18% | 7% |
| 15-20 kb | 15% | 3% |
| 20-50 kb | 7% | 0% (49 contigs) |

Only contigs ≥10 kb show more accessory-like calls than the known members do. These are 394 contigs (7.5 Mb,
6% of bp ≥10 kb, in 72 bins), with ratio_adj median 0.44 (IQR 0.32-0.53), i.e. carried by roughly 30-55% of
cells. This is the population-frequency signal, but at present depths it is a candidate list rather than a
measurement. The off-population class (1,025 contigs, 6.6 Mb) is the clearest actionable output: contigs
binned into a population whose reads come from samples where that population is absent.

## Data and method

- Main Chesterfield co-assembly, 36 nanopore samples, `jgi_summarize_bam_contig_depths`, so contig-level
  depths only. Both focal bins are found essentially in one sample: core about 4x in FAZ82916_barcode04, 34 of
  36 samples below 0.05x. 81 of 93 bins have a core profile that correlates r ≥ 0.8 with another bin's, so
  correlation across samples barely separates populations here. The discriminating evidence comes from the
  following:
  - **off**: depth in samples where the core is absent.
  - **ratio_adj**: sum of contig depth over the informative samples (core ≥ 0.2x) divided by core depth, then
    divided by the median ratio of trusted contigs of the same length class. Raw ratios fall with contig
    length: persistent contigs under 5 kb have a median of 0.72, those of 20 kb or more 1.20.
  - **tnf**: canonical tetranucleotide distance to the core centroid, in units of the core's 95th percentile.
  - **gene-family linkage**: genes on the same contig in PPanGGOLiN families that public genomes have.
- Core: contigs carrying PPanGGOLiN persistent families (focal bins). For all bins: the median profile of each
  bin's contigs of 10 kb or more. This proxy agrees with the persistent-family core (Spearman 1.0 and 0.94 on
  ratio).
- Calibration:
  - Known members are the persistent-family contigs. 93% of the UBA9145 ones pass the
    "core- or accessory-like, tnf ≤ 1.25" test, and 88% for MAG-120531.
  - Contaminant-like contigs are the contigs of bins that co-occur with the focal bin (r ≥ 0.8), scored as if
    they were members. Only 2.3% (47 of 2,072) pass for UBA9145 and 8% for MAG-120531. Coverage alone passes
    most of them; composition rejects them.
- Caveats:
  - Depth is about 4x in a single sample, so the ratio for one contig carries large counting noise.
  - Contig-level averages hide within-contig breakpoints: a chimeric junction would need per-base depth (BAMs).
  - The 21 contaminant modules all come from barcode16, but no bin there was identified as their source. By
    coverage the nearest bin is Spongiibacter (vamb_551_sub, r 0.89), by composition a Halieaceae bin.

## Files

- Script: `experiments/coverage_test/coverage_test.py` (modes `focal`, `allbins`), `marine_check.py`,
  `run.sh` (whole analysis, about 1 min, 1-2 threads).
- Outputs: `/data/scratch/metapanclust_data/coverage_test/`. The `{uba9145,mag120531}/*_{modules,contigs,genes}.tsv`
  tables are per module, per contig and per gene. `allbins/allbins_contigs.tsv` gives the class of every
  binned contig and `allbins/allbins_per_bin.tsv` the per-bin shares. Logs are `*.log`, and
  `marine_check/marine_cov98.tsv` holds the cross-assembly identity.

## Next

- Remove off-population and no-reads contigs from bins before computing surprise. For UBA9145 this alone
  removes 21 of the 36 modules from "surprising".
- Run per-base depth from the BAMs on the flagellar contigs to rule out a chimeric junction, and add more
  samples or deeper sequencing for barcode04 if the population's frequency is to be estimated.
