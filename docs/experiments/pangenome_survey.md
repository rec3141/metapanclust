# Module surprise across Chesterfield species (2026-11-18)

This extends the two-species module-surprise prototype (`docs/first_tests.md`) to every Chesterfield species
with at least 3 GTDB r226 genomes: 13 species with 8-67 public genomes, plus 6 with 3-6.

**Result.** Accessory modules found in our bins and in few or none of the public genomes are common. In
16 of the 19 bins that are at least 50% complete, at least one module occurs in <= 10% of the public
genomes. Most of the surprise, though, comes from **two contaminated bins**:

- **TMED14 sp001438205**: 68 such modules. They carry ribosomal proteins, gyrase, ATP synthase and IF-2,
  so the bin holds a second genome.
- **UBA9145 sp001438145**: 39 such modules, including the flagellar sets reported in
  `first_tests.md`. Most look like a co-binned second organism, not the population's own content.

Two checks separate contamination from genuine accessory content without any curated database:

- **Anchoring:** does the module's contig also carry persistent genes?
- **Coverage:** does the module follow the bin's core across the 36 samples?

In the other species the modules are almost all anchored (47 of 50) and track the core coverage. They are
mostly unannotated genes. Where function can be assigned, they point to:

- carbohydrate and algal-polysaccharide uptake and use: alginate/uronate, mannitol, maltose, glycerol-3-phosphate;
- proteases and Bacteroidota T9SS-sorted secreted proteins;
- regulators and two-component systems;
- a few IS elements and integrases.

They do **not** point to CRISPR or phage. Only 3 of 157 modules sit on geNomad plasmid/virus contigs.

## Method

Code is in `experiments/pangenome_survey/` and data in `/data/scratch/metapanclust_data/pangenome_survey/`.
`run_all.sh` reruns everything at 2 threads, about 1 CPU-hour in total.

1. **Species choice.** These are the GTDB-Tk species of DAS Tool bins (`dastool-*`) in the three
   Chesterfield co-assemblies, counted against GTDB r226 `bac120_metadata`.
   - 13 species have >= 8 GTDB genomes. They are listed in the table below and include the two earlier
     species, which were reused read-only from `metapanclust_data/pangenomes/`.
   - 6 species have 3-6 GTDB genomes.
   - No freshwater-assembly species has more than 2 GTDB genomes, so the freshwater extension was not
     possible at species level.
2. **Pangenomes** (`build_species.sh`). Steps:
   - download the public genomes from NCBI Datasets (`GENOME_FASTA`);
   - add our bins of the species from all three assemblies, as contigs taken from `allbins.fa` via
     `contig2bin.tsv`;
   - run `ppanggolin all -c 2`, then `ppanggolin fasta --prot_families all`.
   
   Each species takes 1-3 minutes.
3. **Surprise** (`survey.py`). This is the same model as `scripts/module_surprise.py`:
   - a module counts as present at completion >= 0.5;
   - prevalence is p = (k+1)/(n+2) over the n public genomes;
   - surprise is -log2 p when the module is present and -log2(1 - p c) when absent, where c is the bin's
     persistent-family completeness.
   
   A module is **rare-present** in a bin if the bin carries it and k/n <= 0.1.
4. **Function.** Each family's representative protein was searched with MMseqs2 against the Bakta proteins
   of the main co-assembly. For genes in main-assembly bins, the Bakta product was taken by exact
   coordinates instead. Bakta calls 75% of all CDS in this assembly "hypothetical protein", so I also ran
   Pfam-A 35 domains (`hmmscan --cut_ga`) on the 553 non-TMED14 rare-module family representatives
   (`pfam_query.py`, `pfam_summary.py`). Products were sorted into coarse keyword classes (`classes.py`).
5. **MGE and anchoring.**
   - MGE: whether the gene's contig appears in the geNomad plasmid/virus/provirus summaries. These exist
     only for the main assembly. For marine-bin genes, the contig of a >= 95% identity Bakta hit is used
     as a proxy, which rarely applies (`n/a`).
   - Anchoring: whether any contig carrying the module also carries a persistent-partition gene.
6. **Coverage** (`coverage.py`). This uses main-assembly bins only, with per-sample depths from
   `assembly_chesterfield/mapping/depths.txt` (36 samples). Steps:
   - The core profile is the length-weighted mean depth of the bin's contigs that carry >= 3 persistent
     genes.
   - For each module I compute Pearson r of log2(depth + 0.1) against the core, and the ratio of summed
     depth (module / core).
   - Calibration: each core contig is compared against the core built from the remaining contigs.

## Species

| Species | Public genomes (median compl.) | Our bins | Families | Pers./shell/cloud | Modules |
|---|---|---|---|---|---|
| MAG-120531 sp000173115 | 67 (0.94) | 2 | 4394 | 1472/267/2655 | 41 |
| UBA9145 sp001438145 | 52 (0.91) | 2 | 5878 | 1818/408/3652 | 105 |
| Pseudothioglobus sp905182815 | 23 (0.81) | 2 | 4011 | 1441/260/2310 | 44 |
| Patiriisocius sp905182755 | 17 (0.90) | 1 | 3149 | 1632/127/1390 | 31 |
| Amylibacter sp905182285 | 13 (0.95) | 2 | 3027 | 1814/223/990 | 29 |
| Planktomarina forsetii | 12 (0.89) | 1 | 4001 | 1908/371/1722 | 71 |
| UBA9320 sp003484585 | 12 (0.87) | 2 | 4383 | 1890/198/2295 | 45 |
| TMED14 sp001438205 | 12 (0.98) | 2 | 4110 | 1536/496/2078 | 76 |
| Pontimonas sp017857245 | 12 (0.96) | 2 | 2470 | 1210/162/1098 | 41 |
| UBA10364 sp905182675 | 11 (0.97) | 3 | 2398 | 1514/173/711 | 23 |
| Marimicrobium sp913061115 | 10 (0.84) | 2 | 5134 | 2164/625/2345 | 88 |
| BACL14 sp913058475 | 8 (0.83) | 2 | 1917 | 1071/196/650 | 26 |
| CACCCN01 sp905182145 | 8 (0.83) | 1 | 1721 | 1089/98/534 | 24 |
| *3-6 genomes:* Thioglobus_A sp028352435 | 6 (0.76) | 2 | 1719 | 1213/97/409 | 17 |
| Hel1-33-131 sp018608625 | 5 (0.93) | 2 | 2554 | 1548/96/910 | 29 |
| MAG-121220-bin8 sp905182765 | 5 (0.88) | 2 | 1926 | 1293/334/299 | 35 |
| Oceaniferula sp037620805 | 5 (0.96) | 1 | 4761 | 3033/0/1728 | 30 |
| Aquiluna sp905182225 | 4 (0.96) | 2 | 1512 | 1121/151/240 | 12 |
| Pseudothioglobus sp913057955 | 3 (0.94) | 2 | 2566 | 1171/675/720 | 45 |

The bins come from overlapping co-assemblies (main and marine), so a species' two bins are mostly the same
population assembled twice. "Carried by both bins" is therefore not independent evidence.

## Surprise per bin

The table below covers the 13 species with >= 8 genomes. Rare-present means k/n <= 0.1; "none" counts
modules found in no public genome. Bits are summed over all of the species' modules. GTDB AF is the
GTDB-Tk alignment fraction to the species representative; CheckM2 is shown as completeness/contamination
(main assembly only).

| Species | Bin | Compl. | Rare-present (none) | Bits total / present / absent | GTDB AF | CheckM2 |
|---|---|---|---|---|---|---|
| tmed14 | maxbin_102 | 0.70 | 68 (68) | 263 / 261 / 2 | 0.68 | 84 / **17.6** |
| tmed14 | marine_maxbin_051 | 0.73 | 68 (68) | 264 / 261 / 2 | 0.71 | |
| uba9145 | maxbin_040 | 0.92 | 39 (36) | 246 / 231 / 14 | 0.76 | 95 / 4.2 |
| uba9145 | marine_maxbin_019 | 0.95 | 39 (36) | 246 / 232 / 14 | 0.78 | |
| marimicrobium | lorbin_000 | 0.92 | 15 (8) | 117 / 104 / 13 | 0.89 | 92 / 2.0 |
| marimicrobium | marine_lorbin_001 | 0.93 | 15 (8) | 120 / 106 / 14 | 0.88 | |
| pontimonas | semibin_040_sub | 0.59 | 14 (14) | 65 / 56 / 9 | 0.85 | 49 / 0.9 |
| pontimonas | marine_semibin_022 | 0.70 | 14 (14) | 69 / 58 / 11 | 0.82 | |
| uba9320 | maxbin_050 | 0.92 | 8 (8) | 60 / 47 / 13 | 0.64 | 97 / 4.5 |
| uba9320 | marine_semibin_014 | 0.88 | 8 (8) | 57 / 45 / 12 | 0.60 | |
| pseudothioglobus | lorbin_018 / marine_lorbin_007 | 0.74 / 0.72 | 3 (1) | 37-38 | 0.88 / 0.91 | 77 / 0.1 |
| amylibacter | lorbin_025 / marine_semibin_004 | 0.81 / 0.99 | 2 (2) | 30-32 | 0.88 / 0.97 | 94 / 0.8 |
| uba10364 | marine_semibin_002 (vamb_171, vamb_193: 0.32, 0.11) | 0.96 | 4 (2) | 27 | 0.93 | |
| mag120531 | semibin_019 (marine_lorbin_024: 0.30) | 0.84 | 2 (1) | 24 | 0.76 | 73 / 6.2 |
| bacl14 | marine_semibin_010 (lorbin_058: 0.18) | 0.79 | 0 | 23 | 0.98 | |
| cacccn01 | marine_semibin_016 | 0.77 | 0 | 22 | 0.77 | |
| planktomarina | lorbin_071 | 0.28 | 0 | 21 | 0.88 | 27 / 0.0 |
| patiriisocius | marine_semibin_023 | 0.83 | 0 | 21 | 0.74 | |

Per bin, over the 15 bins that are >= 50% complete outside TMED14 and UBA9145:

- median surprise is **37 bits**: 28 bits from present modules and 10 from absent ones;
- the median bin carries **3 rare-present modules**;
- 12 of the 15 carry at least one.

With TMED14 and UBA9145 included, the median rises to 57 bits. In the 3-6-genome tier the median is
30 bits and 4 rare-present modules per bin, with 10 of 11 bins carrying at least one. Absence contributes
little everywhere, 2-16 bits per bin.

Bits are not comparable across species. A module present in our bin and in none of n public genomes costs
log2(n+2) bits: 5.75 for UBA9145 (n = 52), 3.8 for n = 12, 2.3 for n = 3. Counts of rare-present modules
travel better than bit totals.

## Contamination vs real accessory content

| Check | TMED14 | UBA9145 | Other 7 species with rare-present modules |
|---|---|---|---|
| Rare-present modules (unique) | 68 | 39 | 50 |
| Anchored on a contig with a persistent gene | 41 | 12 | **47** |
| Housekeeping genes among them | ribosomal L4/L6/L18/L24, gyrase A, ATP synthase beta, IF-2, NusA, EF-P, RadA, DnaB | release factor 1, DNA pol III epsilon, topoisomerase, Glu-tRNA ligase, Asp-tRNA amidotransferase A | GroEL (Thioglobus_A, next to RuBisCO) |
| Coverage r below the core contigs' 5th percentile | 29 / 68 | 20 / 39 | **0 / 50** |
| Module/core depth ratio < 0.5 | 15 | 14 | 4 (all Pontimonas) |

**TMED14.** The bins are chimeric. GTDB-Tk aligns only 68-71% of each bin to the representative (ANI
99.5%), CheckM2 reports 17.6% contamination, and the rare modules are dense in core genes of a second
Flavobacteriales genome. Some of the contaminant's genes cluster with TMED14 persistent families at
PPanGGOLiN's identity threshold, which is why even TMED14 "core" contigs show a wide coverage spread
(r 5th percentile 0.54).

**UBA9145.** The bins are probably also chimeric, but less clearly so:

- 27 of the 39 modules are unanchored, and 20 have decoupled coverage;
- several carry housekeeping genes that none of the 52 public genomes have in that context;
- GTDB AF is 0.76.

Pfam now resolves the flagellar sets into three modules:

- module_6: FliA sigma factor, two FHIPEP (FlhA-family) proteins, a MinD/ParA-type ATPase (FlhG-like);
- module_20: FlgM, FlgA, FlgC, FlaE, the hook cap;
- module_29: flagellin, FlaG.

All three are unanchored. Their coverage correlates with the core (r 0.86-0.96, near the core contigs'
own range) but at 0.38-0.62× the core depth, which fits a co-occurring motile relative at about half the
abundance. A motile half of the population is still
possible, but the unanchored housekeeping genes make contamination the better explanation. The earlier
"motile local population" reading in `first_tests.md` should be treated as unlikely.

**Pseudothioglobus sp913057955** (3-6-genome tier). One rare-present module is likely contamination:
module_12, 18 families, a full NADH:quinone oxidoreductase (*nuo*) set plus Lon and biotin ligase. It is
unanchored, on one contig, at 6.7× the core depth. A coverage ratio far above 1 flags a contig from a more
abundant organism.

The anchoring and coverage checks agree, and both are cheap. GTDB-Tk AF < 0.8 together with housekeeping
genes in rare modules flagged both bad bins before any coverage was examined.

## What the genuine rare modules are

The table below counts the 553 families of rare-present modules from both tiers, TMED14 excluded. Each
family is classed by its Bakta product, or by its Pfam domains where Bakta gives nothing.

| Class | Families (Bakta only) | Bakta, else Pfam |
|---|---|---|
| no annotation (hypothetical / no hit; with Pfam: no informative domain) | 436 | 319 |
| other (product not matching any class) | 17 | 0 |
| metabolism / enzyme | 55 | 126 |
| transport | 16 | 29 |
| CAZymes / glycan use | 5 | 17 |
| regulation / two-component | 6 | 15 |
| motility / pili (mostly UBA9145) | 4 | 15 |
| housekeeping (mostly UBA9145) | 6 | 12 |
| phage / IS / integrase | 6 | 9 |
| secretion (T9SS C-terminal sorting, HlyD) | 1 | 7 |
| cell surface / polysaccharide | 1 | 3 |
| defence / TA | 0 | 1 (AbiEii) |
| CRISPR-Cas | 0 | 0 |

Pfam rescued 134 of the 453 families Bakta left hypothetical. Even so, **58% of rare-module families
have no annotation**, against 75% of all CDS in the assembly under this Bakta run. Rare modules are, if
anything, slightly better annotated than average.

Recurring themes among anchored, coverage-consistent modules (`pfam/summary.md` lists all of them):

- **Algal and storage carbohydrates.** The species and modules involved:
  - UBA10364 module_0: 13 families, alginate lyase, chondroitinase B / heparinase-like, KdgK, hexuronate
    transporter, SusD (k = 1/11);
  - UBA10364 module_4: an alginate lyase with a T9SS sorting domain;
  - UBA9145 module_7 and module_18: alginate / poly-mannuronate lyase, mannitol dehydrogenase (k = 3/52);
  - Pseudothioglobus sp913057955: a mannitol ABC transporter with mannitol dehydrogenase;
  - Aquiluna: alpha-amylase plus regulator, and beta-glucosidase plus solute-binding protein.
  
  This is the strongest pattern: several unrelated species carry polysaccharide (mainly alginate)
  degradation or sugar-alcohol uptake that most public genomes of the species lack. That fits a
  macroalgae- or ice-algae-influenced coastal site.
- **ABC sugar and other transporters.** The main examples are in Pseudothioglobus sp905182815:
  - module_11: maltose/maltodextrin CUT1 transporter with 2-deoxy-scyllo-inosamine dehydrogenase;
  - module_4: sn-glycerol-3-phosphate UgpABCE;
  - module_10: 13 families including arylsulfatase, a biphenyl-diol dioxygenase and a family-7
    solute-binding protein.
  
  Further examples are a Pontimonas MFS/ABC module, an Aquiluna permease, and a Pseudothioglobus
  sp913057955 xanthine permease.
- **Proteases and secreted proteins in Bacteroidota.** UBA9320 has subtilase and caspase-domain modules
  and a T9SS-sorted protein. MAG-121220 has M28/M42 peptidases and a PA-domain T9SS protein. Pontimonas
  has a trypsin-like protease.
- **Regulators and two-component systems.** Amylibacter, Pontimonas (2 modules, one on a geNomad
  plasmid contig) and Marimicrobium (AraC).
- **IS elements and integrases.** Marimicrobium has an IS3 transposase plus integrase (k = 0/10), an
  integrase plus MOSC module, and a transposase plus thrombospondin-repeat module. These are the only
  clearly mobile modules outside the contaminated bins.
- **Single notable cases.**
  - Thioglobus_A module_13 has RuBisCO large chain plus GroEL in a module absent from all 6 public
    genomes. Four of the six are single-cell genomes (median completeness 0.76), so this is probably missed or fragmented
    rather than truly absent.
  - Marimicrobium has a type IV minor-pilin (ComP-like DNA-uptake) module.
  - MAG-121220 has a FabF/ACP module.

Mobile elements are rare among surprising modules. Only 3 of the 157 rare-present modules in the 13-species
tier sit on geNomad contigs:

- Pontimonas module_6 (plasmid): LuxR regulator, trypsin, response regulator;
- UBA9145 module_50 (virus): hypothetical proteins;
- TMED14 module_63 (virus): hypothetical proteins.

No CRISPR-Cas gene appears in any rare-present module. Phage structural genes (capsid, terminase, tail)
are absent too.

## Coverage of genuine modules

In the 8 main-assembly bins outside TMED14 and UBA9145, no module's coverage correlation falls below the
5th percentile of the bin's own core contigs. Module/core depth ratios have a 10th-90th percentile range of 0.62-1.17 (48 modules), so these modules
are carried by most cells of the local population. They are not minority strain content.

The exceptions are 4 Pontimonas modules at 0.03-0.32× the core: module_12 (0.03), module_7 (0.27; an
aldehyde dehydrogenase / deoxyribose-phosphate aldolase module), module_15 and module_16. These are the best
candidates for content that is variable within the population. They need checking with read-level
evidence before calling them that.

The coverage numbers are weak overall. Nanopore depth is low, with the core summing to 3-20× across all
36 samples, and most bins reach >= 1× in only one or two samples. That is why the ratio is computed over
summed depth rather than per sample. Depth was also a binning input, so contigs in a bin are pre-selected
to co-vary; the 5th-percentile calibration against core contigs only partly controls for this.

## Caveats

- **Public genomes are MAGs and SAGs.** None of the 19 species has an isolate genome in GTDB r226. All
  are metagenome-derived or single-cell, except one MAG-120531 genome with no NCBI category. Pseudothioglobus sp905182815 has 20 of 23
  genomes as Boothbay SAGs (median completeness 0.81), and BACL14 and CACCCN01 are half SAGs.
  - Prevalence is not corrected for public-genome completeness, so modules can be "absent" because a
    public genome is incomplete or fragmented. Shell and cloud modules then look rarer than they are, and
    presence is over-scored.
  - PPanGGOLiN modules also need co-located families, which favours long contigs. Our nanopore contigs
    are longer than most public short-read MAG contigs.
- **Environment.** Most public genomes come from temperate coastal waters: Baltic Sea, North Sea,
  Chesapeake Bay, Saanich Inlet, Gulf of Maine (Boothbay). Chesterfield bins being "surprising" partly
  measures geography and habitat, which is the point, but it also measures sampling bias, so it should
  not be read as novelty of the population. The GTDB records carry too little metadata (about half have
  no isolation source) to stratify by environment.
- **Overlapping assemblies.** Main and marine bins of a species are not independent, so k = 0 in both
  counts once.
- **Few public genomes.** With n = 8-12 (most species), and especially n = 3-6, k/n <= 0.1 means k = 0,
  and one public genome changes the answer.
- **Module definitions.** PPanGGOLiN modules are per-pangenome and are not comparable across species.
  Neighbouring modules on one island are counted separately; the 50 genuine modules sit on fewer distinct
  contigs, as several Pontimonas and MAG-121220 modules share one.
- **Function.** Function comes from Bakta products (75% hypothetical in this run) plus Pfam. Keyword
  classes are coarse: one gliding-motility gene can label a 48-family TMED14 module "motility". The
  per-family counts above are more reliable than the per-module labels in `flagged_modules.tsv`.
- **Marine-only bins.** Bins from the marine assembly have no Bakta or geNomad of their own. They get
  function through representative-protein hits and MGE status only through >= 95% identity hits into the
  main assembly.

## For grex

Nothing needs the cluster. The whole survey, 19 pangenomes plus searches, took about 1 CPU-hour at
2 threads. If the coordinator wants more:

- **Bakta on marine bins:** run Bakta (full DB) on the marine and freshwater bins, so their genes get
  products and coordinates directly.
- **Genus-level pangenomes for freshwater:** freshwater species are singletons in GTDB, so this is the
  only way to add them. That needs a different prevalence model (shrinkage up the tree, as in the README).

## Next steps

- Make anchoring, GTDB-Tk AF and coverage decoupling part of the surprise model. A rare-present module
  should be discounted when its contig has no persistent genes and its coverage departs from the core.
- Correct prevalence for public-genome completeness and contig fragmentation, as is already done for our
  bins.
- Recheck the Pontimonas low-ratio modules and the UBA9145 flagellar sets with read-level evidence:
  reads that span from module contigs into core contigs.

## Files

All in `/data/scratch/metapanclust_data/pangenome_survey/`:

| File | Contents |
|---|---|
| `<species>/ppanggolin/` | one pangenome per species |
| `all/`, `small/` | per tier: `per_bin.tsv`, `flagged_modules.tsv`, `flagged_genes.tsv`, `module_coverage.tsv`, `summary.md`, `coverage.md` |
| `pfam/` | `module_functions.tsv` and `summary.md` |
| `public_meta.tsv` | isolation source and country of the public genomes |
