# PanGBank as the source of species pangenomes (2026-11-18)

## Result

- **PanGBank covers only 4 of the 19 Chesterfield species, and for those it holds the same genomes we
  already used.** PanGBank builds a pangenome only for GTDB species with >= 15 genomes. For our species,
  GTDB_all (r226) uses exactly the public genome sets of our own builds: 67/67, 52/52, 23/23 and 17/17
  genomes shared. A bigger, more diverse reference does not exist for these species, because PanGBank
  cannot hold more than GTDB does. None of the 19 is in GTDB_refseq.
- **The latest release (r232) covers only 2 of the 19.** GTDB_all 2.x filters genomes at completeness
  >= 70%, and >= 85% for the species representative. This removes MAG-120531 sp000173115, whose
  representative GCA_000173115.1 is 64% complete by CheckM2. It also removes Patiriisocius sp905182755,
  which falls below 15 genomes. Both species still exist in GTDB r232.
- **Family-level surprise is unchanged.** On the same genomes, the number of a bin's families seen in
  <= 10% of public genomes agrees within 1-20% between projection onto PanGBank and our builds. Examples:
  UBA9145 568 vs 574; Pseudothioglobus 57 vs 66; Patiriisocius 336 vs 352.
- **Module-level surprise drops sharply, but the content has moved, not gone.** UBA9145 falls from
  245 bits to 42 bits per bin, and from 39 to 3 rare-present modules. The cause is that PanGBank modules
  are defined on public genomes only, so a module that only our bins carry cannot appear in them.
  `ppanggolin projection` reports that content as "specific" genes instead:
  - 64 of the 78 survey rare-present (bin, module) pairs in UBA9145 consist entirely of specific genes,
    and all 156 pairs across releases stay rare;
  - grouped into islands (>= 3 consecutive specific genes), UBA9145 has 57-68 islands per bin, 315-382
    genes. At the survey's cost of log2(n+2) bits per unit, that is 328-391 bits, more than the survey's
    module total.
  
  Spurious surprise is therefore not reduced. It is only hidden if specific genes are ignored.
- **The contaminant still stands out, through the specific islands.**

  | Species | Specific genes, % of bin genes | Islands without a persistent gene on their contig |
  |---|---|---|
  | UBA9145 | 16-18% | 40/57 and 53/68 |
  | MAG-120531 | 2-6% | 4/10 |
  | Pseudothioglobus | 2-3% | 0/4 |
  | Patiriisocius | 12% | 3/12 |

  The housekeeping and flagellar modules from the survey are all specific genes here:
  - module_0, DNA pol III;
  - module_26, release factor;
  - module_44, L3;
  - module_6 and module_20, flagellar sets.
- **A sister species does not substitute at default settings.** Projecting onto a congeneric PanGBank
  pangenome at the default 80% identity loses much of the core. Our bins' completeness drops:

  | Our species | Projected onto | Completeness, own build | at 80% identity | at 50% identity |
  |---|---|---|---|---|
  | Amylibacter sp905182285 | Amylibacter sp000153745 | 0.81 / 0.99 | 0.64 / 0.78 | 0.79 / 0.96 |
  | TMED14 sp001438205 | TMED14 sp003507115 | 0.70 / 0.73 | 0.39 / 0.41 | 0.77 / 0.81 |
  | Planktomarina forsetii | P. temperata | 0.28 | 0.18 | 0.25 |

  At 50% identity the core comes back, but 900-1,200 TMED14 genes remain specific. Accessory prevalence
  from a sister species is not this species' prevalence. That makes a sister pangenome a genus-level prior
  for the tree-shrinkage step, not a level-1 species prior.

**Recommendation.** Use PanGBank as an optional source, not as the foundation of the level-1 species
prior:

- fetch PanGBank pangenomes by API where they exist (fast, free, CC BY 4.0, standardised PPanGGOLiN);
- keep building our own pangenomes for the rest, which for marine and freshwater metagenome species is
  most of them (15 of 19 here, and every freshwater species);
- in both cases, score content absent from the reference (projection's specific genes) explicitly as
  k = 0 units, and do not score only PanGBank modules;
- keep anchoring and coverage as the contamination test.

The real benefits are uniform annotation, precomputed partitions, RGPs and spots, and per-genome CheckM2
values in the metadata. It brings no extra genomes.

## Access (verified today)

- **REST API** at `https://pangbank-api.genoscope.cns.fr` (FastAPI; schema at `/openapi.json`). Useful
  endpoints:
  - `/collections/`, with release notes per release;
  - `/pangenomes/?taxon_name=...&substring_taxon_match=true&only_latest_release=false`. Pangenome names
    use underscores (`s__UBA9145_sp001438145`); taxon names use spaces;
  - `/pangenomes/{id}` for metadata;
  - `/pangenomes/{id}/file` for the `.h5` (GET only; HEAD returns 405);
  - `/genomes/?genome_name=GCA_...` for per-genome GTDB r226 and r232 taxonomy plus CheckM/CheckM2.
- **CLI** `pangbank` (PyPI `PanGBank-cli`, also on Bioconda with Mash). Its commands are
  `search-pangenomes --download`, `get-pangenome`, and `match-pangenome`. `match-pangenome` uses Mash
  against a collection sketch to place a genome. It was not installed here; the API was enough. There is
  no bulk archive or Zenodo deposit; bulk download means looping over the API.
- **Web**: https://pangbank.genoscope.cns.fr. **Licence**: data CC BY 4.0. The GTDB taxonomy inside is
  CC BY-SA 4.0.
- **Collections**:

  | Collection | Release | GTDB | Pangenomes | Genomes | PPanGGOLiN |
  |---|---|---|---|---|---|
  | GTDB_all | 1.0.0 | r226 | 4,351 | 334,252 | 2.2.4 |
  | GTDB_all | 2.0.0 / 2.1.0 (latest) | r232 | 4,679 | 382,760 | 2.3.0 |
  | GTDB_refseq | 1.0.0 | r226 | 1,770 | 183,331 | 2.2.4 |
  | GTDB_refseq | 2.0.0 (latest) | r232 | 2,044 | 219,611 | 2.3.0 |

  The `.h5` files from releases 1.0 and 2.x both open in our PPanGGOLiN 2.3.2 (status: annotated,
  clustered, sequences, partitioned, RGPs, spots and modules all present). Clustering used 80% identity
  and 80% coverage, the same as our builds.

## Coverage of the 19 survey species

`api/species_hits.tsv` and `api/genus_hits.tsv` hold the full lists.

| Species | Our public genomes (r226) | PanGBank GTDB_all r226 | PanGBank r232 (2.1.0) |
|---|---|---|---|
| MAG-120531 sp000173115 | 67 | 67 (id 3969) | absent (representative < 85% complete) |
| UBA9145 sp001438145 | 52 | 52 (5902) | 51 (17293); 50 shared with ours |
| Pseudothioglobus sp905182815 | 23 | 23 (4969) | 23 (16229); 19 shared with ours |
| Patiriisocius sp905182755 | 17 | 17 (4508) | absent (< 15 genomes after filter) |
| the other 15 (3-13 genomes) | 3-13 | absent (< 15) | absent |

For 8 of the 15 missing species, PanGBank does hold a pangenome of a sister species in the same genus:

- Amylibacter (2 species);
- Planktomarina (temperata);
- TMED14 (4);
- UBA10364 (4);
- BACL14 (5-6);
- UBA9145 (6 more);
- Thioglobus_A (sp018644675);
- Pseudothioglobus (singularis).

Six genera have nothing: UBA9320, Pontimonas, Marimicrobium, CACCCN01, Hel1-33-131 and
MAG-121220-bin8. Oceaniferula and Aquiluna have nothing either. TMED14 sp001438205, the survey's chimeric
bin, has no pangenome of its own.

## Surprise: PanGBank projection vs the survey builds

The model is the survey's:

- p = (k+1)/(n+2) over the n PanGBank genomes, counting a module at completion >= 0.5;
- surprise is -log2 p when present and -log2(1 - p c) when absent;
- c is the bin completeness reported by projection, computed from persistent families.

"Islands" are runs of >= 3 consecutive specific genes on one contig. Each costs log2(n+2) bits, like a
k = 0 module.

| PanGBank | Bin | n | c (proj / survey) | Rare-present mod. (PB / survey) | Module bits (PB / survey) | Specific genes | Islands (anchored) | Island bits |
|---|---|---|---|---|---|---|---|---|
| MAG-120531 r226 | semibin_019 | 67 | 0.83 / 0.84 | 1 / 2 | 17.5 / 23.6 | 131 | 10 (6) | 61 |
| MAG-120531 r226 | marine_lorbin_024 | 67 | 0.31 / 0.30 | 0 / 1 | 6.1 / 12.2 | 13 | 1 (1) | 6 |
| UBA9145 r226 | maxbin_040 | 52 | 0.92 / 0.92 | 3 / 39 | 42.2 / 245.5 | 519 | 68 (15) | 391 |
| UBA9145 r226 | marine_maxbin_019 | 52 | 0.94 / 0.95 | 3 / 39 | 42.2 / 245.5 | 449 | 57 (17) | 328 |
| UBA9145 r232 | maxbin_040 | 51 | 0.92 / 0.92 | 3 / 39 | 44.8 / 245.5 | 519 | 68 (15) | 390 |
| Pseudothioglobus r226 | lorbin_018 | 23 | 0.71 / 0.74 | 0 / 3 | 20.5 / 36.9 | 38 | 3 (3) | 14 |
| Pseudothioglobus r232 | lorbin_018 | 23 | 0.72 / 0.74 | 0 / 3 | 25.9 / 36.9 | 38 | 3 (3) | 14 |
| Pseudothioglobus r232 | marine_lorbin_007 | 23 | 0.70 / 0.73 | 0 / 3 | 28.3 / 38.2 | 28 | 1 (1) | 5 |
| Patiriisocius r226 | marine_semibin_023 | 17 | 0.83 / 0.83 | 0 / 0 | 21.4 / 20.5 | 320 | 12 (9) | 51 |

- **Projection reproduces the survey's completeness** to within 0.03. This is a useful check that
  projection works on our nanopore bins.
- **Every gene of the survey's flagged modules was matched to a projected gene by coordinates.** Of the
  171 (release, bin, module) rows, 169 stay rare (median family frequency <= 0.1). In r232, the Pseudothioglobus
  Ugp glycerol-3-phosphate module (survey k = 1/23) reaches 0.13 because of the 4 swapped genomes. That
  is the one-genome sensitivity the survey warned about.
- **The genuine anchored modules survive as anchored islands**, for example the Pseudothioglobus
  arylsulfatase/dioxygenase module_10 (13 genes, all specific).
- **Release mismatch matters little within a species.** Between r226 and r232, UBA9145 and
  Pseudothioglobus differ by 1-4 genomes. Bits move by 2-8 per bin.
- **Release mismatch matters a lot for coverage.** The 2.x quality filter drops whole MAG-only species.
  Taxonomy names must be matched to the release: our GTDB-Tk calls are r226, and 2.x is r232. Two
  workarounds: look up by genome accession through `/genomes/`, or pin release 1.0.0.

## Practicality

- **Download.** Each pangenome is 10-29 MB (17-67 genomes, about 0.4-0.6 MB per genome) and takes about
  2 s to fetch. Extrapolating, all of GTDB_all 2.1.0 is about 150-230 GB; that figure is not verified.
  Per-species fetches on demand are the sensible mode.
- **Projection.** It takes 10-30 s wall time per species for 1-2 bins at 3 threads (17-68 CPU-s), with
  about 1 GB peak RSS. Exporting tables from the `.h5` (`write_pangenome --stats --modules --Rtab`)
  takes seconds. Our own builds took 1-3 min per species, plus genome download.
- **Version.** Files from PanGBank (PPanGGOLiN 2.2.4 and 2.3.0) are readable by 2.3.2.
- **Metadata.** PanGBank's metadata has per-genome CheckM2 completeness. This is the input needed for
  the survey's open item "correct prevalence for public-genome completeness".
- **Licence.** CC BY 4.0 data, with GTDB's CC BY-SA 4.0 taxonomy. Citing PanGBank (Mainguy et al., NAR
  2026) and GTDB is enough.

## Method

Code is in `experiments/pangbank/` and data in `/data/scratch/metapanclust_data/pangbank/`. The whole
run is under 5 CPU-minutes.

1. `query_species.py`: queries the API for the 19 species and their genera, across all collections and
   releases. Output is `api/species_hits.tsv` and `api/genus_hits.tsv`.
2. `run_projection.sh`:
   - downloads the 6 `.h5` files (4 species, r226 plus r232 where present);
   - exports module, stats and Rtab tables;
   - runs `ppanggolin projection --table -c 3` on the same bin FASTAs the survey used, from
     `pangenome_survey/<sp>/genomes` and `pangenomes/<sp>/genomes`.
3. `compare.py`: computes per-bin surprise, specific-gene islands and anchoring, and the fate of each
   survey rare-present module. Survey genes are matched to projected genes by contig and coordinates.
   Outputs are in `compare/`: `per_bin.tsv`, `survey_module_fate.tsv` and `specific_islands.tsv`.
4. `cross_species.sh`: projects TMED14, Amylibacter and Planktomarina bins onto a congeneric PanGBank
   pangenome at 80% and 50% identity. Output is in `xsp/`.

## Caveats

- Only 4 species could be tested, and none had more genomes in PanGBank than in our builds. The question
  "does a bigger reference reduce spurious surprise?" therefore could not be tested on Chesterfield data.
  Testing it needs a species with hundreds of PanGBank genomes and a bin of ours, which this site does not
  have.
- Islands of specific genes are a cruder unit than modules:
  - they ignore co-occurrence across genomes;
  - a run of 3 can be split by one gene that hits a PanGBank family;
  - singleton specific genes (most of Patiriisocius' 320) are not scored at all.
  
  Island bits and survey module bits are comparable only in order of magnitude.
- The two bins of a species come from overlapping co-assemblies, as in the survey, and are not
  independent.
- The "specific" call depends on projection's 80% identity and 80% coverage threshold against all
  pangenome genes. Our builds cluster transitively, yet both give nearly the same k = 0 family counts
  (e.g. 519 vs 528).

## For grex

Nothing needs the cluster. If a genus-level prior from PanGBank sister species is wanted for every
Chesterfield genus, `cross_species.sh` scales locally at about 20 s per bin. A full GTDB_all mirror
(about 200 GB) is not justified by this result.
