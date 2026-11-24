# Short-read (Illumina) strain resolution from SNVs and coverage, 2015-2026

Verification legend: **[V]** = author/year/venue confirmed this session from the publisher, PMC/Europe PMC, Crossref or an indexed record. **[P]** = partly confirmed (e.g. venue/year confirmed, details from snippet only). **[U]** = from prior knowledge, NOT verified this session; check before citing.
GitHub "last push" dates were pulled from the GitHub API on 2026-11-24 (`gh api repos/<repo>`), and are given as the maintenance signal.

## Q1. Per-method inventory: input, de novo vs reference, samples, depth, max strains, haplotype scope, accessory genes, compute, maintenance

### Takeaway
Only two published tools both start from your own MAG/co-assembly and use multi-sample covariation to phase strains de novo: DESMAN (SNV frequencies on single-copy core genes, then gene presence/absence) and STRONG/BayesPaths (unitig paths on the assembly graph). Both need many samples (DESMAN: about 30 or more) and give haplotypes for core genes, not whole genomes. Database-matching profilers (StrainGE, StrainScan, StrainEst, StrainPhlAn) work at very low depth but only find strains close to references, which will be rare for novel environmental taxa. inStrain is a within-population diversity and strain-comparison tool, not a haplotype reconstructor.

### Cited Findings

#### DESMAN [V]
- Citation: Quince C, Delmont TO, Raguideau S, Alneberg J, Darling AE, Collins G, Eren AM. "DESMAN: a new tool for de novo extraction of strains from metagenomes." *Genome Biology* 2017; 18:181; doi:10.1186/s13059-017-1309-9 — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5607848/fullTextXML); [DOAJ record](https://doaj.org/article/07347757eb3747c6876a2669a41ee84b)
- Input: a MAG/bin (from co-assembly plus binning, e.g. CONCOCT) and reads mapped from multiple samples. The method finds single-copy core genes (SCGs) present in every strain, counts base frequencies at each position across samples, and uses variant co-occurrence across samples to infer strain haplotypes and abundances. It is de novo: no reference database is needed. — [PMC5607848 summary via search](https://pmc.ncbi.nlm.nih.gov/articles/PMC5607848); [Univ. Galway record](https://researchrepository.universityofgalway.ie/entities/publication/4ed18769-435c-4a34-86c6-c45162b998df)
- Accessory genes: the inferred strain profiles are then used to call non-core gene presence/absence per strain. Average gene presence/absence accuracy was 94.9% across 67 complex-mock haplotypes, rising to 97.39% for haplotypes with SNV error below 1% (the Fig. 2 caption gives 95.7% for the E. coli mock, so text and figure disagree slightly). — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5607848/fullTextXML)
- Benchmarks:
  - E. coli mock (64 samples; 5 E. coli strains plus 15 other gut genomes): all 5 strains recovered. SNV error was 0.17–0.64% per strain, with 99.58% of variable positions correct. Relative frequencies were reproduced with R² = 0.9998. Out of 6,044 variant positions, 1,000 were randomly subsampled for Gibbs sampling.
  - Complex mock (210 genomes / 100 species / 96 samples): used only the 36 universal SCGs.
  - E. coli O104:H4 outbreak: 20 of 53 faecal samples used (mean core-gene coverage above 5), with 440 single-copy core species genes. The outbreak strain was recovered at 99.8% core sequence accuracy (preprint figure). — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5607848/fullTextXML); [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/073825v1.full.pdf)
- Samples needed: error rises below about 30 samples and reaches roughly 15% with 10 samples. — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5607848/fullTextXML)
- Depth: there is no hard minimum. Samples with mean core-gene coverage of 5 or less were excluded in the outbreak data. In Tara Oceans (32 MAGs, 93 samples), the chance of detecting a strain rose with that strain's coverage (logistic regression p = 0.0035). Some strains were resolved from as few as 10 SNVs. — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5607848/fullTextXML)
- Number of strains G: chosen heuristically, from the elbow of posterior mean deviance over a range of G (3–8 tested on the mock, 2–10 on the outbreak). An automated heuristic combining deviance and mean SNV uncertainty was used for complex data. — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5607848/fullTextXML)
- Compute: one G=5 run on the 20-genome mock took about 117 min on 1 core. Runtime scales roughly linearly with the number of samples, and replicate runs and different G values can run in parallel. — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5607848/fullTextXML)
- Maintenance: github.com/chrisquince/DESMAN was last pushed 2020-06-17 (72 stars). Effectively unmaintained. — GitHub API, 2026-11-24

#### STRONG / BayesPaths [V]
- Citation: Quince C, Nurk S, Raguideau S, James R, Soyer OS, Summers JK, Limasset A, Eren AM, Chikhi R, Darling AE. "STRONG: metagenomics strain resolution on assembly graphs." *Genome Biology* 2021; 22:214; doi:10.1186/s13059-021-02419-7 — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8311964/fullTextXML); [DOAJ](https://doaj.org/article/c74ea889eddf4483a41b363af8e7877a)
- Input: raw reads from multiple samples. The pipeline co-assembles with metaSPAdes and bins with CONCOCT (default) or MetaBAT2. It keeps the high-resolution assembly graph from before variant simplification, extracts the subgraph for each single-copy core gene in each MAG, and threads reads onto it to get per-sample unitig coverages. BayesPaths, a variational Bayesian algorithm, then infers the number of strains, each strain's coverage in each sample, and each strain's haplotype path through every SCG subgraph. — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8311964/fullTextXML)
- Important for a Flye/metaMDBG workflow: STRONG needs its own metaSPAdes de Bruijn graph that has not been simplified. It cannot take a MAG from an external assembler directly. (Inference; follows from the pipeline description above.)
- Synthetic benchmark (vs DESMAN and mixtureS):

  | Metric | STRONG | DESMAN | mixtureS |
  |---|---|---|---|
  | Strains found | 216 | 208 | 154 |
  | Strains repeated | 6 | 30 | 64 |
  | Per-base error | 0.052% | 0.176% | 48.7% |
  | Abundance adjusted R² | 0.85 | 0.81 | — |

  STRONG was the best at predicting strain number but slightly underestimated it. — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8311964/fullTextXML)
- Depth: no hard minimum is given. Strains above about 20x coverage were resolved very reliably in the low-diversity synthetic set. Missed strains fell below the minimum coverage of any detected strain (5.68x). — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8311964/fullTextXML)
- Max strains: performance drops at high strain diversity, and the authors give about 8–9 strains per species as the probable upper limit. — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8311964/fullTextXML)
- Scope: resolves only SCGs (36 bacterial/archaeal COGs by default). It does not assign accessory genes to strains, is not extended to whole genomes, and does not support eukaryotes. — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8311964/fullTextXML)
- Compute: BayesPaths took up to 36 h on the synthetic data. The anaerobic digester (AD) dataset took over 16 days end to end on 64 threads, about 60% of it in BayesPaths. Only the metaSPAdes co-assembly needed substantial memory. — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8311964/fullTextXML)
- Nanopore validation (AD time series, Bin_72, COG0532): two strains found in Nanopore reads matched STRONG haplotypes at 100% identity over 2,313 bp. The predicted frequency was 28.0% against 27.6% observed. Agreement was weaker for lower-coverage MAGs. — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8311964/fullTextXML)
- Maintenance: github.com/chrisquince/STRONG was last pushed 2023-07-18 (50 stars). — GitHub API, 2026-11-24

#### ConStrains [P]
- Citation: Luo C, et al. "ConStrains identifies microbial strains in metagenomic datasets." *Nature Biotechnology* 2015 (33:1045–1052 [U]); doi:10.1038/nbt.3319 — [Springer link](https://link.springer.com/10.1038/nbt.3319); [PMC4676274](https://pmc.ncbi.nlm.nih.gov/articles/PMC4676274)
- Method: reference-guided. It uses SNP patterns in universal bacterial/archaeal marker genes (MetaPhlAn markers) to infer conspecific strains and their phylogeny. Each strain is encoded as a concatenated SNP string ("uniGcode"), so it works when a species has only one reference genome. Output is marker-gene haplotypes, not whole genomes. — [search summary of nbt.3319 / PMC4676274](https://pmc.ncbi.nlm.nih.gov/articles/PMC4676274)
- [U] Depth: commonly cited as needing roughly 10x or more species coverage. Not verified this session.
- Maintenance: no GitHub repo found under metagenomics-strains/constrains (404). The code was distributed via Bitbucket [U]. Treat as unmaintained.

#### StrainFacts [V]
- Citation: Smith BJ, Li X, Shi ZJ, Abate A, Pollard KS. "Scalable Microbial Strain Inference in Metagenomic Data Using StrainFacts." *Frontiers in Bioinformatics* 2022; 2:867386; doi:10.3389/fbinf.2022.867386 — [Frontiers full text](https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2022.867386/full)
- Input: per-sample allele counts at biallelic SNP sites for one species, from GT-Pro (reference SNP sets derived from UHGG). Inference of strain genotypes is de novo, but the SNP sites are fixed by the reference database. Samples with under 5% horizontal coverage are dropped, and species with over 5,000 SNPs are downsampled to 5,000. — [Frontiers](https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2022.867386/full)
- Model: genotypes are "fuzzy" (continuous rather than discrete), making the model differentiable so it can be fit by gradient descent. The authors draw the analogy to non-negative matrix factorization (NMF). — [Frontiers](https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2022.867386/full)
- Scale: simulations used 200–10,000 samples and 40–500 strains. With 400 samples, 120 strains and 250 SNPs, StrainFacts took 17.2 min on CPU or 5.1 min on GPU, against 6.4 h for Strain Finder. Fitting 10,000 samples × 400 strains used 7.7 GB of memory. It was applied to more than 10,000 human stool metagenomes. — [Frontiers](https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2022.867386/full)
- Accuracy: similar to Strain Finder on UniFrac and Hamming distance, slightly worse on Bray-Curtis, and substantially better on Shannon-entropy error. MixtureS was worse than both. — [Frontiers](https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2022.867386/full)
- Output and limits: genotypes at SNP sites only, essentially core-genome SNPs (non-core SNPs are listed as future work). Only biallelic sites are modelled. Strain number and hyperparameters were set heuristically. Under-specifying strain number (0.8× the true count) "degraded performance dramatically". A strain at 3% abundance had poor concordance with single-cell genotypes. — [Frontiers](https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2022.867386/full)
- Maintenance: github.com/bsmith89/StrainFacts was last pushed 2024-09-29 (12 stars). — GitHub API, 2026-11-24
- Related, same group: StrainPGC (github.com/bsmith89/StrainPGC, last pushed 2024-11-09) assigns gene content to strains using MIDAS v3. Citation not verified [U]: Smith BJ et al., "Accurate estimation of intraspecific microbial gene content variation in metagenomic data with MIDAS v3 and StrainPGC", reportedly *Genome Research* 2025.

#### PStrain [V]
- Citation: Wang S (first author per Crossref), et al. "PStrain: an iterative microbial strains profiling algorithm for shotgun metagenomic sequencing data." *Bioinformatics* 2020; 36(22-23):5499–5506; doi:10.1093/bioinformatics/btaa1056 — [Crossref](https://api.crossref.org/works/10.1093%2FBIOINFORMATICS%2FBTAA1056); [CityU record](https://scholars.cityu.edu.hk/en/publications/pstrain-an-iterative-microbial-strains-profiling-algorithm-for-sh/)
- Method: reference/marker-based, working on SNVs in MetaPhlAn2 marker genes. It iteratively optimizes strain profiles using genotype frequencies plus reads that span multiple SNV loci (read-level phasing). It runs on single samples. The authors report improvements of 87.75% for strain abundance and 59.45% for genotypes over then state-of-the-art tools. — [CityU record](https://scholars.cityu.edu.hk/en/publications/pstrain-an-iterative-microbial-strains-profiling-algorithm-for-sh/)
- Code: github.com/wshuai294/PStrain per the paper. Push date not checked (the wanghaoyu1994 path returned 404). — [CityU record](https://scholars.cityu.edu.hk/en/publications/pstrain-an-iterative-microbial-strains-profiling-algorithm-for-sh/)

#### inStrain [V]
- Citation: Olm MR, Crits-Christoph A, Bouma-Gregson K, Firek BA, Morowitz MJ, Banfield JF. "inStrain profiles population microdiversity from metagenomic data and sensitively detects shared microbial strains." *Nature Biotechnology* 2021; 39:727–736; doi:10.1038/s41587-020-00797-0 — [CU Boulder record](https://experts.colorado.edu/display/pubid_512638); [PMC9223867](https://pmc.ncbi.nlm.nih.gov/articles/PMC9223867)
- Input: a BAM (bowtie2) of reads mapped to a set of representative genomes, which can be your own dereplicated MAGs (`--database_mode` for many genomes). Per sample it computes nucleotide diversity, SNSs/SNVs, linkage, pN/pS and iRep. The recommended rule for the database is genomes no more than 98% and at least 95% ANI apart, e.g. dereplicated with dRep. — [inStrain docs, Important concepts](https://instrain.readthedocs.io/en/latest/important_concepts.html)
- Comparisons between samples:
  - conANI compares consensus alleles only.
  - popANI counts a site as a difference only if the two samples share no allele, major or minor.
  - Compared bases must have at least 5x coverage in both samples.
  - The suggested same-strain call is popANI ≥ 99.999% with ≥ 50% of the genome compared. Presence is called at ≥ 50% breadth, and breadth ≈ 1 − e^(−0.883 × coverage).
  - [inStrain docs](https://instrain.readthedocs.io/en/latest/important_concepts.html); [PMC9223867](https://pmc.ncbi.nlm.nih.gov/articles/PMC9223867)
- Scope: it describes and compares populations (sample consensus versus population) and does not phase multiple strains within a sample into haplotypes. Linkage is reported as a metric. — [inStrain docs](https://instrain.readthedocs.io/en/latest/important_concepts.html)
- Maintenance: github.com/MrOlm/inStrain was last pushed 2026-10-20 (194 stars). Active. — GitHub API, 2026-11-24

#### StrainPhlAn 4 / MetaPhlAn 4 [P]
- The original StrainPhlAn is Truong DT, Tett A, Pasolli E, Huttenhower C, Segata N, *Genome Research* 2017 (authors from memory [U], venue/year [P]). It reconstructs a consensus (dominant-strain) sequence in species-specific marker genes and builds strain phylogenies. It was applied to over 1,500 gut metagenomes covering more than 125 species. — [PMC5378180](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5378180/)
- [U] StrainPhlAn 4 / MetaPhlAn 4: Blanco-Míguez A, et al., "Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4", *Nature Biotechnology* 2023 (41:1633–1644). Not retrieved this session. Its marker database includes SGBs (species-level genome bins) built from MAGs, so some uncultured taxa have markers. Output is still one dominant-strain consensus per sample, not a mixture.
- Maintenance: github.com/biobakery/MetaPhlAn was last pushed 2026-11-16 (436 stars). Active. — GitHub API, 2026-11-24

#### StrainGE [V]
- Citation: van Dijk LR, Walker BJ, Straub TJ, Worby CJ, Grote A, Schreiber HL, ... Earl AM. "StrainGE: a toolkit to track and characterize low-abundance strains in complex microbial communities." *Genome Biology* 2022; 23:74 (volume and article number [P]); doi:10.1186/s13059-022-02630-0 — [PubMed 35255937](https://pubmed.ncbi.nlm.nih.gov/35255937/); [DOAJ](https://doaj.org/article/ebea40720d8445e09bfcaeec17a670df)
- Method: reference-database based, with two parts. StrainGST identifies the closest reference strain(s) by k-mers, and StrainGR calls variants against them. It can identify strains at 0.1x coverage and call variants for multiple strains of one species from as little as 0.5x. Demonstrated on E. coli and Enterococcus. — [PubMed 35255937](https://pubmed.ncbi.nlm.nih.gov/35255937/)
- Requires a curated reference panel per species, clustered so that only one representative is kept per cluster (the StrainScan paper criticizes this). — [StrainScan PMC10433603](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10433603/)
- Maintenance: github.com/broadinstitute/StrainGE was last pushed 2024-10-21. — GitHub API, 2026-11-24

#### StrainScan [V]
- Citation: Liao H, Ji Y, Sun Y. "High-resolution strain-level microbiome composition analysis from short reads." *Microbiome* 2023 (11:183 [U on the number]); doi:10.1186/s40168-023-01615-w — [PMC10433603](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10433603/); [Microbiome PDF](https://microbiomejournal.biomedcentral.com/counter/pdf/10.1186/s40168-023-01615-w)
- Method: a hierarchical tree-based k-mer index over known strain genomes. It resolves strains within clusters that StrainGE and StrainEst collapse to one representative, and reports about 20% higher F1 for identifying multiple strains. It works on single samples and is reference-only, so it cannot find novel strains. — [PMC10433603](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10433603/)
- Maintenance: github.com/liaoherui/StrainScan was last pushed 2024-04-08. — GitHub API, 2026-11-24

#### StrainEst [P]
- Citation: Albanese D, Donati C. "Strain profiling and epidemiology of bacterial species from metagenomic sequencing." *Nature Communications* 2017 (8:2260 [U on the number]); doi:10.1038/s41467-017-02209-5 — [PMC5741664](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5741664/); [Springer](https://link.springer.com/article/10.1038/s41467-017-02209-5)
- Method: reference-based. It clusters SNV profiles of a species' reference genomes, then estimates the number of strains present and their abundances by fitting sample SNV frequencies with Lasso regression [Lasso is U]. It works on single samples. — [PMC5741664](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5741664/)
- Maintenance: github.com/compmetagen/strainest was last pushed 2020-01-17. Unmaintained. — GitHub API, 2026-11-24

#### Newer methods, 2022–2026
- **MetaStrainer [V]:** Sharaf H, Bobay L-M. *Bioinformatics* 2026; 42(6):btag340; doi:10.1093/bioinformatics/btag340 (published 24 May 2026).
  - Input and method: one reference genome (GenBank; can be a contig-level assembly, so a MAG would fit) plus one paired-end sample. Reads are mapped with Bowtie2, variants are linked into linkage groups, and an MCMC (emcee) estimates strain composition and phasing. Strains at ≥ 99.5% ANI are merged. It runs on a single sample.
  - Accuracy on simulated *Gilliamella apicola* data (CAMISIM): 92.1% of variants called correctly, against 39.3% for mixtureS. Strain number was correct in 53/56 simulations, against 4/54. Abundance accuracy was 99.1% against 71.3%. Accuracy held down to about 8x (about 5% frequency) for some strains. Each allele needs at least 20 supporting reads.
  - Limits: at most 3 strains, and accuracy drops with 4 or more or with similar frequencies (64.2% in the equal-frequency case). Only core and reference genes are covered. Cross-mapping from related species is a risk. Validated on one simulated species only.
  - Code: github.com/lbobay/MetaStrainer (pushed 2026-11-24).
  - — [Bioinformatics](https://academic.oup.com/bioinformatics/article/42/6/btag340/8691842)
- **Floria [P]:** Shaw J, Yu YW. "Floria: fast and accurate strain haplotyping in metagenomes." *Bioinformatics* 2024 (ISMB issue; PMC11211831). It clusters reads into strain-level haplotype sets, works on short or long reads, and is reference/MAG-based. On synthetic metagenomes it was more than 3× faster than Strainberry and recovered 21% more strain content. Phasing takes minutes on a workstation. Repo: bluenote-1577/floria, last pushed 2025-04-10. — [PMC11211831](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11211831/); [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.01.28.577669v1.full)
- **StrainXpress [P]:** a strain-aware de novo assembler for short reads, *Nucleic Acids Research* 2022; 50(17):e101. Its benchmark included ConStrains, StrainFinder, Gretel and DESMAN. Author list not verified (Kang X et al. [U]). — [NAR](https://academic.oup.com/nar/article/50/17/e101/6625806)
- **E. coli short-read profiler benchmark (bioRxiv 2026; v1 posted 2026-05-19) [P]:** compared PanTax, PathoScope 2.0, StrainGE, Strainify, StrainR2 and StrainScan. StrainR2 had the top F1 (0.995), and the tools gave noticeably different strain compositions. Full text could not be fetched (HTTP 429), so authors and details are unverified. — [bioRxiv v2](https://www.biorxiv.org/content/10.64898/2026.05.19.726160v2)
- **mixtureS [U]:** a single-sample, reference-based mixture model used as the weak baseline in STRONG, StrainFacts and MetaStrainer. Citation not verified this session (Li X, Saadat S, Hu H, Li X, *Bioinformatics* 2020 [U]). — performance numbers from [STRONG](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8311964/fullTextXML) and [MetaStrainer](https://academic.oup.com/bioinformatics/article/42/6/btag340/8691842)

### Inferences
- How each method could fit a "bin first, then resolve strains" workflow built on a Flye/metaMDBG + Illumina co-assembly:
  - **DESMAN** fits directly. Map Illumina reads from every sample to the MAG's contigs, then run it on SCG variant frequencies.
  - **MetaStrainer and Floria** fit directly, but run per sample; Floria can also use the long reads.
  - **inStrain** fits directly for microdiversity and for same-versus-different strain calls between samples, not for haplotypes.
  - **STRONG** does not fit without re-co-assembling the short reads in metaSPAdes.
  - **StrainFacts** needs a GT-Pro SNP database. Building one for novel MAG species is possible in principle but not documented here.
- For taxa whose strains differ mainly by SNVs and that are sampled across many samples, DESMAN is the only maintained-enough de novo option that also assigns accessory genes, though its repo has been idle since 2020. STRONG is more accurate (per-base error 0.052% against 0.176%) but is limited to SCGs and costs days of compute.
- Reference-DB profilers (StrainGE, StrainScan, StrainEst, PStrain, ConStrains, StrainPhlAn) will mostly fail on novel environmental lineages that lack close references.

### Gaps
- Not retrieved this session: exact ConStrains depth thresholds; StrainPhlAn 4 paper details; the StrainScan article number; StrainEst numerical accuracy; author list of the 2026 E. coli benchmark.
- No CAMI-based head-to-head benchmark of the de novo haplotypers (DESMAN, STRONG, StrainFacts) was found. CAMI challenges assess strain-level binning and profiling, not haplotype phasing.
- Reviews named in the brief were not fetched, so these citations are unverified [U]:
  - Anyansi C, Straub TJ, Manson AL, Earl AM, Abeel T, "Computational methods for strain-level microbial detection in colony and metagenome sequencing data", *Front Microbiol* 2020; 11:1689.
  - Van Rossum T, Ferretti P, Maistrenko OM, Bork P, "Diversity within species: interpreting strains in microbiomes", *Nat Rev Microbiol* 2020; 18:491–506.

## Q2. How the methods use multi-sample coverage covariation; the shared-core assumption

### Takeaway
The de novo methods (DESMAN, STRONG/BayesPaths, StrainFacts, and Strain Finder before them) all treat each sample as a different mixture of the same G strains. Strain haplotypes are the latent factors that explain how allele (or unitig) frequencies covary across samples, so they need strain proportions to vary between samples. All of them work on a shared core (SCGs or core SNP sites) and handle accessory genes, if at all, in a second step.

### Cited Findings
- DESMAN works on genes "present as a single copy in every strain" and links variants into haplotypes by their co-occurrence across samples. Strain profiles then assign non-core genes. — [PMC5607848](https://pmc.ncbi.nlm.nih.gov/articles/PMC5607848); [Univ. Galway record](https://researchrepository.universityofgalway.ie/entities/publication/4ed18769-435c-4a34-86c6-c45162b998df)
- DESMAN assumes variant positions are independent and does not use variant co-occurrence within reads (no read-level phasing). — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5607848/fullTextXML)
- In STRONG, strain separation depends on differences in strain proportions between samples. If proportions are stable, haplotype accuracy falls, although abundances are still predicted well. — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8311964/fullTextXML)
- StrainFacts factorizes a samples × SNP allele-count matrix into strain fractions times fuzzy genotypes, an NMF analogue, using biallelic core SNPs. — [Frontiers](https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2022.867386/full)
- Single-sample methods use read-level linkage instead of cross-sample covariation:
  - PStrain uses reads that cover several SNV loci — [CityU](https://scholars.cityu.edu.hk/en/publications/pstrain-an-iterative-microbial-strains-profiling-algorithm-for-sh/)
  - MetaStrainer uses linkage groups plus MCMC — [Bioinformatics](https://academic.oup.com/bioinformatics/article/42/6/btag340/8691842)
  - Floria clusters reads — [PMC11211831](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11211831/)

### Inferences
- Co-assembling several samples from one environment gives the cross-sample variation these models need. It only helps if strain ratios actually differ between samples. Replicates or samples taken close together in time with constant ratios will make DESMAN and STRONG unstable.
- With 150-bp reads, read-level phasing can only link SNVs within about 150–500 bp (insert size). In low-diversity taxa (few SNVs per kb), the single-sample methods will produce fragmented haplotypes.

### Gaps
- None of the sources I retrieved quantify how much between-sample variance in strain proportions is needed for DESMAN or STRONG to work.

## Q3. Reported accuracy and known failure modes

### Takeaway
On synthetic data, the best de novo tools reach per-base haplotype errors of 0.05–0.2% for strains above about 20x coverage. Quality drops with fewer than about 30 samples, coverage below 5–6x, more than 8–9 strains, similar or stable abundances, too few SNVs, recombination, and binning or mapping errors.

### Cited Findings
- **DESMAN:**
  - About 15% error with only 10 samples.
  - Failed mapping and fragmented (multi-contig) assemblies caused false negatives.
  - Contaminated bins caused errors.
  - Very few variants make deconvolution hard.
  - Variation below the SNV error rate is not resolved, and multi-copy genes are not handled.
  - [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5607848/fullTextXML)
- **STRONG:**
  - Strains below about 5.7x were missed.
  - Practical upper limit of about 8–9 strains per species.
  - Haplotypes are "metagenome strains" rather than isolate genomes, and results for rapidly recombining organisms may not correspond to real individuals.
  - [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8311964/fullTextXML)
- **StrainFacts:**
  - Under-specifying strain number sharply degrades results.
  - A strain at 3% abundance was poorly reconstructed.
  - Low coverage produces missing values.
  - [Frontiers](https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2022.867386/full)
- **MetaStrainer:**
  - 4 or more strains, or strains at similar frequencies, reduce accuracy (64.2% in the equal-frequency case).
  - Cross-mapping from related species is a risk.
  - [Bioinformatics](https://academic.oup.com/bioinformatics/article/42/6/btag340/8691842)
- **mixtureS** is a consistently poor baseline: 48.7% per-base error in the STRONG benchmark, 39.3% genotype accuracy in the MetaStrainer benchmark, and its accuracy falls with lower coverage. — [STRONG](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8311964/fullTextXML); [MetaStrainer](https://academic.oup.com/bioinformatics/article/42/6/btag340/8691842)
- Reference profilers disagree substantially even on E. coli, one of the best-referenced species: six tools gave noticeably different compositions, and StrainR2 had F1 = 0.995. — [bioRxiv 2026 benchmark](https://www.biorxiv.org/content/10.64898/2026.05.19.726160v2)
- Short reads struggle to resolve similarity between strains, and most haplotyping methods have been designed and tested only on short reads. — [Floria bioRxiv](https://www.biorxiv.org/content/10.1101/2024.01.28.577669v1.full)

### Inferences
- If two strains have near-identical abundance profiles across all samples, the covariation methods cannot separate them. They will be merged or mis-phased, and only read- or long-read linkage can split them.
- For this project, long reads offer the cross-check STRONG used: validate short-read haplotypes against Nanopore reads for SCGs of high-coverage MAGs.

### Gaps
- No systematic benchmark of recombination effects was found. Performance on CAMI II strain-madness data for DESMAN, STRONG or StrainFacts was not found.

## Q4. Newer methods (2023–2026), comparative reviews, and tools that take a MAG/bin as input

### Takeaway
Since 2022 the new short-read tools have been single-sample, reference/MAG-anchored haplotypers (Floria 2024, MetaStrainer 2026) and reference-DB profilers (StrainScan 2023, with StrainR2, Strainify and PanTax in a 2026 benchmark). No newer multi-sample de novo method replacing DESMAN or STRONG turned up. The tools that accept an existing MAG are DESMAN, inStrain, Floria and MetaStrainer; STRONG accepts only its own MAGs.

### Cited Findings
- MetaStrainer accepts a "complete genome or contig-level assembly" in GenBank format as the reference, so a MAG would work after annotation. The authors suggest estimating strain number with Floria first. — [Bioinformatics 2026](https://academic.oup.com/bioinformatics/article/42/6/btag340/8691842)
- Floria phases reads mapped to references or MAGs, from short or long reads. — [PMC11211831](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11211831/)
- inStrain runs on any genome set, including MAGs, given a BAM, and recommends a database of genomes at 95–98% ANI. — [inStrain docs](https://instrain.readthedocs.io/en/latest/important_concepts.html)
- DESMAN was demonstrated on CONCOCT bins from Tara Oceans (32 MAGs, 93 samples). — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5607848/fullTextXML)
- StrainScan (2023) improves on cluster-level reference profilers by about 20% F1. — [PMC10433603](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10433603/)
- Other tools that appeared in searches but were not examined: MAGinator (de novo MAG profiling with strain-level phylogenies, PMC11231285) and StrainFLAIR (variation graphs, PMC8388557). — [MAGinator](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11231285/); [StrainFLAIR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8388557/)

### Inferences
- Practical shortlist for the project:
  1. inStrain on dereplicated MAGs across all samples, to find which taxa have within-sample microdiversity and which strains are shared between samples. Cheap and actively maintained.
  2. DESMAN on the MAGs flagged in step 1, if enough samples exist (about 30 or more) and strain ratios vary between them. This gives SCG haplotypes and accessory-gene presence/absence.
  3. Floria (short and long reads) or MetaStrainer for per-sample phasing where there are few samples. Note MetaStrainer's cap of 3 strains.
  4. STRONG only if you are willing to run a separate metaSPAdes Illumina co-assembly.

### Gaps
- No review dated 2023–2026 that compares the de novo haplotypers head-to-head was found in this session. The 2020 reviews (Anyansi et al.; Van Rossum et al.) were not retrieved and their citations are unverified.
- Not covered: Strainberry and Strainy (long-read), and newer hybrid tools.
