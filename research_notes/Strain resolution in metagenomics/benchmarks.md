# Benchmark evidence for strain resolution and binning of closely related genomes (2017-2026)

Scope note: PMC, Nature and Springer pages were partly blocked (reCAPTCHA or paywall redirects), and bioRxiv returned HTTP 429 several times. The CAMI I, CAMI II, Strainberry and StrainScan full texts were read. Strainy, myloasm, STRONG and the 2025 binning-evaluation preprint were seen only through abstracts or search snippets, and the notes below mark them that way.

## Q1. CAMI I and CAMI II: assembly and binning of "common strains" (ANI >= 95%) vs unique strains

### Takeaway
Both CAMI rounds found the same thing. Short-read assemblers recover unique strains almost completely but recover only a small fraction of genomes that have a relative at >= 95% ANI. In CAMI I the best median genome fraction for those genomes was 22.5%, against about 98% for unique strains. In CAMI II strain-madness the drop in genome fraction was about 76 percentage points. Binning also degrades on common strains, and no binner wins on all metrics. CAMI II concludes that most assemblers "did not resolve strain variation" and that long or hybrid reads help.

### Cited Findings
**CAMI I** (Sczyrba A, Hofmann P, Belmann P, et al., *Nat Methods* 14(11):1063-1071, 2017, doi:10.1038/nmeth.4458; verified on PMC5903868) — [PMC5903868](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903868/)
- High-complexity dataset, unique strains (<95% ANI): Minia and MEGAHIT had the highest genome fraction (median 98.2%), then Meraga (96.0%) and VELOUR_k31_C2.0 (62.9%) — [PMC5903868](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903868/)
- Common strains (>= 95% ANI): recovery fell for every assembler. The best was MEGAHIT (ep_mtl200 setting) at a median of 22.5%, then Meraga 12.0%, Minia 11.6% and VELOUR 4.1% — [PMC5903868](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903868/)
- Binning, unique strains: MaxBin 2.0 was best on all three datasets (>90% purity, >=70% completeness). MetaBAT, CONCOCT and MetaWatt-3.5 came close on two datasets — [PMC5903868](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903868/)
- Binning, common strains: completeness fell substantially for all programs and purity fell for most. MaxBin 2.0 was the exception and kept >90% purity on all datasets. The page gives no numeric completeness values for common strains — [PMC5903868](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903868/)
- Conclusions: strain-level diversity is "a substantial challenge to all evaluated programs". Closely related strains are "a notable hurdle" for binning. Deconvolving strains may need "substantially larger numbers of replicate samples than those analyzed here (up to 5)" — [PMC5903868](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903868/)
- Inconsistency in the paper: the Fig. 2b caption defines common strains as ">90% ANI", while the Results text uses 95% — [PMC5903868](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903868/)

**CAMI II** (Meyer F, Fritz A, Deng Z-L, et al., *Nat Methods* 19(4):429-440, 2022, doi:10.1038/s41592-022-01431-4; verified on PMC9007738) — [PMC9007738](https://pmc.ncbi.nlm.nih.gov/articles/PMC9007738/)
- Marine dataset, unique vs common genomes: unique genomes did better by +9.7% strain recall, +19.3% genome fraction, 7-fold NGA50 and +6.5% strain precision. Mismatches were 548 (unique) vs 486 (common) per 100 kb. The misassembled-contig fraction was 2.6% (unique) vs 3.1% (common) — [PMC9007738](https://pmc.ncbi.nlm.nih.gov/articles/PMC9007738/)
- Strain-madness dataset, unique vs common: the differences were 79.1% in strain recall, 75.9% in genome fraction and 20.6% in strain precision, with a 50-fold NGA50 gap. Mismatches were 199 (unique) vs 511 (common) per 100 kb. Duplication ratio differed by +0.01 (marine) and -0.08 (strain madness) — [PMC9007738](https://pmc.ncbi.nlm.nih.gov/articles/PMC9007738/)
- Whole-dataset genome fraction for short reads: MEGAHIT reached 41.1% (marine) and 10.4% (strain madness). The gold-standard assemblies reached 76.9% and 90.8%. The hybrid assembler A-STAR reached 44.1% and 30.9% — [PMC9007738](https://pmc.ncbi.nlm.nih.gov/articles/PMC9007738/)
- Best assemblers on common strains:
  - Marine: HipMer ranked best, with strain recall of 4.5% on common genomes vs 20.4% on unique genomes, both at 100% precision. OPERA-MS (hybrid) had the highest NGA50.
  - Strain madness: GATB ranked best. A-STAR had the highest strain recall at only 1.5% (23.1% precision) and the highest NGA50. GATB had the fewest misassemblies.
  - A-STAR ranked first on genome fraction in all four partitions.
  — [PMC9007738](https://pmc.ncbi.nlm.nih.gov/articles/PMC9007738/)
- Unique genomes: GATB, HipMer, MEGAHIT and OPERA-MS reached 100% strain recall and precision on unique strain-madness genomes. SPAdes had the highest NGA50 on unique genomes in both datasets — [PMC9007738](https://pmc.ncbi.nlm.nih.gov/articles/PMC9007738/)
- Binning on gold-standard assemblies: unique strains binned "substantially better" than common strains in both marine and strain madness. Unique-strain bin purity in strain madness was 97.9 ± 0.4%; common-strain values are in the figures and were not extracted. UltraBinner ranked best for unique genomes, CONCOCT ranked best for common strains and across all partitions, Vamb had the highest purity overall, and UltraBinner had the best ARI — [PMC9007738](https://pmc.ncbi.nlm.nih.gov/articles/PMC9007738/)
- Profiling in strain madness: at genus rank, MetaPhlAn 2.9.22 reached 89.2% completeness and 92.8% purity. At species rank, "no method excelled" — [PMC9007738](https://pmc.ncbi.nlm.nih.gov/articles/PMC9007738/)
- Conclusions:
  - "Related strains still were challenging for assembly and genome recovery through binning."
  - Most assemblers "sometimes intentionally, did not resolve strain variation."
  - "Long reads help to distinguish strains", and hybrid assemblers were less affected by closely related strains.
  - For HipMer, single-sample assembly may be preferable "if coverage is sufficient and closely related strains are expected."
  — [PMC9007738](https://pmc.ncbi.nlm.nih.gov/articles/PMC9007738/)

### Inferences
- The 95% ANI threshold that defines "common strains" is coarse. CAMI does not report results binned by finer ANI levels (99%, 99.9%), so CAMI alone cannot locate the exact divergence at which collapse happens.
- Among short-read assemblers, CAMI rewards different behaviours on different metrics. Conservative collapsing (GATB) wins on low misassembly, and aggressive separation (A-STAR) wins on recall at low precision. No method achieves both on strain-rich data.

### Gaps
- No CAMI III benchmark results were found in this session (searches were limited). CAMI III status is unverified.
- Per-binner completeness and purity numbers for common strains in CAMI II are in the figures and supplementary data, which were not extracted.
- The second half of both CAMI papers (~75 kb of CAMI II and ~35 kb of CAMI I) was not read. It may contain more strain-level profiling detail.

## Q2. Strain-level profilers and haplotypers on synthetic or mock data

### Takeaway
Reference-based strain identifiers work at very low depth: StrainScan works down to about 1x and StrainGE claims 0.1x. They fail to separate strains inside a reference cluster at roughly >= 99.9% ANI. StrainScan addresses this with a hierarchical search and gains about 20% in F1 on multi-strain samples. De novo haplotypers (DESMAN, STRONG) need many samples and give only moderate accuracy on strain number.

### Cited Findings
- StrainScan (Liao H, Ji Y, Sun Y, *Microbiome* 11:183, 2023, doi:10.1186/s40168-023-01615-w):
  - Single-strain benchmark: 1,400 simulated datasets at 10x, 5x, 3x and 1x depth. F1 was near-perfect above 1x and dropped slightly at 1x for *A. muciniphila* and *M. tuberculosis*. StrainEst could not identify strains below 5x — [PMC10433603](https://pmc.ncbi.nlm.nih.gov/articles/PMC10433603/)
  - Multi-strain benchmark: 360 datasets with 2, 3 or 5 strains at 10-100x. StrainScan improved F1 by >20% over the other tools. StrainGE's cluster-level recall on 5-strain *A. muciniphila* data fell from 0.94 to 0.42. For 2 *C. acnes* strains at 10x/1x, F1 was 0.98 (StrainScan), 0.88 (KrakenUniq) and 0.81 (StrainGE) — [PMC10433603](https://pmc.ncbi.nlm.nih.gov/articles/PMC10433603/)
  - The default cluster threshold (k-mer Jaccard >= 0.95) corresponds to about 99.89% ANI. Within-cluster strain identity is the case that StrainGE and StrainEst cannot resolve, because they return only the cluster representative — [PMC10433603](https://pmc.ncbi.nlm.nih.gov/articles/PMC10433603/); [search summary of StrainScan preprint](https://www.biorxiv.org/content/10.1101/2022.01.26.477962.full.pdf)
  - Strains missing from the database (90 new *E. coli*, 10x): StrainScan F1 was 95% vs 78% for StrainGE on the 2-strain test — [PMC10433603](https://pmc.ncbi.nlm.nih.gov/articles/PMC10433603/)
  - The comparison set was KrakenUniq, StrainSeeker, Pathoscope2, Sigma, StrainGE and StrainEst. inStrain and StrainPhlAn were not benchmarked — [PMC10433603](https://pmc.ncbi.nlm.nih.gov/articles/PMC10433603/)
  - The StrainScan paper notes that many strain tools need >10x strain coverage, and cites StrainGE's claim of strain identification at 0.1x and variant detection at 0.5x — [StrainScan via search](https://microbiomejournal.biomedcentral.com/counter/pdf/10.1186/s40168-023-01615-w)
- DESMAN (Quince et al., *Genome Biol* 2017, per search snippet): the synthetic test used 64 samples, 5 *E. coli* strains and 15 other strains. Inferred strain frequencies matched relative coverage with R² = 0.9998, and gene presence/absence was 95.7% accurate — [search result, bioRxiv 073825](https://www.biorxiv.org/content/10.1101/073825.full.pdf). Venue and year are from memory and not re-verified this session.
- STRONG (Quince et al., *Genome Biol* 2021, doi:10.1186/s13059-021-02419-7): on high-diversity synthetic sets, the correlation between actual and predicted strain number per MAG was r = 0.52-0.62 (compared against DESMAN and mixtureS). This comes from a search snippet; the full text was blocked — [Genome Biology fig. 3](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02419-7/figures/3)
- Marker-gene tools such as StrainPhlAn do not untangle SNVs from strain mixtures; they report the dominant strain. Assembly-based tools need higher coverage than low-abundance members usually have. This is a search-snippet characterization from the StrainFacts and related literature — [escholarship StrainFacts PDF](https://escholarship.org/content/qt99s137n5/qt99s137n5.pdf?t=rjachl)

### Inferences
- There are two regimes:
  - Reference-based classifiers (StrainScan, StrainGE, StrainEst) work at 0.1-1x but are limited by database clustering at about 99.9% ANI.
  - De novo SNV haplotypers (DESMAN, STRONG, StrainFacts) need multi-sample designs and, from prior literature, several-fold per-strain coverage.

### Gaps
- No quantitative benchmark of inStrain, StrainPhlAn 4 or StrainFacts against each other on ANI-controlled mixtures was retrieved.
- From memory and unverified: inStrain uses a popANI threshold of 99.999% to call "same strain" and a minimum coverage of about 5x for SNV calling. StrainFacts (Smith et al., Front Bioinform 2022) reports accuracy comparable to DESMAN/StrainFinder at much greater speed. Check these against the primary papers.
- STRONG haplotype SNV error rates and minimum coverage were not extracted (Springer redirect).

## Q3. Long-read assemblers with strain mixes (metaFlye, hifiasm-meta, metaMDBG, Strainy, Strainberry, myloasm)

### Takeaway
Standard long-read metagenome assemblers collapse strains at about 99% ANI or higher. On the Zymo D6331 gut mock (5 *E. coli* strains), no assembler resolves all of them. Phasing add-ons separate strains well:
- Strainberry needs >= 0.5% divergence (99.5% ANI) and about 20x per strain, plateauing at 30x.
- Strainy reports better results than Nanopore-based alternatives (numbers not retrieved).

myloasm (2025 preprint) targets this gap by building the graph from polymorphic k-mers.

### Cited Findings
**Strainberry** (Vicedomini R, Quince C, Darling AE, Chikhi R, *Nat Commun* 2021, doi:10.1038/s41467-021-24515-9; full PDF read) — [UEA eprint PDF](https://ueaeprints.uea.ac.uk/88125/1/Vicedomini_etal_2021_NatureComms.pdf)
- Mock3 (*E. coli* K-12 + W, 98.65% identity to each other, plus *B. cereus*):
  - Flye collapsed the two *E. coli* into two consensus sequences (2 Mbp and 2.8 Mbp).
  - Canu had imperfect identity (99.64% and 99.74%) and a high duplication ratio.
  - After Strainberry, both strains were near-completely covered at 99.97% ANI (ssFlye).
  — [PDF](https://ueaeprints.uea.ac.uk/88125/1/Vicedomini_etal_2021_NatureComms.pdf)
- Mock9 (includes two *S. aureus* and two *E. coli*): Canu's *S. aureus* identity was 99.86% and 99.74%. Strainberry separation gave near-complete strains at 99.9% identity — [PDF](https://ueaeprints.uea.ac.uk/88125/1/Vicedomini_etal_2021_NatureComms.pdf)
- Coverage series (5x-50x): quality plateaus after 30x and is "mediocre at 10x". More than 95% reference representation and more than 99.8% identity require at least 20x. Uneven coverage between strains did not matter; only absolute depth did — [PDF](https://ueaeprints.uea.ac.uk/88125/1/Vicedomini_etal_2021_NatureComms.pdf)
- Divergence series (*E. coli* K-12 vs 7 strains, simulated PacBio): Strainberry fails at low divergence. Reference coverage improves significantly at 0.39% divergence, and assembly is near-complete at 0.50% divergence (about 99.5% ANI) — [PDF](https://ueaeprints.uea.ac.uk/88125/1/Vicedomini_etal_2021_NatureComms.pdf)
- Number of strains (2-5 *E. coli*, pairwise divergence 0.7-1.4%): average reference coverage fell from 95% (3 strains) to 75% (5 strains), and ANI fell from 99.9% to 99.6% — [PDF](https://ueaeprints.uea.ac.uk/88125/1/Vicedomini_etal_2021_NatureComms.pdf)
- Recombination: Flye gave strain-oblivious assemblies for species with average and high recombination rates (*E. coli*, *H. pylori*, *N. meningitidis*). Strainberry separated them, with slightly lower ANI — [PDF](https://ueaeprints.uea.ac.uk/88125/1/Vicedomini_etal_2021_NatureComms.pdf)
- Overall: "an average depth of coverage of 60-80x suffices" for near-complete strains at >99.9% identity, and "any strain present at coverage above 10x ends up being correctly separated". PacBio gives better identity and Nanopore better contiguity. On Nanopore the identity of the *Lactobacillus* strains was 98.78-99.76%, vs 99.69-99.97% with PacBio — [PDF](https://ueaeprints.uea.ac.uk/88125/1/Vicedomini_etal_2021_NatureComms.pdf)

**Strainy** (Kazantseva E, Donmez A, Frolova M, Pop M, Kolmogorov M, *Nat Methods* 21(11):2034-2043, 2024, doi:10.1038/s41592-024-02424-1; PMID 39327484; metadata verified via Europe PMC)
- The abstract states that Strainy "outperform[s] current Nanopore-based methods and [is] comparable with PacBio-based algorithms in completeness and accuracy" on simulated and mock data. Long-read assemblers "suppress strain-level variation in favor of a species-level consensus". The abstract contains no numbers — [Europe PMC record](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:39327484%20AND%20SRC:MED&resultType=core&format=json)
- The preprint benchmarked against metaFlye, Strainberry and iGDA on mock bacterial communities — [bioRxiv 2023.01.31.526521](https://www.biorxiv.org/content/10.1101/2023.01.31.526521.full.pdf)

**Zymo D6331 gut mock (5 *E. coli* strains), HiFi:**
- The figures below come from a search summary. Attribution to specific papers is inferred from the URLs: PMC9343089 is probably hifiasm-meta (Feng et al., *Nat Methods* 2022) and PMC10541625 is probably metaMDBG (Benoit et al., *Nat Biotechnol* 2024). Verify both.
  - hifiasm-meta assembled *E. coli* B766 as one circular contig, B3008 as 2 contigs and the other strains as fragments. HiCanu circularized both B766 and B3008. metaFlye circularized none of the 5 — [search result, PMC9343089](https://pmc.ncbi.nlm.nih.gov/articles/PMC9343089)
  - The same source says high-divergence strains separate into disconnected contigs, while low-divergence strains collapse into mosaic contigs — [PMC9343089](https://pmc.ncbi.nlm.nih.gov/articles/PMC9343089)
  - metaMDBG and hifiasm-meta each circularized one *E. coli* strain, and no assembler resolved all of the *E. coli* diversity. Across the Zymo mock, metaMDBG and hifiasm-meta each produced 10 circular genomes and metaFlye 9. metaMDBG also produced 2 near-complete (>99.8%) linear genomes. The metaMDBG authors note that sufficiently similar strains may merge into a consensus because of low minimizer density — [search result, PMC10541625](https://pmc.ncbi.nlm.nih.gov/articles/PMC10541625)
- Trigodet et al., *Nat Biotechnol* 2026 (per search snippet; the AWI PDF was not opened) compared HiCanu, hifiasm-meta, metaFlye and metaMDBG on 21 HiFi metagenomes. metaFlye and metaMDBG produced up to 180x more read-clipping events (a proxy for assembly errors) than HiCanu and hifiasm-meta. metaMDBG produced about 2x more circular contigs than hifiasm-meta, which the authors flag as possible premature circularization — [AWI PDF](https://epic.awi.de/60560/1/Trigodet_2026_nature_biotechnology.pdf)

**myloasm** (Shaw J, Marin MG, Li H, bioRxiv 2025, doi:10.1101/2025.09.05.674543; preprint, not peer-reviewed as of indexing)
- myloasm builds the string graph from polymorphic k-mers ("SNPmers") and simplifies it by differential abundance.
- On real ONT metagenomes it reports 3x more complete circular contigs than the next-best assembler.
- It recovered 6 complete single-contig *Prevotella copri* genomes from one gut sample and 8 complete TM7 contigs with >93% similarity from an oral sample.
- It targets HiFi and ONT R10.4. The docs list v0.4.0 (Jan 2026) and say it uses more RAM than metaMDBG.
— [bioRxiv PDF](https://www.biorxiv.org/content/10.1101/2025.09.05.674543.full.pdf); [myloasm docs](https://myloasm-docs.github.io); [results page](https://myloasm-docs.github.io/results)

### Inferences
- Strainberry's divergence series gives the best-documented threshold:
  - Long-read phasing separates strains down to about 0.4-0.5% divergence (99.5-99.6% ANI) at >= 20-30x per strain.
  - Below that, default settings leave strains collapsed.
  - Plain Flye or metaFlye collapses even pairs at 98.65% ANI (*E. coli* K-12 vs W).
- The 5-strain *E. coli* Zymo mock is the de facto test of strain collapse for HiFi assemblers. No assembler reported in these sources resolves all 5.

### Gaps
- No Strainy numbers (genome fraction, NGA50, phasing errors, ANI and coverage limits) were retrieved; the bioRxiv full text and PMC were rate-limited or blocked. Get them from the paper's Fig. 2-3.
- No myloasm benchmarks on synthetic ANI-controlled strain mixes were extracted. The myloasm results page may have them.
- No head-to-head of Strainy vs Strainberry vs myloasm on the same ANI ladder was found.
- This session found no long-read strain benchmark on ATCC MSA-1003 / MSA-2006.

## Q4. Binners on strain-rich data and strain-aware binning

### Takeaway
No binner separates near-identical genomes well, and no published benchmark was found that sweeps binners across controlled ANI (e.g. 97%, 99%, 99.9%). In CAMI II, CONCOCT ranked best on common strains and Vamb had the highest overall purity. A 2025 evaluation finds SemiBin2 and COMEBin best overall. On a strain-associated CAMI2 set, embedding accuracy was only 36-52%.

### Cited Findings
- CAMI I: MaxBin 2.0 alone kept >90% purity on common strains. Completeness collapsed for all binners — [PMC5903868](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903868/)
- CAMI II: CONCOCT ranked best for common strains and across all partitions, UltraBinner for unique genomes, and Vamb had the highest purity — [PMC9007738](https://pmc.ncbi.nlm.nih.gov/articles/PMC9007738/)
- "Evaluation of metagenome binning: advances and challenges" (bioRxiv 2025.02.21.639465; later in PMC12636519, venue not verified because the PMC page was blocked):
  - On a CAMI2 strain-associated dataset, embedding accuracy was McDevol 52.2%, GenomeFace 46.1%, COMEBin 41.6% and Vamb 35.9%.
  - SemiBin2 and COMEBin gave the best overall binning, and MetaBAT2 and GenomeFace were the fastest.
  - These figures come from a search snippet.
  — [bioRxiv PDF](https://www.biorxiv.org/content/10.1101/2025.02.21.639465.full.pdf); [PMC12636519](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12636519/)

### Inferences
- Binners rely on coverage profiles and tetranucleotide frequency. Strains at >= 99% ANI have essentially identical composition, so only differential coverage across samples can separate them. This matches CAMI I's remark that many more replicate samples are needed. A practical design is to bin at species level and resolve strains afterwards (STRONG, Strainy, Strainberry, inStrain).

### Gaps
- No dedicated "strain-aware binning" method with benchmark numbers was found in this session.
- No VAMB, SemiBin2, COMEBin or MetaBAT2 benchmark on ANI-stratified strain mixes was found. A search assistant said none exists for 99% ANI; that is plausible but not exhaustively verified.

## Q5. General thresholds (ANI collapse, depth for phasing)

### Takeaway
Reported thresholds:

| Approach | Behaviour |
|---|---|
| Short-read assemblers (CAMI) | Fail broadly once a genome has a relative at >= 95% ANI |
| Standard long-read assemblers | Collapse pairs at about 98.6% ANI or higher (Flye on *E. coli* K-12/W) |
| Long-read phasing (Strainberry) | Separates down to about 99.5% ANI; needs >= 20x per strain, best at >= 30x |
| Reference-based identifiers | Detect strains at 0.1-1x but merge strains above about 99.9% ANI unless they search within clusters (StrainScan) |

### Cited Findings
- Common-strain threshold: 95% ANI, with short-read genome fraction dropping from about 98% to <= 22.5% (CAMI I) and a 76-point drop in CAMI II strain madness — [PMC5903868](https://pmc.ncbi.nlm.nih.gov/articles/PMC5903868/); [PMC9007738](https://pmc.ncbi.nlm.nih.gov/articles/PMC9007738/)
- Flye collapses *E. coli* K-12 and W at 98.65% identity — [Strainberry PDF](https://ueaeprints.uea.ac.uk/88125/1/Vicedomini_etal_2021_NatureComms.pdf)
- Long-read phasing limits: 0.39-0.50% divergence, >= 20x per strain (plateau at 30x), and up to 5 strains with loss of quality — [Strainberry PDF](https://ueaeprints.uea.ac.uk/88125/1/Vicedomini_etal_2021_NatureComms.pdf)
- Reference-based identification: StrainScan works down to 1x and StrainEst needs >= 5x. StrainGE claims 0.1x for identification and 0.5x for variants. The cluster threshold is about 99.89% ANI — [PMC10433603](https://pmc.ncbi.nlm.nih.gov/articles/PMC10433603/)

### Inferences
- There are no clean published "fail at X% ANI with short reads" numbers finer than CAMI's 95% split.
- From these results, short-read de Bruijn assemblers break up or collapse genomes somewhere between 95% and 99%+ ANI. Long-read phasing pushes the separable range to about 99.5% ANI at 20-30x. Above 99.9% ANI, only reference-based or SNV-haplotyping approaches distinguish strains.

### Gaps
- No dedicated short-read ANI-ladder benchmark (e.g. assembly at 97/98/99/99.5/99.9% ANI) was retrieved. Such studies exist in prior literature (e.g. Quince, Olm), but no numbers were obtained in this session.
- The minimum coverage for short-read SNV phasing (DESMAN, STRONG, inStrain) was not extracted from primary sources.
