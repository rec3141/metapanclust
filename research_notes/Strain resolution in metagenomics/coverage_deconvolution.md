# Coverage-based deconvolution of strain and species mixtures across multiple metagenome samples (2012-2026)

Scope note: the primary full texts of DESMAN, STRONG, StrainPanDA, MSPminer and StrainFacts were read via Europe PMC full-text XML (PMC5607848, PMC8311964, PMC10989911, PMC6499236, PMC9580935). Other items come from abstracts or search snippets and are marked as such. Items marked [UNVERIFIED] come from background knowledge and were not checked against a source in this session.

## 1. Formal models and the tools that implement them

### Takeaway
Every coverage-deconvolution tool fits the same linear mixing model, `X (features x samples) ≈ H (features x strains, 0/1 or copy number) · Γ (strains x samples, abundance)`. They differ in what a "feature" is: SNV allele counts (DESMAN, StrainFacts), assembly-graph unitigs (STRONG/BayesPaths), or gene/gene-family depth (DESMAN's accessory step, StrainPanDA, MSPminer, MIDAS/StrainPGC). They also differ in noise model and inference: NMF/NTF initialisation followed by a Gibbs sampler, variational Bayes, MAP gradient descent, or sparse NMF. Accessory-gene-to-strain assignment works in two stages in every tool: (a) estimate strain abundances Γ from features known to be single-copy core, then (b) hold Γ fixed and fit a binary presence vector per gene or contig.

### Cited Findings
**DESMAN** (Quince, Delmont, Raguideau, Alneberg, Darling, Collins, Eren; *Genome Biology* 18:181, 2017; verified) — [PMC5607848](https://pmc.ncbi.nlm.nih.gov/articles/PMC5607848); [GitHub](https://github.com/chrisquince/DESMAN)
- Input: a co-assembly whose contigs have been binned, with CONCOCT used as the binner. Variants in core genes are linked by their co-occurrence across samples into haplotypes and abundance profiles. "These are then searched for against non-core genes to determine the accessory genome of each strain." — [PMC5607848](https://pmc.ncbi.nlm.nih.gov/articles/PMC5607848)
- Core definition: DESMAN uses genes "present in all strains as a single copy". With no references these are the 36 single-copy core COGs of Alneberg et al. With reference genomes, it also uses "single-copy core species genes (SCSGs)". Variant positions are called by a likelihood-ratio test on base frequencies summed across samples. — [PMC5607848](https://pmc.ncbi.nlm.nih.gov/articles/PMC5607848)
- Haplotype step: a non-negative tensor factorisation (multiplicative NMF updates on a 4V x S stacked matrix under KL divergence) gives the initial solution. A Gibbs sampler then refines it. The NMF is justified because observed base proportions are treated as "a sum of independent Poisson-distributed components from each haplotype". — [PMC5607848](https://pmc.ncbi.nlm.nih.gov/articles/PMC5607848)
- Accessory-gene model (the exact "coverage = Σ abundance × indicator" model):
  - Per-strain coverage in sample s is γ_{g,s} = π_{g,s}·Y_s, where Y_s is the mean total coverage on core genes.
  - Gene coverage x_{f,s} ~ Poisson(λ_{f,s}) with λ_{f,s} = Σ_g η_{f,g} γ_{g,s}, where η_{f,g} ∈ {0,1} is the presence of gene f in strain g. The authors note the model "could be extended to multi-copy genes" with a geometric prior penalising copies.
  - This coverage likelihood is multiplied by a multinomial variant likelihood. The variant likelihood renormalises strain frequencies over only those strains that carry the gene, π′ = π_{g,s}/Σ_{g∈G_f} π_{g,s}.
  - η is Gibbs-sampled gene by gene with π and the error matrix fixed. — [PMC5607848](https://pmc.ncbi.nlm.nih.gov/articles/PMC5607848)
- The binary assumption is reasonable because 98.8% of genes in the STEC genome are present in zero or one copy. The authors also say the strategy "could equally well be applied to either contigs or genes". — [PMC5607848](https://pmc.ncbi.nlm.nih.gov/articles/PMC5607848)
- Choosing the number of strains G: G is chosen to maximise the number of haplotypes with mean SNV uncertainty below 10% and mean relative abundance above 5%. — [PMC5607848](https://pmc.ncbi.nlm.nih.gov/articles/PMC5607848)
- Performance:
  - Mock community (50 species, 210 genomes, 96 samples): 75 haplotypes were inferred, of which 3 were false positives. Of the 67 from true multi-strain clusters, 50.7% were exact and 79.1% were within 5 SNVs.
  - Gene presence/absence accuracy was 94.9% on average, and 97.39% for haplotypes with SNV error below 1%. Gene-assignment error tracks haplotype error (adjusted R² = 0.697).
  - Predicted versus true strain frequency regression: slope 0.82, i.e. a bias towards underestimation. — [PMC5607848](https://pmc.ncbi.nlm.nih.gov/articles/PMC5607848)
- E. coli STEC outbreak data (Loman et al., 20 samples): eight strains were inferred across 8,566 genes. The outbreak-matching strain had 91.8% of its gene calls identical to the O104:H4 reference. — [PMC5607848](https://pmc.ncbi.nlm.nih.gov/articles/PMC5607848)

**STRONG / BayesPaths** (Quince, Nurk, Raguideau, James, Soyer, Summers, Limasset, Eren, Chikhi, Darling; *Genome Biology* 22:214, 2021; verified) — [PMC8311964](https://pmc.ncbi.nlm.nih.gov/articles/PMC8311964); [Genome Biology record](https://doaj.org/article/c74ea889eddf4483a41b363af8e7877a)
- STRONG keeps the co-assembly graph before bubble and variant simplification, extracts the single-copy core gene (SCG) subgraphs for each MAG, and decomposes their unitig coverages across samples into strain paths. — [PMC8311964](https://pmc.ncbi.nlm.nih.gov/articles/PMC8311964)
- Model: λ_{v,s} = θ_v Σ_g γ_{g,s} φ_{v,g}. Here φ_{v,g} indicates that strain g's path traverses unitig v, γ_{g,s} is strain coverage in sample s, and θ_v is a per-node bias.
  - Noise is Gaussian with empirically smoothed per-node precision, using a GAM on log x. This is analogous to voom, because "the true distribution of read counts is Poisson".
  - Inference is coordinate-ascent variational inference, initialised by variational Bayesian NMF.
  - The number of strains is chosen by 10-fold cross-validation of the predictive likelihood, with weights τ′ = 1/(0.5 + x) to mimic Poisson noise. — [PMC8311964](https://pmc.ncbi.nlm.nih.gov/articles/PMC8311964)
- Resolution limits: with k-mers of about 100, binning "best case" resolves about 1% divergence. Strains with divergence close to the reciprocal of read length (about 98-99.5% identity for 75-150 bp reads) produce very tangled graphs. Shared contigs in a bin become "a consensus of the strains present". — [PMC8311964](https://pmc.ncbi.nlm.nih.gov/articles/PMC8311964)
- STRONG resolves strains on SCGs only. It does not assign accessory genes to strains in the published version. This follows from the paper's description of output on single-copy core genes. — [PMC8311964](https://pmc.ncbi.nlm.nih.gov/articles/PMC8311964)

**StrainPanDA** (Hu, Tan, Li, … Dai; *iMeta* 2022; verified as iMeta) — [PMC10989911](https://pmc.ncbi.nlm.nih.gov/articles/PMC10989911); [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.02.15.480535.full.pdf)
- Model: D = P·S. D is the N gene families x S samples normalised count matrix, P is the binary N x K gene-content matrix, and S is the K x S strain composition matrix with columns summing to 1.
- Fitting uses sparse NMF ("snmf/r", R package NMF 0.21), and "the addition of sparsity constraints … ensures the uniqueness of factorization". P′ is binarised with a threshold θ and a confidence score.
- Rank K is the smallest value in 1-12 that satisfies three conditions: every strain has mean abundance above 0.1, every strain has more than 0.5·g_min gene families, and every pair of strains has Jaccard distance above 0.1 in gene content. — [PMC10989911](https://pmc.ncbi.nlm.nih.gov/articles/PMC10989911)
- StrainPanDA is reference-pangenome based (PanPhlAn-style mapping to a species pangenome), not de novo. — [PMC10989911](https://pmc.ncbi.nlm.nih.gov/articles/PMC10989911)

**StrainFacts** (Smith, Li, Shi, Abate, Pollard; *Frontiers in Bioinformatics* 2:867386, 2022; verified) — [PMC9580935](https://pmc.ncbi.nlm.nih.gov/articles/PMC9580935); [Frontiers](https://www.frontiersin.org/articles/10.3389/fbinf.2022.867386)
- Model: p_{ig} = Σ_s γ_{sg} π_{is}, i.e. P = ΓΠ, where p is the allele frequency at SNP g in sample i. The authors "highlight parallels to non-negative matrix factorization".
- Genotypes are relaxed to fuzzy (continuous) values so gradient MAP fitting works. Counts are modelled as BetaBinomial, i.e. overdispersed.
- The number of strains was set to 30% of the number of samples, or 200 for 20,550 metagenomes. The method scales to tens of thousands of samples on GPU and is about 100x faster than earlier tools. — [PMC9580935](https://pmc.ncbi.nlm.nih.gov/articles/PMC9580935)
- StrainFacts works on SNVs (GT-Pro metagenotypes) and does not assign accessory genes. — [PMC9580935](https://pmc.ncbi.nlm.nih.gov/articles/PMC9580935)

**Co-abundance gene binning (species level, with accessory genes)**
- **Nielsen et al.** (*Nature Biotechnology* 32:822-828, 2014; verified bibliographically).
  - Canopy clustering of gene abundance profiles across 396 gut samples produced 7,381 co-abundance gene groups (CAGs), of which 741 were metagenomic species (MGS), and 238 genomes were assembled from them.
  - The search summary says species could be separated with "as few as 18 samples". This was not verified in the full text. — [Nature Biotech PDF](https://www.nature.com/nbt/journal/v32/n8/pdf/nbt.2939.pdf); [DTU dataset page](https://services.healthtech.dtu.dk/datasets/MetaGenomicSpecies)
- **MSPminer** (Plaza Oñate, Le Chatelier, Almeida, Cervino, Gauthier, Magoulès, Ehrlich, Pichaud; *Bioinformatics* 35(9):1544-1552, 2019, online 2018; verified).
  - Bins co-abundant genes using "a new robust measure of proportionality coupled with an empirical classifier to group and distinguish not only species core genes but accessory genes also". Applied to 9.9 M genes x 1,267 stool samples, it produced 1,661 MSPs.
  - Default input filter: genes need more than 6 reads in at least 3 samples.
  - Simulations reduced the number of samples carrying the species (200, 100, 50) to test completeness. They also simulated two-strain samples with the dominant strain 5-10x more abundant than the subdominant one.
  - The paper notes that earlier pipelines grouped core genes with highly prevalent accessory genes and missed low-prevalence accessory genes. Linking accessory to core clusters had relied on Fisher's exact test on presence/absence (Nielsen 2014), which "does not account for the co-abundance of genes". — [PMC6499236](https://pmc.ncbi.nlm.nih.gov/articles/PMC6499236)

**Reference-pangenome depth-ratio methods (accessory genes against core/species depth)**
- **MIDAS** (Nayfach, Rodriguez-Mueller, Garud, Pollard; *Genome Research* 2016).
  - Species coverage comes from universal single-copy marker genes, then pangenome genes give gene coverage and copy number.
  - A figure caption states that presence/absence "is accurately predicted when genome coverage is above 1×, and a gene copy number cutoff of 0.35 is used". — [PMC5088602 via search](https://pmc.ncbi.nlm.nih.gov/articles/PMC5088602)
- **MIDAS v3 + StrainPGC** (Smith et al., *Genome Research* 35:1247-1260, 2025; abstract-level verification).
  - Strain tracking first splits samples into strain-pure subsets. A gene is then called present in a strain if it passes thresholds on both the depth ratio (gene depth / species depth) and a correlation score.
  - The paper lists a "decreased depth ratio for strain-specific genes when other strains of the same species are also abundant" as a key failure mode of depth ratio alone. — [Genome Research](https://genome.cshlp.org/content/35/5/1247); [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.04.10.588779.full.pdf)
- PanPhlAn (Scholz et al., *Nature Methods* 2016; PanPhlAn 3 in Beghini et al., *eLife* 2021) calls gene presence from gene coverage relative to the species' median coverage, dominant strain only. [UNVERIFIED in this session; StrainPanDA's own text describes earlier pangenome methods as able to infer only the dominant strain — [PMC10989911](https://pmc.ncbi.nlm.nih.gov/articles/PMC10989911)]

**Contig coverage binners and hybrid assembly**
- **CONCOCT** (Alneberg et al., *Nature Methods* 11:1144-1146, 2014; verified bibliographically) clusters contigs on composition plus coverage across multiple samples. — [PubMed 25218180](https://pubmed.ncbi.nlm.nih.gov/25218180/); [arXiv 1312.4038](https://arxiv.org/abs/1312.4038)
- **OPERA-MS** (Bertrand et al., *Nature Biotechnology* 37:937-944, 2019; verified bibliographically).
  - Bayesian clustering of contigs on coverage plus long-read contiguity, giving "strain-resolved assembly in the presence of multiple genomes of the same species".
  - It is a single-sample method: strains are separated by coverage within one sample, not across samples.
  - Recommends more than 15x short-read coverage. Long reads from about 9x give near-complete genomes for rare species. — [Nature Biotech](https://www.nature.com/articles/s41587-019-0191-2); [GitHub](https://github.com/CSB5/OPERA-MS)
- MetaGEM is a workflow wrapper (CONCOCT/MetaBAT2/MaxBin2 plus SMETANA). I found no evidence that it adds a strain-deconvolution model. [UNVERIFIED; not searched]

### Inferences
- The model the project wants (contig coverage vector = Σ_strains abundance vector x indicator) is exactly DESMAN equation 7, λ_{f,s} = Σ_g η_{f,g} γ_{g,s}, applied to contigs instead of genes. DESMAN states that the same method applies to contigs. On a co-assembly with near-identical strains, the practical recipe is:
  1. Estimate γ (strains x samples) from SNVs on SCGs or SCSGs, or from graph bubbles (STRONG).
  2. Fit a binary η per accessory contig by Poisson or negative-binomial likelihood with γ fixed.
- When strains are too similar for SNV phasing but differ in gene content, StrainPanDA's sparse-NMF formulation on gene or contig depth is the closest de novo analogue. It needs a pangenome; a co-assembly's contig catalogue can play that role.
- Gene-content-based tools (StrainPanDA, MSPminer) cannot separate strains with identical accessory content. SNV-based tools (DESMAN, StrainFacts) cannot separate strains identical on core genes. The two signals are complementary.

### Gaps
- I found no metagenomics tool that is explicitly framed as archetypal analysis or convex-hull vertex finding on contig coverage. The geometric machinery exists in the NMF literature (see section 4), but no metagenomic application turned up.
- The CONCOCT and MetaBAT2 papers' sample-number-versus-binning-accuracy curves were not retrieved (paywalled). No numbers are given here.
- PanPhlAn thresholds, ConStrains (Luo et al., *Nat Biotech* 2015), Lineage (O'Brien et al., 2014) and mixtureS were not checked. PanSN is a sequence-naming convention for pangenome graphs, not a deconvolution method [UNVERIFIED].

## 2. Identifiability: samples, abundance variation, noise, detection limits

### Takeaway
Empirical requirements from benchmarks:
- **Samples:** about 10 or more, with performance still improving up to roughly 30.
- **Coverage:** about 5x or more per strain. Above 1x genome coverage is enough for gene presence calls on a dominant strain.
- **Abundance:** species relative abundance of at least 1%; minor strains of about 5-10% or more.
- **Divergence:** for SNV methods, at least about 0.05-0.4% divergence on the SCGs.

Constant strain ratios across samples are the core non-identifiable case. Coverage alone cannot separate them, and SNV methods lose accuracy as ratios stabilise. In theory, NMF is unique only under separability or "sufficient scatteredness" of the abundance vectors.

### Cited Findings
- **DESMAN, samples.**
  - Sample subsets of size 1-64 were drawn from the 64-sample E. coli mock: "The SNV error rate starts to increase when the sample number is below about 30; however, the average error is still around 15%, even with just ten samples."
  - "At low sample number … typically some of the strains are resolved accurately and others are missed completely." — [PMC5607848](https://pmc.ncbi.nlm.nih.gov/articles/PMC5607848)
- **DESMAN, divergence.**
  - Strains were resolvable at more than 100 SNPs, which is about 0.38% divergence over the SCGs. Detection sometimes worked with "as few as ten SNVs".
  - The Gibbs sampler for η converged after about 20 iterations using 20 random variant positions. — [PMC5607848](https://pmc.ncbi.nlm.nih.gov/articles/PMC5607848)
- **STRONG, samples.**
  - Synthetic designs used 3, 5, 10 and 15 samples (Synth_S03 to S15).
  - "Performance does improve with increasing sample number … even when only a small number of samples are available we are still able to accurately predict (with 0.068% per base error rate) strains, and when ten or more samples are available we obtain error rates below 0.05%."
  - Across all data sets STRONG found 216 strains versus 208 for DESMAN. — [PMC8311964](https://pmc.ncbi.nlm.nih.gov/articles/PMC8311964)
- **STRONG, divergence and coverage.**
  - Synthetic strains were constrained to 0.05-5% SCG divergence, i.e. at least about 15 nt over about 30 kbp of SCGs.
  - Of the 280 strains, 64 (22.8%) were missed. Only 7 of those 64 were below the minimum coverage of detected strains (5.68x), so most were missed because they were "not sufficiently divergent in terms of nucleotides or coverage profiles".
  - Graph edges below max(2.0, 1% of MAG coverage) are removed during simplification. — [PMC8311964](https://pmc.ncbi.nlm.nih.gov/articles/PMC8311964)
- **STRONG, constant abundance ratios.**
  - Dirichlet concentration θ was varied while keeping 10 samples. The coefficient of variation of strain proportions fell from 1.45 to 0.65 to 0.23.
  - Strain number was still found correctly, but haplotype error rose from "practically zero" to 0.086% to 0.38% for the most stable configuration.
  - The authors caution about "highly diverse or stable MAGs". In a real anaerobic-digester time series (10 samples, about 5-weekly), about 40% of MAGs showed strain variation, but it was "relatively stable compared to the MAG dynamics". — [PMC8311964](https://pmc.ncbi.nlm.nih.gov/articles/PMC8311964)
- **StrainPanDA, samples and depth.**
  - Accuracy "improves with sample size … and sequencing depth". The authors recommend "at least 10 samples and the relative abundance of the species of interest to be above 1%".
  - Hard cases are "genetic mosaic or highly similar gene content profiles (i.e., lacking strain-unique features)", so the tool is "most suitable for … subspecies".
  - Simulations used a minimum strain frequency of 5% and samples at 20x total species depth. Sample-size runs used 40 samples downsampled to 20, 15, 10 and 5. — [PMC10989911](https://pmc.ncbi.nlm.nih.gov/articles/PMC10989911)
- **StrainFacts, simulation and filters.**
  - Simulations used sequencing depth m = 10 per site, error 0.01, and BetaBinomial noise to model overdispersion.
  - Real-data filters drop SNPs with less than 5% minor-allele occurrence and samples with less than 5% horizontal coverage. — [PMC9580935](https://pmc.ncbi.nlm.nih.gov/articles/PMC9580935)
- **MIDAS, depth threshold.** Gene presence/absence is accurate above 1x genome coverage with a copy-number cutoff of 0.35 (figure caption via search). — [PMC5088602](https://pmc.ncbi.nlm.nih.gov/articles/PMC5088602)
- **MSPminer, sample threshold.** Genes need at least 3 samples with counts above 6. Completeness was tested as species prevalence dropped from 200 to 50 samples. — [PMC6499236](https://pmc.ncbi.nlm.nih.gov/articles/PMC6499236)
- **Nielsen 2014, samples.** Species separation was reported possible with about 18 samples (search summary only). — [Nature Biotech PDF](https://www.nature.com/nbt/journal/v32/n8/pdf/nbt.2939.pdf)
- **NMF identifiability theory.**
  - Separability requires one factor to contain a scaled permutation submatrix. For strains, that means each strain is alone in some sample, i.e. there are "anchor" samples.
  - The sufficiently scattered condition (SSC) only needs the cone of the rows to contain the second-order (ice-cream) cone tangent to the facets of the orthant. For rank ≤ 2, separability and SSC coincide.
  - Fu, Huang & Sidiropoulos, "On identifiability of NMF", *IEEE Signal Processing Letters* 2018, give an identification criterion that needs only one factor to be sufficiently scattered. — [arXiv 1709.00614](https://arxiv.org/abs/1709.00614); review: [arXiv 1803.01257](https://arxiv.org/pdf/1803.01257)
  - StrainPanDA relies on sparsity regularisation for uniqueness. — [PMC10989911](https://pmc.ncbi.nlm.nih.gov/articles/PMC10989911)

### Inferences
- Linear-algebra minimum: separating k strains by coverage needs the k strain abundance vectors (length S) to be linearly independent, so S ≥ k. If all strains keep a constant ratio across samples, rank(Γ) = 1 and core and accessory contigs fall on the same ray. Only magnitude then differs: an accessory contig sits at fraction f = a_g/Σa of the core coverage. You can still see that a contig is in only a subset of the population, but not which strain carries it.
  - With k = 2 strains, a constant ratio r still separates "carried by strain 1" (coverage fraction r/(1+r)) from "strain 2" (1/(1+r)) if r ≠ 1 and noise is small. This is magnitude-based assignment from a single ray, and it fails for r ≈ 1.
- For a unique (non-rotatable) NMF, the abundance vectors need to be "scattered". In practice that means samples where each strain is relatively dominant. Strong anti-correlation, or samples where one strain is absent, is the best case. That geometry explains why time series with stable ratios (STRONG's AD data) are hard.
- Noise: Poisson variance equals the mean, so the coefficient of variation per contig is about 1/√(coverage x length/read length) for read counts. Long contigs average out noise, and gene-level (about 1 kb) features are much noisier than contig-level ones. Empirically, coverage is overdispersed, which is why STRONG uses empirical variances and StrainFacts uses BetaBinomial. A rule of thumb from the benchmarks is that strains below about 5x per sample, or below about 5% relative frequency, are frequently missed.
- Sample-count guidance from the literature: ≥10 samples (StrainPanDA, STRONG); ≤30 samples raises DESMAN error; Nielsen's MGS reportedly needs about 18 or more. All of these assume abundance actually varies across those samples.

### Gaps
- I found no closed-form sample-size or power analysis for coverage-only strain deconvolution, e.g. "S samples at depth d separate k strains with ratio CV c". All quantitative guidance is empirical, from simulations.
- I did not retrieve the StrainPanDA supplementary figures S3/S18 numbers (depth and sample-size curves) or the CONCOCT/MetaBAT sample curves.
- Long-read (nanopore) coverage noise models for this problem were not found.

## 3. Accessory-to-core coverage ratio as population frequency (fraction of cells carrying an element)

### Takeaway
The gene-depth / core-depth ratio is the standard estimator. It is MIDAS's copy number, StrainPGC's "depth ratio", and DESMAN's λ normalised by Y_s. For a single-copy element it estimates the fraction of the species' cells (or coverage) carrying it, Σ_{g carries} π_g. Three things confound it: copy number (plasmids, multicopy MGEs), strain mixing (it lowers the ratio for strain-specific genes), and cross-mapping. Plasmid copy-number studies use the same ratio but read it as copies per cell, not as carriage fraction. The two readings cannot be separated from coverage alone.

### Cited Findings
- In DESMAN, a gene's expected coverage is λ = Σ_g η_g π_{g,s} Y_s, so λ/Y_s = Σ_{g: η=1} π_{g,s}. This is the fraction of species coverage from the strains carrying the gene. Multi-copy genes would scale this by copy number (η > 1). — [PMC5607848](https://pmc.ncbi.nlm.nih.gov/articles/PMC5607848)
- MIDAS estimates gene copy number as gene coverage relative to species marker-gene coverage and calls presence at a ratio of at least 0.35 when genome coverage is above 1x. — [PMC5088602](https://pmc.ncbi.nlm.nih.gov/articles/PMC5088602)
- StrainPGC calls the depth ratio "the key criterion for the assignment of genes to a species". It notes this ratio is decreased for strain-specific genes when other strains of the same species are abundant, and adds a cross-sample correlation score to compensate. — [Genome Research 2025](https://genome.cshlp.org/content/35/5/1247)
- Plasmid copy number (PCN) is the length-normalised plasmid coverage divided by the length-normalised chromosome coverage, giving "copies per the longest chromosome per cell" (isolate-based work found via search). Metagenomes are harder because of a wide spread of chromosome and plasmid coverages. — [PMC12222811 (via search)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12222811); [Genome Research plasmid detection](https://genome.cshlp.org/content/29/6/961)
- Example of extreme plasmid abundance: pBI143 recruited 0.1-3.5% of gut metagenomic reads at a median coverage above 7,000x. A thesis section on "Approximate copy number ratio calculation in metagenomes" exists but was not read. — [Fogarty thesis, UChicago](https://knowledge.uchicago.edu/record/7635/files/Fogarty_uchicago_0330D_16852.pdf) (via search snippet only)
- Abundance-correlation host assignment of resistance genes (ARGs): compare many samples with varying composition and test whether a gene's abundance tracks host taxonomic markers. The markers need not be on the same fragment, and multiple samples are required. — [ES&T/PMC10055365 via search](https://pmc.ncbi.nlm.nih.gov/articles/PMC10055365)

### Inferences
- For a co-assembly bin, the ratio r_{c,s} = cov(contig c, sample s) / median cov(core contigs, sample s) estimates the carrier fraction f_{c,s}, assuming single copy.
  - If r is constant across samples while the strain ratios vary, the contig is either core with a constant copy number, or carried by all strains with ratio f (i.e. multicopy or partially present).
  - If r tracks one strain's π_{g,s} across samples, the contig belongs to that strain. This "ratio follows a strain profile" test is the depth-ratio-plus-correlation logic of StrainPGC and DESMAN.
- Ratios above 1 imply multi-copy elements (plasmids, IS elements, rRNA operons). Values between 0 and 1 are ambiguous between "carried by a subset at 1 copy" and copy number below 1, which is biologically uncommon except under plasmid loss or heterogeneity. A plasmid with PCN > 1 carried by a subset cannot be resolved from coverage alone. Hi-C, long reads, or the multi-sample covariation of the ratio with a strain profile are needed.

### Gaps
- I found no dedicated paper that formally estimates the "fraction of host cells carrying an MGE" from multi-sample metagenome coverage ratios, as opposed to copy number. Gene-frequency studies on this point (e.g. anvi'o pangenome "detection" or "MGE prevalence" studies) were not located in this session.
- PlasX, MOB-suite host assignment, and Hi-C methods were not reviewed in depth.

## 4. Geometric intuition: rays, cones and subspaces in coverage space

### Takeaway
Under the linear mixing model, each contig's S-dimensional coverage vector lies in the cone spanned by the k strain abundance vectors.
- **Core contigs** (carried by all strains) lie on one ray: the sum vector Σ_g a_g, scaled by copy number.
- **Strain-specific contigs** lie on the k extreme rays (the a_g themselves).
- **Contigs shared by a subset** lie on rays Σ_{g∈subset} a_g, i.e. on faces of the cone.

Finding strains is finding the extreme rays (vertices after normalisation onto the simplex). This is exactly separable NMF / archetypal analysis / conical-hull vertex finding. Accessory assignment is projecting each contig onto the face it lies on. Metagenomic tools implement this implicitly via NMF (DESMAN, StrainPanDA) or Bayesian factorisation (STRONG, StrainFacts). No tool I found frames it explicitly as convex-hull or archetype fitting.

### Cited Findings
- Separable NMF can be posed as "finding the extreme rays of the conical hull of a set of vectors". Archetypal analysis (convex NMF) constrains archetypes to convex combinations of data points; it is more interpretable but has higher reconstruction error. SNPA-type vertex-finding algorithms are sensitive to outliers. — [arXiv 1405.4275](https://arxiv.org/pdf/1405.4275); [Javadi & Montanari, arXiv 1705.02994](https://arxiv.org/pdf/1705.02994); [Near-Convex AA, arXiv 1910.00821](https://arxiv.org/pdf/1910.00821)
- Identifiability conditions are conical-geometry statements:
  - Separability means some data points lie exactly on the extreme rays. In this problem that means strain-specific contigs exist for every strain, or each strain is alone in some sample.
  - SSC only requires the data cone to contain the ice-cream cone. — [Fu et al. 2018, arXiv 1709.00614](https://arxiv.org/abs/1709.00614)
- DESMAN's NTF and StrainPanDA's sparse NMF are the operational implementations of this factorisation in metagenomics. — [PMC5607848](https://pmc.ncbi.nlm.nih.gov/articles/PMC5607848); [PMC10989911](https://pmc.ncbi.nlm.nih.gov/articles/PMC10989911)
- Co-abundance binners (CONCOCT, MSPminer, Nielsen canopies) use correlation or proportionality, i.e. the angle between coverage vectors, which is a direction/ray criterion. MSPminer explicitly uses proportionality, so genes on the same ray cluster together. — [PMC6499236](https://pmc.ncbi.nlm.nih.gov/articles/PMC6499236)

### Inferences
- **Why binners merge strains.** Correlation-based binners put the core ray and all strain rays in one bin whenever the angles between strain vectors are small. That is the case when the strains co-vary, which is the usual situation for near-identical strains. MSPminer-style proportionality would put strain-specific contigs on separate rays and might split them off as "accessory" clusters linked to the species core.
- **Proposed algorithm** (for the project; not from a published tool):
  1. Normalise each contig's coverage vector by the core (SCG) median per sample. This gives r_c ∈ ℝ^S, the carrier-fraction profile.
  2. Core contigs then cluster at r ≈ 1 (the all-ones vector). Strain-specific contigs lie on k extreme rays whose profiles π_g sum to 1.
  3. Find the extreme rays with separable NMF (SPA/SNPA, XRAY) or archetypal analysis on the non-core contigs.
  4. Assign each contig to a subset of strains by non-negative least squares or a Poisson/negative-binomial likelihood with binary weights. This is DESMAN eq. 7 with γ estimated geometrically instead of from SNVs.
  5. Check consistency with SNV-based π from DESMAN, STRONG or StrainFacts where SNVs exist.
- **Degenerate geometry.**
  - When the strain vectors are nearly collinear (constant ratios), the cone collapses to a thin wedge. The extreme rays are then poorly determined, and noise of about 1/√(reads) per contig can exceed the angular separation.
  - Expect needed S and depth to scale with 1/(angle between strain profiles)². This scaling is an inference from standard estimation theory and was not verified in a source.

### Gaps
- No published metagenomics method explicitly fits rays, cones or archetypes to contig coverage for strain-level accessory assignment. The search for "archetypal analysis / convex hull / separable NMF metagenomics strain" returned only general NMF literature.
- No source was found that gives angular-separation or noise-based sample-size formulas for this geometry in metagenomics.
