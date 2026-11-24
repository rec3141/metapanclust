# Point-cloud / LiDAR processing techniques for metagenomic binning and contig/bin embeddings

Scope note: about 20 searches and fetches. PMC pages blocked WebFetch (browser check), so the binny and GroopM full texts were pulled through the PubMed MCP full-text tool. Items marked **[unverified]** are standard references from background knowledge that were not re-checked in this session. Check them before citing.

## Q1. Has LiDAR / 3D point-cloud processing been applied to metagenomic binning or sequence embeddings?

### Takeaway
I found no published paper or preprint that explicitly brings LiDAR or 3D point-cloud pipelines (PCL-style filtering, RANSAC ground removal, Euclidean cluster extraction, PointNet) into metagenome binning. The closest precedents are of four kinds:
- binners that cluster contigs as points in a low-dimensional density space (binny: openTSNE + HDBSCAN; VizBin: BH-SNE scatter);
- GroopM's "transformed coverage profile space", with length-scaled point plots and interactive split/merge;
- point-cloud deep nets (PointNet descendants) applied to genomics data other than contigs (gene expression, single cells, smFISH spots);
- TDA (persistent homology, Mapper) applied to metagenome reads and microbiome profiles, but not to contig binning.

### Cited Findings
- A search for "point cloud metagenomic binning" turned up only conventional binners (MaxBin, MetaCluster 5.0, MyCC, BusyBee Web, a cellular-automata binner) and no point-cloud work — [MaxBin PMC4129434](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4129434/), [BusyBee Web PMC5570254](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5570254/), [MetaCluster 5.0](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3436824/)
- A search for "LiDAR algorithm applied to genomics embedding clustering" found no genomics method borrowing from LiDAR. The only LiDAR hit was a genetic-algorithm clustering paper for airborne laser scanning, with no genomics — [arXiv 1401.4848](https://arxiv.org/pdf/1401.4848)
- **binny** (Hickl, Queirós, Wilmes, May, Heintz-Buschart; Briefings in Bioinformatics 2022, doi 10.1093/bib/bbac431) is the closest thing to "point-cloud segmentation of contigs" — [PMC9677464](https://pmc.ncbi.nlm.nih.gov/articles/PMC9677464), [ORBilu record](https://orbilu.uni.lu/handle/10993/52754). Verified details from the full text (via PubMed MCP):
  - Features: k-mer frequencies for k = 2, 3, 4 plus average read coverage per sample, "both centered log-ratio transformed".
  - Embedding: FFT-accelerated openTSNE to 2D, with a multi-scale kernel (perplexity ranges 10–20 and 100–130) following Kobak & Berens. Early exaggeration and learning rate scale with the number of unbinned contigs.
  - Clustering: "HDBSCAN on the generated two-dimensional embedding, in iterations". min_samples is cycled over {1, 5, 10}, minimum cluster size = ln(n contigs), cluster_selection_epsilon cycles 0.25 → 0.0 in steps of 0.125, Manhattan distance.
  - Each cluster is scored for completeness and purity with lineage-specific CheckM marker sets. Good clusters are extracted as MAGs and the remaining unbinned contigs are re-embedded. Clusters below threshold are re-split with raw coverage added "as additional dimension(s)".
  - Minimum contig length starts high and "is decreased by 500 bp if less than half of the iterative clustering steps returned MAGs, until a minimum size of 500 bp". This is a coarse-to-fine schedule.
  - Single-contig genomes (≥40 markers, ≥90% pure, ≥92.5% complete) are pulled out first so they are not "indistinguishable from noise during clustering".
  - rRNA and CRISPR regions are masked before k-mer counting.
  - Rationale stated in the paper: "Clustering algorithms perform better in fewer dimensions, because distance information becomes increasingly imprecise at higher dimensions."
- **VizBin** (Laczny et al., Microbiome 2015, doi 10.1186/s40168-014-0066-1): parallel Barnes-Hut t-SNE on CLR-transformed genomic signatures gives a 2D scatter of contigs for "human-augmented" (manual) binning. It handles about 30,000 points in under 4 minutes on a 2011 laptop — [PMC4305225](https://pmc.ncbi.nlm.nih.gov/articles/PMC4305225), [VizBin site](https://claczny.github.io/VizBin)
- **GroopM** (Imelfort et al., PeerJ 2014, doi 10.7717/peerj.603) bins mainly on differential coverage across related metagenomes. Its workflow is parse → core → refine → recruit → extract. Its figures position contigs "according to their x and y coordinates in GroopM transformed coverage profile space" with "diameter of each circle … proportional to the length its respective contig". The software "allows users to view and merge bins, as well as split bins based on composition, coverage or contig length profiles". The authors note "GroopM's inability to reliably resolve bins at the strain level" (two S. epidermidis strains were merged into one bin) — [PMC4183954](https://pmc.ncbi.nlm.nih.gov/articles/PMC4183954/), [PeerJ](https://peerj.com/articles/603)
  - **[unverified]** My recollection is that GroopM projects multi-sample coverage onto a 3D space and uses Hough-transform-like or density "coring" in that space. That would make it the closest LiDAR-like precedent. The main-text sentences I extracted do not describe the transform (it is probably in the supplementary notes), so confirm before citing.
- **Point-cloud deep learning in genomics (not binning):**
  - GenePointNet / GPNet models gene-expression data as a point cloud with a PointNet backbone for tumour classification (bioRxiv 2024) — [bioRxiv 2024.06.02.597020](https://www.biorxiv.org/content/10.1101/2024.06.02.597020.full.pdf)
  - HiPoNet (NeurIPS 2025) is a multi-view simplicial-complex network for high-dimensional point clouds, with cells as points. It argues PointNet-style models do not suit single-cell data — [arXiv 2502.07746](https://pdf.arxiv.org/abs/2502.07746), [NeurIPS poster](https://neurips.cc/virtual/2025/loc/san-diego/poster/117742)
  - PointFISH uses PointNet-style shared MLPs and max-pooling on smFISH spot coordinates to classify RNA localisation patterns — [arXiv 2302.10923](https://arxiv.org/pdf/2302.10923)
- **TDA in metagenomics:**
  - Guzmán-Sáenz, Haiminen, Basu & Parida (BMC Genomics 2019, APBC 2019) used persistent-homology bar diagrams and a Čech-complex model for strain-level signal enrichment of multi-mapping reads. This is read-level, not contig binning — [BMC Genomics](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-5490-y), [PMC6456948](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6456948/)
  - tmap (Liao et al., Genome Biology 2019) applies Mapper to population-scale microbiome profiles for stratification. Samples are the points, not contigs — [Genome Biology](https://link.springer.com/doi/10.1186/s13059-019-1871-4), [GitHub](https://github.com/gpz-bioinfo/tmap)
  - Chang, VanInsberghe & Kelly (npj Biofilms Microbiomes 2020) used Mapper plus kNN density to define community states in gut and marine time series — [npj](https://www.nature.com/articles/s41522-020-00145-9)
  - A 2022 arXiv preprint does topology-based assembly-graph post-processing (cycle finding) to complement binning in HiFi assemblies. It does not use persistent homology — [arXiv 2210.00098](https://arxiv.org/pdf/2210.00098)
- **Current density / graph clustering in binners:**
  - COMEBin (Wang et al., Nat Commun 2024, 15:585) learns contrastive multi-view embeddings from augmented contig fragments and clusters them with Leiden — [PMC10794208](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10794208/), [GitHub](https://github.com/ziyewang/COMEBin)
  - A benchmark review found SemiBin2 and COMEBin best overall, with contrastive models on top — [PMC12636519](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12636519/)
  - VAMB encodes co-abundance and k-mers with a VAE, then clusters the latent space. The clustering is iterative medoid clustering per Nissen et al. 2021, Nat Biotechnol **[unverified this session]** — [PMC12636519](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12636519/)
  - MetaDecoder: my search returned nothing. Its method (reported as a DPGMM plus semi-supervised k-means/SVM) is **[unverified]**.

### Inferences
- The "glimpse of LiDAR applied to bins" most likely refers to the visual likeness between length-weighted contig scatters (GroopM coverage space, VizBin/binny t-SNE maps, anvi'o or BlobTools plots) and point-cloud renders. binny is already a working "segmentation of a 2D point cloud" pipeline: density clustering, an outlier/noise class, iterative extraction of segments, and coarse-to-fine scale. So the LiDAR vocabulary maps onto existing practice. The un-borrowed parts are robust geometric model fitting (RANSAC), local shape descriptors (linearity/planarity), registration (ICP), and supervoxel over-segmentation.
- The gap is genuine. Saying "no direct LiDAR-in-metagenomics literature exists" is defensible on these searches. There may still be blog posts, talks or unindexed preprints.

### Gaps
- Could not find the source of the user's "glimpse". Searches on LiDAR + bins / metagenomics returned nothing.
- GroopM's exact coverage transform (3D projection, Hough?) is not confirmed. Check the PeerJ supplementary notes.
- SemiBin2's clustering (I believe Leiden/Infomap on a kNN graph of the embedding) and MetaDecoder's method were not verified here.
- I did not search for anvi'o/BlobTools 3D views or Hi-C contact-based binners as "multi-channel" analogues.

## Q2. Standard LiDAR / point-cloud pipeline steps and their analogues for binning (contigs as length-weighted points in coverage × k-mer space or an embedding)

### Takeaway
Nearly every stage of the PCL pipeline has a usable analogue. The most promising transfers for strain-aware binning are:
- RANSAC / robust line fitting in log multi-sample coverage space, to peel off a species' dominant "ray" and examine the residuals;
- local eigenvalue shape features (linearity/planarity/scattering) computed on kNN neighbourhoods, to flag contigs sitting on elongated (strain-mixture) structures versus compact blobs;
- statistical outlier removal and length-weighted voxel downsampling before clustering;
- region growing / label propagation along assembly-graph adjacency, the analogue of scan-line connectivity (already done by GraphBin, METAMVGL, BinSPreader, GraphMB);
- Procrustes/ICP alignment of embeddings across runs or samples.

### Cited Findings
- PCL (Rusu & Cousins, "3D is here: Point Cloud Library (PCL)", IEEE ICRA 2011, doi 10.1109/ICRA.2011.5980567) groups its algorithms into "filtering, feature estimation, surface reconstruction, registration, model fitting and segmentation". Filtering is used, e.g., "to filter outliers from noisy data" — [PCL ICRA 2011 PDF](https://www.cs.utexas.edu/~jsinapov/teaching/cs378/readings/W10/pcl_icra2011.pdf), [pointclouds.org](https://pointclouds.org/about/)
- Weinmann, Jutzi, Hinz & Mallet (ISPRS J. Photogramm. Remote Sens. 105:286–304, 2015) describe a four-stage pipeline: neighbourhood selection → feature extraction → feature selection → classification. They evaluate 7 neighbourhood definitions, 21 geometric features and 10 classifiers, and find that choosing an *optimal per-point neighbourhood* (eigenentropy-minimising scale) matters most — [TUM portal](https://portal.fis.tum.de/en/publications/semantic-point-cloud-interpretation-based-on-optimal-neighborhood/), [DGPF 2017 follow-up](https://dgpf.de/src/tagung/jt2017/proceedings/proceedings/papers/18_DGPF2017_Weinmann_et_al.pdf)
- The eigenvalue features (λ1 ≥ λ2 ≥ λ3 of the local covariance) are linearity, planarity, sphericity/scattering, omnivariance and eigenentropy. Exact formula variants differ between papers — [Thomas et al. 2018, arXiv 1808.00495](https://arxiv.org/pdf/1808.00495)
- Differential coverage premise: contigs whose coverage co-varies across related samples come from the same population. CoMet notes that clustering is "prone to outliers" and uses robust outlier filtering for precision, at the cost of fewer contigs binned — [GroopM PeerJ](https://peerj.com/articles/603), [CoMet PMC5751405](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5751405/)
- Assembly-graph adjacency already serves as a "connectivity" channel:
  - GraphBin: label propagation on the assembly graph to correct and extend an existing binning.
  - METAMVGL: assembly graph + paired-end graph with learned weights. It "could strikingly rescue the short contigs and correct the binning errors from dead ends".
  - BinSPreader: propagation over graph edges, allowing multi-binning of shared edges.
  - GraphMB: GNNs fed the assembly graph during binning itself.
  - Sources — [METAMVGL PMC8296540](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8296540/), [GraphBin GitHub](https://github.com/Vini2/GraphBin), [GraphMB bioRxiv 2022.02.25.481923](https://vbn.aau.dk/ws/files/469239833/2022.02.25.481923v1.full_1_.pdf)
- binny already applies these point-cloud stages, though not under those names (see Q1):
  - a multi-scale neighbourhood (multi-perplexity kernel);
  - density segmentation with a noise class (HDBSCAN);
  - iterative "extract segment, re-process remainder", similar to sequential RANSAC / iterative ground removal;
  - a coarse-to-fine length threshold, similar to voxel-size schedules;
  - masking of rRNA/CRISPR regions, similar to removing spurious returns — [PMC9677464](https://pmc.ncbi.nlm.nih.gov/articles/PMC9677464)

**Mapping table.** The analogues are my synthesis. Algorithm references marked [unverified] are classic citations from background knowledge.

| LiDAR / point-cloud step | Classic reference | Binning analogue | Notes / existing precedent |
|---|---|---|---|
| Voxel-grid downsampling (centroid per voxel) | PCL VoxelGrid [unverified] | Collapse near-duplicate contigs, or tile a long contig into fixed-length chunks so each point carries comparable mass. Alternatively build a length-weighted grid in an embedding for density estimation. | Length weighting: one 100 kb contig should count for more than one 1 kb contig. Voxelising (chunking long contigs, as COMEBin's fragment "views" do) turns weight into point density. |
| Statistical outlier removal (mean kNN distance > μ + kσ) | PCL StatisticalOutlierRemoval [unverified] | Drop or defer contigs whose kNN distance in coverage+TNF space is anomalous: short noisy contigs, plasmids, phage, mobile elements, chimeras. | Same as CoMet's robust filtering and HDBSCAN's noise label. Run it before clustering, then recruit the outliers back afterwards (GroopM's "recruit" stage). |
| Radius outlier removal | PCL [unverified] | Minimum-neighbour rule at a fixed radius in CLR-coverage space | Simple pre-filter |
| Ground-plane removal by RANSAC | Fischler & Bolles 1981, CACM [unverified] | Fit and remove the dominant structure before segmenting the rest, e.g. the main ray of a dominant species in log-coverage space, or the "low-coverage floor" of rare contigs. | Sequential RANSAC = iteratively extract the best-supported model, remove its inliers, repeat. This matches binny's extract-and-re-embed loop. |
| RANSAC line/plane fitting (model fitting) | Fischler & Bolles 1981 [unverified] | **Strain rays.** In multi-sample coverage space, contigs of one genome lie on a ray through the origin (coverage ∝ abundance vector, with scale set by copy number and mapping biases). In log space the ray becomes a line with slope 1 offset by a per-contig constant. Contigs shared by two co-occurring strains lie near the *sum* of two rays, a 2D plane/cone spanned by the strain abundance vectors. Strain-specific accessory contigs lie on the individual rays. RANSAC with a ray model (2 parameters per sample, after normalising out contig scale) finds the core ray. A 2-ray plane model finds a strain pair. | I found no published binner doing this explicitly (gap). Non-negative matrix factorisation of coverage (e.g. DESMAN/STRONG-style approaches for strains) is the statistical cousin [unverified]. |
| Hough transform | Hough 1962; Duda & Hart 1972 [unverified] | Vote over ray directions (points on the unit simplex of normalised coverage). Peaks are species/strain abundance vectors. Length-weight the votes. | Normalised coverage direction = CLR/proportion of coverage across samples, so a Hough accumulator on the simplex amounts to a length-weighted density of coverage directions. GroopM may have done something like this [unverified]. |
| Euclidean cluster extraction / region growing | PCL EuclideanClusterExtraction; region growing with normal/curvature smoothness [unverified] | Connected components on a radius graph in feature space. Region growing seeded at high-marker, long contigs and grown while coverage-ratio smoothness holds. Growth along **assembly-graph edges**. | GraphBin / METAMVGL / BinSPreader are label propagation along assembly adjacency. |
| DBSCAN / HDBSCAN / OPTICS segmentation | Ester et al. 1996 KDD; Campello et al. 2013; Ankerst et al. 1999 [unverified] | Density segmentation of embedding | Already used by binny (HDBSCAN on openTSNE). Length-weighted HDBSCAN can be done by chunk-replicating contigs or with sample weights. |
| Normals / curvature / eigenvalue shape features (linearity, planarity, scattering, eigenentropy) | Weinmann et al. 2015; Demantké et al. 2011 [unverified] | For each contig, take the PCA of its kNN neighbourhood in log-coverage (+TNF) space. High **linearity** means a 1D ray structure: a single genome with copy-number or GC bias spread, or a strain gradient. **Planarity** means a mixture of two abundance vectors (shared core of two strains). **Scattering** means a noise/chimera zone. | Novel for binning as far as I found. Could flag bins that are strain composites before CheckM contamination shows it (shared-core contamination is invisible to single-copy marker checks). |
| Multi-scale / optimal neighbourhood (eigenentropy-minimising k) | Weinmann et al. 2015 | Choose kNN size per contig, or use a multi-scale t-SNE kernel | binny's multi-perplexity kernel (Kobak & Berens) is the embedding-side version |
| Supervoxels (VCCS) / over-segmentation then merge | Papon et al. 2013 CVPR [unverified] | Over-segment into small, high-purity "micro-bins" (e.g. HDBSCAN leaf clusters or graph components), then merge micro-bins using marker-gene completeness/redundancy and coverage-ray agreement. | Matches binny's cycling of cluster_selection_epsilon from 0.25 to 0 (micro-cluster merging), and DAS Tool / MetaWRAP refinement. |
| Multi-return / intensity channels | LiDAR multi-echo, intensity | Extra per-contig channels beyond geometry: Hi-C contacts, methylation motifs (nanopore 4mC/5mC/6mA), taxonomy votes, marker genes, read-level SNV density. | Methylation-motif profiles are a strong "intensity" channel for separating co-abundant genomes [not searched this session]. |
| Scan-line ordering / sensor topology | LiDAR ring/scan-line adjacency | Assembly graph (GFA) adjacency, paired-end/long-read links, contig-end overlaps | GraphMB, METAMVGL, GraphBin, BinSPreader |
| ICP / Procrustes registration | Besl & McKay 1992 PAMI [unverified] | Align embeddings from (a) different t-SNE/UMAP seeds, (b) different samples or co-assemblies, (c) successive runs, to track contigs or bins that move between clusters. Contig IDs give known correspondences, so plain Procrustes is enough within one assembly. ICP-style nearest-neighbour correspondence is needed across different assemblies (e.g. by bin centroid). | scStability uses Procrustes for exactly this between embedding seeds (Q4). |
| Ground-truth semantic segmentation nets (PointNet/PointNet++, KPConv) | Qi et al. CVPR 2017 / NeurIPS 2017 [unverified]; Thomas et al. 2019 | Learn per-contig labels from neighbourhood point sets (CAMI-trained) | GenePointNet, HiPoNet, PointFISH show the architecture transfers to genomics point sets (Q1). Not applied to binning. |

### Inferences
- **Length weighting** is the main technical difference from LiDAR, where points are roughly equal-mass. Length enters naturally as a point weight in density estimates (weighted KDE / weighted HDBSCAN core distances), as a vote weight in RANSAC/Hough, and as noise variance. Coverage and TNF variance of a contig shrink roughly as 1/length, so short contigs should get looser inlier thresholds, analogous to range-dependent noise in LiDAR. A practical trick is to chunk long contigs into fixed windows (voxelisation along the sequence). Each window becomes a point, the window cloud of one contig gives its own local spread, and contigs with internally discordant windows are chimeras or misassemblies.
- **Strain rays** are the most "LiDAR-native" idea. With more than 3 samples, a genome is a ray in coverage space and two co-occurring strains span a 2D cone. RANSAC with ray and cone models plus local planarity features could identify shared-core vs strain-specific contigs. A dominant-ray removal step (fit the species ray, subtract or project out, recluster residuals) is the "ground removal" analogue. This needs enough samples with varying strain ratios; with one sample all genomes collapse onto a 1D coverage axis.
- Of the steps above, only the shape features and RANSAC ray/cone fitting appear to have no precedent in binners. The others already exist under binning names.

### Gaps
- I found no published use of RANSAC, Hough transforms or eigenvalue shape features on metagenomic coverage profiles. GroopM's coring method is the one to check.
- Exact PCL defaults for SOR/voxel filters and the original RANSAC/ICP/DBSCAN/PointNet citations were not re-verified this session.
- I did not verify strain-deconvolution tools that model coverage as mixtures (DESMAN, STRONG, MetaCarvel).

## Q3. RANSAC / robust line-plane fitting, Hough transforms and manifold detection in multi-sample coverage space

### Takeaway
I found no paper applying RANSAC or Hough transforms to contig coverage for strain or core detection. Differential-coverage binners (GroopM, CoMet, CONCOCT, MetaBAT2) treat co-varying coverage as the signal and use robust outlier filtering. Fitting a geometric ray/cone model of strain mixtures seems to be an open niche.

### Cited Findings
- GroopM bins "primarily [on] differential coverage profiles" and cannot reliably separate strains: two S. epidermidis strains were combined into one bin — [PMC4183954](https://pmc.ncbi.nlm.nih.gov/articles/PMC4183954/)
- CoMet adopts "robust outlier filtering strategies … to improve precision in binning", which reduces the number of contigs binned, and reports better precision in a multi-strain sample — [PMC5751405](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5751405/)
- A search for coverage-space robust regression and strain ratios found no relevant method (only general differential-coverage papers) — [GroopM PeerJ](https://peerj.com/articles/603)
- binny re-splits impure clusters by adding raw coverage depth as extra dimensions. This is a heuristic way of separating along coverage magnitude (the "position along the ray") as well as direction — [PMC9677464](https://pmc.ncbi.nlm.nih.gov/articles/PMC9677464)

### Inferences
- Proposed recipe (untested synthesis):
  1. Compute per-contig coverage vectors c over S samples. Normalise direction u = c/‖c‖ (or CLR) and keep the magnitude separately.
  2. Run a length-weighted Hough/density accumulator over u to find candidate abundance directions (genomes).
  3. For each candidate, run RANSAC with a ray model in log space, with inlier tolerance scaled by 1/√(length·depth).
  4. For residual contigs, fit 2-ray cones (planes through the origin). Inliers to a cone but not to either ray are candidate shared-core contigs of a strain pair, and ray inliers are strain-specific.
  5. Score with marker genes, as binny does.
- In t-SNE/UMAP space the rays become curved filaments or elongated clusters. Local linearity features (the PCA eigenvalue ratio of kNN neighbourhoods) or principal-curve fitting can detect them. The original coverage space is better for line fitting because the embedding distorts the geometry.

### Gaps
- No empirical evidence found that ray/cone models outperform NMF-style strain deconvolution. Coverage bias (GC, mapping to shared regions, ambiguous multi-mapping at shared cores) may break the linear model, and this was not investigated.

## Q4. Stabilising t-SNE/UMAP clustering across random initialisations (consensus) — detecting fragments that do not co-cluster with their genome

### Takeaway
Established tools and papers aggregate many embeddings or clusterings across seeds: Procrustes-aligned stability scores, co-association / averaged pairwise distances, and geometric-median consensus embeddings. Applied to binning, the per-contig co-assignment frequency across seeds is a direct "this fragment wanders" score. I found no binner that does this explicitly. binny's cycling of HDBSCAN parameters and iterative re-embedding is the nearest analogue.

### Cited Findings
- **Median Consensus Embedding (MCE)** takes the geometric median of multiple embeddings to offset random-initialisation variability. It "converges rapidly and significantly reduces instability". The paper notes the earlier consensus-embedding framework of Viswanath & Madabhushi lacks quantitative evaluation of seed dependence — [arXiv 2503.08103](https://arxiv.org/pdf/2503.08103)
- **scStability** (CRAN R package) generates many UMAP/t-SNE embeddings from different seeds. It scores embedding stability with pairwise Procrustes alignment plus Kendall's τ, and cluster stability with NMI between Louvain/Leiden assignments across runs — [CRAN scStability PDF](https://stat.ethz.ch/CRAN/web/packages/scStability/scStability.pdf)
- **tSNEd** (informal GitHub project, not peer reviewed): 1000 t-SNE runs found pairwise sample distances "mostly consistent across runs". Averaging the pairwise distances and then running hierarchical clustering gives a consensus — [tSNEd README](https://rdrr.io/github/zhezhangsh/tSNEd/f/README.md)
- One comparative arXiv study on MNIST over 100 seeds found K-means on t-SNE most variable (ARI SD 0.045 vs 0.028 for UMAP) — [arXiv 2503.23215](https://arxiv.org/pdf/2503.23215)
- UMAP's visual output can change substantially with the random seed — [arXiv 2103.14608](https://arxiv.org/pdf/2103.14608)
- binny's own run-to-run handling: HDBSCAN is repeated over min_samples ∈ {1, 5, 10} and a cycled cluster_selection_epsilon, and the embedding is recomputed on the remaining contigs each iteration — [PMC9677464](https://pmc.ncbi.nlm.nih.gov/articles/PMC9677464)

### Inferences
- Concrete binning use:
  1. Run N embeddings (seeds and perplexities) plus HDBSCAN.
  2. Build a length-weighted co-association matrix restricted to kNN pairs for scale.
  3. For each contig, compute the fraction of runs in which it lands in the same cluster as the majority of its current bin (weighting bin-mates by length).
  4. Contigs with low agreement are candidates for mis-binning: mobile elements, strain-variable regions, chimeras. Bins whose members split consistently along the same partition are candidates for strain splits.
- This is the "ICP registration" analogue: contig identities provide correspondences, so Procrustes aligns runs exactly without ICP's nearest-neighbour matching.
- Consensus over seeds averages out embedding noise but not systematic biases (GC, coverage), so it cannot replace checks in the original feature space.

### Gaps
- I found no metagenomics-specific paper evaluating seed consensus for binning. The Viswanath & Madabhushi consensus-embedding paper and the classic consensus-clustering work (Monti et al. 2003, Mach Learn) were not fetched — **[unverified]**.
