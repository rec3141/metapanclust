# Related work: what to reuse, what to beat (2026-11-18)

**Bottom line.** Nobody computes what this project proposes: a calibrated probability of each conserved
gene cluster's state in a genome or population bin, given its clade and the rest of the genome, with
absence corrected for completeness, contig ends and depth, and then used to score bins. The pieces exist
separately:

- **Cluster discovery:** Spacedust, PPanGGOLiN panModule, Panorama, anvi'o pangenome graphs.
- **Within-species prevalence:** PPanGGOLiN partitions, precomputed for more than 4,600 GTDB species in
  PanGBank (2026).
- **Completeness-aware presence:** CheckM collocated sets, mOTUpan, and the 2026 phylogenetic occupancy
  model.
- **Gene-gene dependence:** Coinfinder, Goldfinder, and the random-forest pangenome predictor of Beavan
  et al.
- **Contextual prediction:** gLM, Bacformer, MetaPathPredict.
- **Coverage-based gene presence in metagenomes:** MIDAS2, PanPhlAn, inStrain.

None of them joins prevalence across the tree, detection and context into one likelihood, and none uses
gene content for merge or move decisions beyond single-copy markers (CheckM, DAS Tool). I found no prior
use of surprisal (-log P) as a genome-content or binning score (searches below).

**Two threats to plan around:**

- **PanGBank** already gives per-species prevalence, partitions and modules for 393k GTDB genomes. Do
  not rebuild species pangenomes; reuse them.
- **The phylogenetic occupancy model** (Mattick et al. 2026) already combines phylogeny, presence and
  completeness for single genes. Level 1 has to beat it on held-out data, or adopt it.

Verification legend:

- **[v]** Checked against the paper, the official docs or the code today.
- **[a]** Abstract only.
- **[u]** Unverified: from memory or a secondary source.

---

## 1. Conserved gene cluster discovery

### Spacedust (Zhang, Mirdita & Söding, Nat Methods 22:2065-2073, 2025)
[doi:10.1038/s41592-025-02816-x](https://doi.org/10.1038/s41592-025-02816-x)

**What it computes** [v, PMC12510874]:

- **Unit of analysis.** Each query and target genome pair separately (sets Q and T, possibly the same).
- **Homology search.** Proteins are matched with Foldseek (MMseqs2 for proteins without a structure),
  e-value 1e-3 and coverage 0.8. The best hit per query protein and target genome is kept if p ≤ 1e-7.
- **Clustering P value.** The probability, under a null where any query protein hits any target protein
  with probability q0 = 1e-3, of at least k hits within a window of span m genes.
- **Order-conservation P value.** The probability of at least n neighbouring pairs that keep adjacency and
  relative strand.
- **Score and clustering.** The two P values are multiplied as if independent. Clusters are grown by
  greedy agglomerative merging on a score built from their negative logs.
- **Output.** Pairwise cluster matches. These are grouped afterwards into "nonredundant clusters", but not
  into named families with presence/absence matrices.

**Scale** [v]:

- 1,308 genomes (4.19 M genes).
- Search and filtering took 72 h on 2 × 2 × 64-core EPYC (about 150 ms per genome pair). Clustering took
  51 min.
- Output: 72,483 nonredundant clusters covering 58% of genes. The abstract says 72,843; the paper is
  internally inconsistent.
- Recovered 95% of PADLOC defence clusters (5,255/5,520).
- Mean F1 0.61 on 207 hand-annotated BGCs.
- Runtime is quadratic in genome number.

**Limitations the authors state** [v]: partial conservation predicts functional association "with only
moderate precision"; functional modules split across loci cannot be found; and genes hitchhike with
clusters without being part of them.

**What we can reuse:**

- The binary and its pairwise hit tables (already running: `scripts/run_spacedust.sh`).
- The p_clu formula as a per-pair co-location statistic.
- The PADLOC recall test as a catalogue benchmark.

**Where we differ:**

- **The null is wrong for close relatives.** It is a uniform-hit null with no phylogeny. Between two
  genomes of one genus nearly every protein hits, so q0 = 1e-3 makes every syntenic block "significant".
- **Conservation is not separated from inheritance.** Spacedust measures co-location conservation, not
  conservation in excess of shared ancestry.
- **No cross-genome families.** It produces no families, so every merge into a family is our problem. The
  first test shows single linkage collapsing into one giant family (`docs/first_tests.md`).

We need families defined by gene-family content (Panorama's GHGF/GFRR approach, below) and a co-location
propensity per family that is corrected for divergence time. The Spacedust authors list a
divergence-aware score as future work.

### PPanGGOLiN, panRGP, panModule (Gautreau et al. PLoS Comput Biol 2020; Bazin et al. Bioinformatics 2020; panModule preprint 2021)

- PPanGGOLiN: [doi:10.1371/journal.pcbi.1007732](https://doi.org/10.1371/journal.pcbi.1007732)
- panRGP: [doi:10.1093/bioinformatics/btaa792](https://doi.org/10.1093/bioinformatics/btaa792)
- panModule preprint: [doi:10.1101/2021.12.06.471380](https://doi.org/10.1101/2021.12.06.471380)

**What it computes:**

- **Partitions** [a]. A graph with gene-family nodes and genomic-neighbour edges. The nodes are partitioned
  into persistent, shell(s) and cloud by EM on a multivariate Bernoulli mixture coupled to a Markov
  random field over the graph. The paper tested 439 isolate species and MAGs of 78 species.
- **RGPs** [a]. panRGP calls runs of shell and cloud families as regions of genome plasticity, and groups
  them into insertion spots by their flanking persistent families.
- **Modules** [v, `ppanggolin/mod/module.py`, master today]:
  - Build a graph with a transitive closure of t = 4 genes.
  - Drop families present in fewer than 2 genomes.
  - Drop edges whose Jaccard (co-occurrence within the window) is below 0.85.
  - Keep connected components with at least 3 families.
  - **Discard any component containing a non-multigenic persistent family.** Modules are therefore
    accessory-only by construction.
- **Projection** [v, docs/user/projection.md]. Projects new genomes onto an existing pangenome, assigning
  families, partitions, RGPs, spots and modules. RGPs on contig borders are left unassigned to spots.

**PanGBank (Mainguy et al., NAR 2026)** [a, [doi:10.1093/nar/gkag948](https://doi.org/10.1093/nar/gkag948)]:

- Precomputed PPanGGOLiN pangenomes for more than 4,600 GTDB species, covering more than 393,000 genomes.
- Two collections: GTDB_all (with MAGs and SAGs) and GTDB_refseq.
- Includes panRGP and panModule results, with REST API and CLI access.

**What we can reuse:**

- Species-level presence/absence matrices and partitions from PanGBank, which saves stage 2's per-species
  work.
- `ppanggolin projection` to place our bins into those pangenomes.
- Module definitions as one cluster source for accessory content.

**Where we differ:**

1. **Modules exclude persistent families.** The most informative absences (a core operon missing from a
   bin) are therefore invisible to panModule. Our first prototype scored only modules, so it could only
   ever see accessory content. Core operons need their own units: Spacedust clusters, or co-located runs
   of persistent families.
2. **Families do not carry across species.** Shrinking up the tree needs cross-species cluster identity.
3. **Partition is a 3-class label, not a probability.** It is not shrunk across species and not corrected
   per genome for completeness. PPanGGOLiN does handle MAG incompleteness inside the mixture model
   [u: details of how MAG missingness is treated were not re-read].

### PANORAMA (Arnoux et al., PLoS Comput Biol 22:e1013856, 2026)
[doi:10.1371/journal.pcbi.1013856](https://doi.org/10.1371/journal.pcbi.1013856)

**What it computes** [v, PMC13379101]:

- **System prediction.** Rule-based "systems" on top of PPanGGOLiN graphs, using MacSyFinder-style
  grammar plus a Functional Unit level:
  - Family roles are mandatory, accessory, neutral or forbidden, with quorum rules.
  - Families are annotated with HMMs. 667 models and 6,272 HMMs were converted from DefenseFinder,
    PADLOC and CasFinder.
  - Components may be at most a `transitivity` distance apart.
  - Graph edges with Jaccard below 0.8 (computed over genomes carrying the families) are removed.
- **Per-genome calls.** Each genome is called *strict*, *extended* or *split*. Split means separated by
  other genes or contigs.
- **Cross-species comparison:**
  - Families from different pangenomes are grouped into GHGFs (MMseqs2, 50% identity, 80% coverage).
  - Spots are linked when the flanking-family overlap is at least 0.6 (minGFRR), then grouped by Louvain.

**Benchmarks and scale:**

- Family-level F1 of 98% against PADLOC and 96.6% against DefenseFinder on 941 *P. aeruginosa* genomes.
- More than 6,000 Enterobacteriaceae genomes across four species.

**What we can reuse:**

- The GHGF step is the simplest working recipe for cross-species cluster identity: cluster family
  representatives at 50/80, then define a cluster by its GHGF set.
- GFRR (shared fraction over the smaller or larger set) is a ready content-similarity measure for
  merging pairwise Spacedust loci into families. That is exactly the "next step" in `first_tests.md`.
- The strict/extended/split labels already encode a three-state observation model, matching our
  intact / partial / absent.

**Where we differ:** systems come from curated rules. We want rule-free clusters, and prevalence plus
surprise rather than detection calls.

### Synteny-aware pangenome graphs in anvi'o (Henoch et al., bioRxiv 2026)
[doi:10.64898/2026.07.03.736256](https://doi.org/10.64898/2026.07.03.736256) (posted 2026-07-03);
tutorial at <https://merenlab.org/tutorials/undatipelagibacter-pangenome-graph/>.

**What it computes** [a]:

- Network-pruning and layout algorithms for gene-cluster graphs that keep synteny.
- Applied to 29 *Undatipelagibacter* genomes (SAR11 Ia.3.VI).
- Finds that variability is "a structured continuum, not a few hypervariable islands". Some variable
  regions have a "propensity to vary [that] is conserved across genera".
- Construction details were not read; the full text was not accessible.

**What we can reuse:**

- Anvi'o graph exports (`anvi-pan-genome-graph`) for SAR11-like marine clades.
- The "conserved propensity to vary" finding supports sharing each cluster's dispersion parameter (κ_k in
  `model_spec.md`) across genera instead of fitting it per genus.

**Where we differ:** it is a within-genus visual and analytical tool with no probabilistic expectation.
The tools are on anvi'o's development branch only [v, tutorial page].

### Operon and neighbourhood databases (evaluation and priors, not methods to beat)

- **DOOR / DOOR 2.0** (Mao et al. NAR 2009, 2014; Cao et al. Brief Bioinform 2019):
  - 6,975,454 predicted operons in 2,072 complete genomes, plus RNA-seq-derived TUs for 24 genomes [a].
  - [doi:10.1093/bib/bbx088](https://doi.org/10.1093/bib/bbx088),
    [doi:10.1093/nar/gkt1048](https://doi.org/10.1093/nar/gkt1048).
  - Whether the server is still online was not checked [u].
- **ProOpDB** (Taboada et al. NAR 2012):
  - Operon predictions for more than 1,200 genomes [a].
  - [doi:10.1093/nar/gkr1020](https://doi.org/10.1093/nar/gkr1020).
- **Uber-operons** (Che et al. NAR 2006):
  - Groups of operons linked through gene-set conservation across 90 reference genomes. Found 158
    uber-operons covering 1,830 genes in *E. coli* K-12 [a].
  - This is the earliest version of "clusters of clusters" and is the right granularity for checking our
    family-of-loci step.
  - [doi:10.1093/nar/gkl294](https://doi.org/10.1093/nar/gkl294).
- **STRING 12** (Szklarczyk et al. NAR 2023):
  - "Conserved genomic context" is one evidence channel. Scores are transferred to other organisms through
    hierarchical orthology, and the network can be built for user-submitted genomes [a].
  - The exact neighbourhood-score formula was not re-read [u].
  - [doi:10.1093/nar/gkac1000](https://doi.org/10.1093/nar/gkac1000).

**Use:**

- **DOOR and ProOpDB:** silver-standard operon pairs across about 2,000 genomes, beyond the EcoCyc pairs
  used now. They are predictions, so use them to test recall only, never as truth for precision.
- **STRING neighbourhood scores:** a pairwise co-location prior to compare against our learned
  co-location propensity.

Review background: Rogozin et al. 2004, gene neighbourhood methods,
[doi:10.1093/bib/5.2.131](https://doi.org/10.1093/bib/5.2.131) [a].

---

## 2. Co-occurrence beyond inheritance (level-2 context)

### Coinfinder (Whelan, Rusilowicz & McInerney, Microb Genom 2020)
[doi:10.1099/mgen.0.000338](https://doi.org/10.1099/mgen.0.000338)

**What it computes:**

- Tests gene-family pairs for association or dissociation beyond chance [a].
- Uses the tree to measure each gene's lineage dependence (Fritz & Purvis D) and excludes genes with low
  D [v via the Goldfinder text].
- The pair test is reportedly binomial on terminal co-occurrence counts with a phylogeny-independent null
  [u].

### Goldfinder (Gavriilidou, Paulitz, Resl, Ziemert, Kupczok & Baumdicker, bioRxiv 2024; extension in Phil Trans B 2025)

- Preprint: [doi:10.1101/2024.04.29.591652](https://doi.org/10.1101/2024.04.29.591652)
- Phil Trans B extension (Kupczok et al.): [doi:10.1098/rstb.2024.0070](https://doi.org/10.1098/rstb.2024.0070)
- Code: <https://github.com/fbaumdicker/goldfinder>

**What it computes** [v, preprint full text]:

- **Ancestral states.** Fitch parsimony on the rooted tree. Each change is a gain or loss.
- **Default statistic.** The *simultaneous score* A_S = B_S - B_D: branches where both genes change in the
  same direction, minus branches where they change in opposite directions. Terminal and Coinfinder-style
  scores are also offered.
- **Null distribution:**
  - Each gene is simulated independently on the tree.
  - The root state is drawn from the empirical root frequency.
  - The number of changes comes from the empirical homoplasy distribution, placed on branches in
    proportion to branch length (constant rate).
  - The simulated genes go through the same reconstruction, which gives null scores and empirical P
    values with FDR.
- **Networks.** MCL clustering of the significant pairs.
- **Result.** Applied to the Coinfinder *S. pneumoniae* data, it finds 18,443 significant pairs versus
  Coinfinder's 76,416.
- Runtime and genome counts are not given in the preprint.
- No journal version of the original method was found. The 2025 Phil Trans B paper extends it to
  phage-host matrices [a].

### Beavan, Domingo-Sananes & McInerney (PNAS 121:e2304934120, 2024)
[doi:10.1073/pnas.2304934120](https://doi.org/10.1073/pnas.2304934120)

**What it computes** [a]:

- A random forest per accessory gene predicts the gene's presence from the presence or absence of the
  other accessory genes in *E. coli*.
- Many genes are "highly predictable from other genes alone".
- Secondary sources say about 2,500 complete genomes [u].

This is the **closest precedent for our level-2 masked model**, at single-species scale.

### What we can reuse and where we differ

**Reuse:**

- Goldfinder's simulation null is the cheapest correct way to filter level-2 features. Only let cluster
  k' predict k when their association survives a phylogeny-conditioned null; otherwise shared clade
  membership leaks in twice, once through level 1 and once through the context.
- Beavan et al. is the design and baseline for level 2.

**Where we differ:**

- **Units.** Our units are clusters, not genes. Most "linkage" co-occurrence (genes adjacent on the
  chromosome) is absorbed into the cluster itself, so what remains between clusters is mostly
  co-selection or co-transfer. That is the signal these tools hunt for.
- **Prediction, not testing.** We need predictions for MAGs from new clades, not significance of pairs
  within one species.
- **Offset.** Our model must be conditioned on the level-1 clade expectation (logit offset), which
  neither tool provides.

---

## 3. Genome language models (an alternative level 2)

### gLM (Hwang, Cornman, Kellogg, Ovchinnikov & Girguis, Nat Commun 15:2880, 2024)
[doi:10.1038/s41467-024-46947-9](https://doi.org/10.1038/s41467-024-46947-9)

**What it computes** [v, PMC10991518]:

- **Input.** Subcontigs of 15-30 genes. Each gene is a 1,280-d ESM2 embedding plus an orientation
  feature.
- **Training objective.** 15% of genes are masked. The model predicts four candidate embeddings per masked
  gene with likelihoods.
- **Training data.** 7,324,684 MGnify subcontigs.
- **Model.** RoBERTa-style, 19 layers, about 1B parameters.
- **Operons.** One attention head correlates with *E. coli* operons (ρ = 0.44). A logistic probe on all
  heads gives operon mAP 0.775.
- **Enzyme function.** EC linear probe accuracy is 51.6%, a 4.6-point gain over context-free embeddings.

### Bacformer (Wiatrak et al., bioRxiv 2025, v3 2026-08-19)
[doi:10.1101/2025.07.20.665723](https://doi.org/10.1101/2025.07.20.665723)

**What it computes** [v for v1 full text; v3 abstract only]:

- **Input.** Whole genomes as ordered lists of proteins: mean-pooled ESM-2 35M embeddings plus contig
  embeddings. Contig order is shuffled per epoch for MAGs.
- **Training objective.** Masked (15%) or autoregressive prediction over 50,000 protein-family "tokens"
  (k-means centroids of ESM-2 embeddings).
- **Training data.** 1,314,344 MAG-derived genomes (MGnify v2, SPIRE; 2.7 B proteins), then 48,932
  complete genomes.
- **Model.** 12 layers, hidden size 480.
- **Tasks.** Operons (zero-shot), STRING PPI, essentiality (median AUROC 0.81), 139 phenotypes and
  "genome completion".
- **Robustness to completeness.** Reported as "relatively robust until completeness drops below 75%".
- **Not addressed in v1:** MAG contamination and presence/absence probabilities.
- v3 adds authors and experimental operon validation [a]. Whether v3 changed the vocabulary or size was
  not checked.

### What we can reuse and where we differ

**Reuse:**

- Both are baselines for level 2. Bacformer's masked-family objective is formally the same as our masked
  per-cluster prediction, with learned rather than observed units.
- Pretrained weights are public on Hugging Face (`macwiatrak/*`) [v, search listing]. Bacformer
  embeddings of a bin could also serve as a context feature.

**Where we differ:**

- **Calibration.** Neither outputs calibrated P(cluster present | clade, rest of genome), so surprise in
  bits is undefined without recalibration.
- **No explicit phylogenetic correction or detection model.** Trained on MAGs, they learn the MAG
  artefacts (missing rRNA and ribosomal-protein contigs, below) as "normal".
- **Opaque and costly.** They are hard to interpret, and their cost scales with model size. Our
  per-cluster logistic models are auditable.

---

## 4. Genome quality from gene content (the narrow existing form of "surprise")

### CheckM (Parks et al., Genome Res 25:1043, 2015)
[doi:10.1101/gr.186072.114](https://doi.org/10.1101/gr.186072.114)

**What it computes** [v, PMC4484387]:

- **Markers.** Lineage-specific markers are Pfam/TIGRFAM families single-copy in more than 97% of genomes
  of a lineage, computed at every node of a 5,656-genome tree built from 43 genes.
- **Collocated sets.** Markers within 5 kb of each other in more than 95% of the lineage's genomes are
  merged into collocated sets by chaining shared genes.
- **Completeness.** (1/|M|) Σ_s |s ∩ G_M| / |s|, so every collocated set has equal weight.
- **Contamination.** Extra copies, normalised the same way.
- **Choosing the lineage node.** The node is chosen by simulation (10 kb fragments, 50-100% completeness,
  0-20% contamination) to minimise error. If fewer than 10 tree genes are found, the universal set is
  used.

**Read as surprise:**

- CheckM's markers are clusters with level-1 prevalence of at least 0.97 at one tree node, weighted 0/1.
- Collocated sets exist so that co-located genes, which go missing together, are not counted as
  independent absences. That is exactly why our unit is the cluster.

Our model generalises CheckM in three ways:

- Every prevalence, not only ≥ 0.97, with log-loss weights instead of equal weights.
- Prevalence shrunk across nodes instead of one chosen node.
- Presence surprise (foreign content) as well as absence surprise.

### CheckM2 (Chklovski, Parks, Woodcroft & Tyson, Nat Methods 20:1203, 2023)
[doi:10.1038/s41592-023-01940-w](https://doi.org/10.1038/s41592-023-01940-w)

- **What it computes** [a]: machine-learned completeness and contamination, accurate for novel and
  reduced-genome lineages (Patescibacteria, DPANN), and updatable with new genomes.
- **Use:** our detection model takes c and r from CheckM2.
- **Caveat:** CheckM2 itself learned from gene content, so a bin's surprise and its CheckM2 completeness
  are not independent. A bin missing an unusual amount of expected content gets both low c and fewer
  absence bits. Report absence surprise both with c and with c fixed at 1, the latter as an unconditional
  check.

### GUNC (Orakov et al., Genome Biol 22:178, 2021)
[doi:10.1186/s13059-021-02393-0](https://doi.org/10.1186/s13059-021-02393-0)

**What it computes** [v, PMC8201837]:

- **Gene labels.** Every gene gets its best DIAMOND hit in proGenomes 2.1 and a taxonomy at 7 ranks.
- **Clade separation score.** CSS = 1 - Ĥ(T|C)/Ĥ(T|R): the conditional entropy of taxon labels given
  contig, relative to a contig-size-matched null. CSS > 0.45 at any rank flags a chimera.
- **Reference representation.** RRS (fraction mapped × fraction retained × mean identity) flags genomes
  outside the reference.
- **Contamination measures.** A minority-label fraction and Teff (inverse Simpson index minus 1).
- **Prevalence of chimeras.** 5.7% of GenBank and 15-30% of filtered "HQ" MAG sets look chimeric.

**What we can reuse:** the entropy-with-contig-null construction. Our per-contig presence-surprise
density can be tested for contig structure the same way.

**Where we differ:** GUNC asks whether contigs carry different *lineage labels*. We ask whether a bin's
*content* is improbable for its clade. A same-genus contaminant passes GUNC (species-level CSS is weak at
low contamination [v]) but can carry clusters with zero prevalence in the host species.

### DAS Tool (Sieber et al., Nat Microbiol 3:836, 2018) and MAGpurify (Nayfach et al., Nature 2019)

- DAS Tool: [doi:10.1038/s41564-018-0171-1](https://doi.org/10.1038/s41564-018-0171-1)
- Nayfach et al. (MAGpurify): [doi:10.1038/s41586-019-1058-x](https://doi.org/10.1038/s41586-019-1058-x)

**DAS Tool** [v, README and `DAS_Tool.R`]:

- Picks bins across binners by a single-copy-gene score with weights b (duplicates) and c (megabins,
  default 0.5).
- Keeps selecting while the score is at least 0.5.
- This is the bin-selection rule our surprise-based score should replace or extend.

**MAGpurify** [v, README]:

- Removes contigs flagged by phylo-markers, clade-markers, **conspecific** (contigs that fail to align to
  closely related genomes), tetranucleotide, GC, coverage and known-contaminant modules.
- The conspecific module is a crude, unweighted ancestor of our "surprisingly present" score at contig
  level.

Both are baselines for stage 5.

### Evidence that MAG absences are not missing at random (detection model inputs)

- **Ribosomal protein genes** (Mise & Iwasaki, ISME Commun 2:118, 2022)
  [doi:10.1038/s43705-022-00204-6](https://doi.org/10.1038/s43705-022-00204-6) [a]:
  - Ribosomal protein genes are missing from 20-40% of MAGs with completeness of at least 90%, across
    more than 190,000 MAGs.
  - The cause is codon usage and tetranucleotide-based binning, worst in fast growers.
  - So the detection probability depends on the cluster (its composition), not only on c.
- **Genomic islands and plasmids** (Maguire et al., Microb Genom 6:e000436, 2020)
  [doi:10.1099/mgen.0.000436](https://doi.org/10.1099/mgen.0.000436) [a]:
  - In a simulated 30-genome community with 12 pipelines, 82-94% of chromosomes were recovered but only
    38-44% of genomic islands and 1-29% of plasmid sequence.
  - Accessory and MGE clusters are therefore under-detected in MAGs. Prevalence learned from MAG-heavy
    GTDB species is biased low for them.
- **Contaminant genes in pangenomes** (Majernik, Beaver & Bradley, mSphere 10:e00857-24, 2025)
  [doi:10.1128/msphere.00857-24](https://doi.org/10.1128/msphere.00857-24) [v]:
  - 33 of 86 (38%) differentially prevalent genes in a MAG-based Clostridia pangenome were contaminants,
    despite MAG contamination below 5%.
  - Flagged by gene-to-species abundance correlation (lineage-match AUROC 0.97). Tool: PanSweep.
  - Two consequences for us. Reference prevalence must be robust to rare contaminant presences, which a
    floor on ε handles. And the same correlation idea is our population-level "coverage follows the core"
    test.

---

## 5. Metapangenomics: gene presence from coverage in metagenomes

### Anvi'o metapangenome (Delmont & Eren, PeerJ 6:e4320, 2018)
[doi:10.7717/peerj.4320](https://doi.org/10.7717/peerj.4320)

**What it computes:**

- Pangenome gene clusters combined with read recruitment from TARA metagenomes [a].
- Some core gene clusters sit in hypervariable islands that systematically lack read recruitment.
- Splits gene clusters into "environmental core" and "environmental accessory" by per-gene coverage
  relative to the genome. The exact threshold (we recall a fraction of the genome's median coverage) was
  not re-read; PMC was blocked by a captcha [u].

### MIDAS2 (Zhao et al., Bioinformatics 39:btac713, 2023)
[doi:10.1093/bioinformatics/btac713](https://doi.org/10.1093/bioinformatics/btac713)

**What it computes** [v, `docs/cnv_module.rst` in czbiohub-sf/MIDAS]:

- Copy number = gene mean coverage / median coverage of 15 universal single-copy marker genes.
- A gene is present if copy number ≥ 0.35. Species with mean coverage below 1× are filtered.
- Databases: UHGG (4,644 species, 286,997 genomes) and GTDB r202 (47,893 species) [v, paper].

### PanPhlAn 3 (Beghini et al., eLife 10:e65088, 2021)
[doi:10.7554/eLife.65088](https://doi.org/10.7554/eLife.65088)

**What it computes** [v, `panphlan_profiling.py` defaults]:

- Coverage normalised per family.
- Present if above 0.5, absent if below 0.25, multicopy if above 1.5.
- Strain detection requires at least 2× coverage and a coverage plateau between 0.75 and 1.25.

### inStrain (Olm et al., Nat Biotechnol 39:727, 2021)
[doi:10.1038/s41587-020-00797-0](https://doi.org/10.1038/s41587-020-00797-0)

**What it computes** [v, docs]:

- Per-scaffold and per-gene coverage and breadth.
- **Expected breadth** = 1 - exp(-0.883 × coverage). Actual breadth well below expected breadth means
  only part of the reference is in the sample.
- Microdiversity (nucleotide diversity) per gene.

### What we can reuse and where we differ

**Reuse:**

- These three supply the ingredients of the population-level detection model: normalisation by
  single-copy-core coverage (MIDAS2) and the coverage-to-breadth relation (inStrain).
- The fixed thresholds (0.35; 0.25/0.5) are the baselines our posterior on frequency f must beat.

**Where we differ:**

- **Calls, not posteriors.** All three threshold per sample and ignore the across-sample constancy of the
  ratio f = gene coverage / core coverage. That constancy is the signal we use for binning.
- **No prior on the gene's frequency.** None puts a prior on f from the clade.
- **Reference-based.** All three need a reference pangenome. Our co-assembly bins are their own
  reference.

---

## 6. Module and pathway presence in incomplete genomes

### gutSMASH (Pascal Andreu et al., Nat Biotechnol 41:1416, 2023; NAR web server 2021; gutSMASH 2.0, J Mol Biol 2026)

- gutSMASH: [doi:10.1038/s41587-023-01675-1](https://doi.org/10.1038/s41587-023-01675-1)
- Web server: [doi:10.1093/nar/gkab353](https://doi.org/10.1093/nar/gkab353)
- gutSMASH 2.0: [doi:10.1016/j.jmb.2025.169567](https://doi.org/10.1016/j.jmb.2025.169567)

**What it computes** [a]:

- Rule-based detection of primary metabolic gene clusters (MGCs). The web server covers 41 pathways.
- 19,890 MGCs in 4,240 gut genomes.
- Version 2.0 adds 14 rules (12 MGC types) and comparisons against 15,024 CGR2 clusters.

**Use:** MGC classes as an external label set for catalogue validation (named in the README plan). It
covers gut anaerobes only, so it is of little use for our marine and freshwater bins.

### MetaPathPredict (Geller-McGrath et al., eLife 13, 2024)
[doi:10.7554/eLife.85749](https://doi.org/10.7554/eLife.85749)

**What it computes** [a]:

- Deep-learning models predict KEGG module presence from a genome's gene annotations.
- Robust "within highly incomplete genomes".

**Relation to us:** this is our level-2 idea on curated modules, with a learned imputation of what is
missing. It has no explicit clade prior and no surprise output. Our `docs/background.md` pathway
benchmark is where it should be compared.

### mOTUpan (Buck, Mehrshad & Bertilsson, NAR Genom Bioinform 4:lqac060, 2022)
[doi:10.1093/nargab/lqac060](https://doi.org/10.1093/nargab/lqac060)

**What it computes** [v, paper]:

- **Core likelihood.** A trait is core if present in all genomes, with P(absent | core) = 1 - c_g.
- **Accessory likelihood.** A trait is accessory if drawn independently with frequency |trait|/|T|.
- **Decision.** Classified by log-likelihood ratio. Completeness is re-estimated iteratively as
  |core ∩ g|/|core|.
- **Validation.** False-discovery rates by bootstrap with synthetic rarefied genomes.

**Relation to us:**

- mOTUpan's core likelihood is the special case of our detection model with θ = 1, π = 1, ε = 0.
- Its iterative completeness re-estimation is worth copying when CheckM2 is unreliable, for example for
  small genomes.

### Phylogenetic occupancy model (Mattick, DeMontigny & Delwiche, bioRxiv 2026, v3 2026-05-20)
[doi:10.64898/2026.02.27.708499](https://doi.org/10.64898/2026.02.27.708499)

**What it computes** [a only; bioRxiv rate-limited the full text and Europe PMC (PPR1221547) has no full
text]:

- An ecological occupancy model (latent true presence, imperfect detection) combined with evolutionary
  modelling.
- Gives "a probability distribution over gene presence that accounts for shared information across
  related genomes".
- Jointly estimates genome completeness and the probability that a gene is present but unobserved.
- Outperforms mOTUpan for core-genome inference. Supports single-gene presence and ancestral states.
- Released as a Python package.

**This is the nearest existing method to our level 1 plus detection.**

- **Must do:** read the full text before implementing level 1. If its phylogenetic part is a gain/loss
  Markov chain on a tree, it may scale worse than the taxonomy-rank shrinkage proposed here, but it is
  the right comparison.
- **Where we differ (from the abstract):** our units are clusters, we add contig-edge and depth terms, a
  context level, population frequencies and binning. Its scale and phylogenetic model are unknown.

### Phylogenetic trait prediction priors

**PICRUSt2** (Douglas et al., Nat Biotechnol 38:685, 2020)
[doi:10.1038/s41587-020-0548-6](https://doi.org/10.1038/s41587-020-0548-6):

- Predicts gene-family copy numbers for unsampled tips by hidden-state prediction on a reference tree,
  using castor [u: the method details are from memory].

**castor** (Louca & Doebeli, Bioinformatics 34:1053, 2018)
[doi:10.1093/bioinformatics/btx701](https://doi.org/10.1093/bioinformatics/btx701):

- Comparative phylogenetics on trees with millions of tips [a].
- castor's `hsp_*` and Mk-model fitting are an off-the-shelf tree-based alternative to our rank-based
  shrinkage.

**Mishra & Lercher** (Mol Biol Evol 42:msaf166, 2025)
[doi:10.1093/molbev/msaf166](https://doi.org/10.1093/molbev/msaf166) [a]:

- Presence/absence methods on species trees beat gene-tree reconciliation for HGT inference.
- They benchmark by the tendency of HGT to move neighbouring genes together, which is the same
  co-transfer signal that defines our clusters.

---

## 7. Surprisal as a genome-content or binning score: prior use

**Searches run:**

- Web search: "surprisal OR information content, unexpected gene content, MAG contamination, binning,
  co-occurrence" (extended mode).
- Web search: "phylogenetic prior gene presence absence hierarchical Bayesian taxonomy shrinkage".
- PubMed: gene co-occurrence and pangenome tools (as cited above).

**Result:** no paper scores genomes or bins by -log P(cluster state | clade, context). The closest are:

- CheckM's 0/1 marker sets.
- GUNC's entropy of lineage labels.
- Mise & Iwasaki's "unexpected absence" (descriptive, no score).
- Generic anomaly detection by surprisal: arXiv 2603.09318, not genomic [u].
- SurprisalAnalysis (PMC12707977) [u]. This is "surprisal analysis" in the thermodynamic, Levine sense,
  applied to expression data. It is unrelated despite the name.

**Claim to keep modest:** information-theoretic framing of genome completeness is implicit in every
marker-set method. What is new is the calibrated, clade- and context-conditional probability over all
conserved clusters, and its use as a binning likelihood.

---

## 8. Benchmarks

- **CAMI II** (Meyer et al., Nat Methods 19:429, 2022)
  [doi:10.1038/s41592-022-01431-4](https://doi.org/10.1038/s41592-022-01431-4) [a]:
  - About 1,700 genomes and 600 plasmids and viruses.
  - "Related strains still were challenging for assembly and genome recovery through binning". The
    strain-madness and marine datasets have known contig membership.
- **CAMI III** [v, search snippets of cami-challenge.org only, page not opened]:
  - Longitudinal Human Gut challenges (short and long reads) opened 2026-07-14.
  - The genome-binning challenge supplies reads, gold-standard assemblies and per-sample assemblies.
  - A toy human-gut dataset (10 individuals × 2 samples, public genomes) was released 2026-03-19.
  - No dedicated strain-resolved CAMI III dataset was found. Gold standards for the challenge data are
    presumably released after closing [u].
- **Curated operons:** EcoCyc TUs (in use), SubtiWiki, DOOR/ProOpDB predictions, PADLOC/DefenseFinder
  systems (Spacedust's and Panorama's benchmarks).

---

## Summary table

| Work | Unit | Prevalence prior | Detection model | Context | Scale | Reuse |
|---|---|---|---|---|---|---|
| Spacedust | pairwise conserved loci | none | none | none | 1.3k genomes, quadratic | discovery, p_clu, PADLOC test |
| PPanGGOLiN / PanGBank | families; accessory modules | within-species partition | in mixture model [u] | graph MRF | 4.6k species / 393k genomes | species matrices, projection |
| PANORAMA | rule-defined systems | none | strict/extended/split | rules | >6k genomes | GHGF + GFRR for cross-species identity |
| anvi'o graphs (Henoch) | synteny graph | none | none | none | 29 genomes | shared κ across genera |
| Coinfinder / Goldfinder | gene pairs | none | none | pairwise test, phylo null | one species | feature filter for level 2 |
| Beavan et al. 2024 | genes | none | none | random forest | one species | level-2 design and baseline |
| gLM / Bacformer | genes in context | implicit | none (trained on MAGs) | transformer | 7M contigs / 1.3M genomes | level-2 baseline |
| CheckM / CheckM2 | markers / all genes | ≥97% at one node / ML | collocated sets | none | ~5.7k ref genomes | c, r inputs; special case |
| GUNC | all genes, lineage labels | none | contig-size null | none | proGenomes | entropy-with-null construction |
| mOTUpan | genes | core vs accessory | c_g | none | 1000s genomes | completeness re-estimation |
| Phylo occupancy (Mattick 2026) | genes | phylogenetic | joint completeness | none | ? | **level-1 competitor** |
| MIDAS2 / PanPhlAn / inStrain | genes in metagenomes | none | coverage thresholds | none | 100k+ genomes | detection ingredients |
| MetaPathPredict | KEGG modules | none | learned | DL classifier | ? | pathway-level comparison |
| DAS Tool / MAGpurify | SCGs / contigs | none | none | none | | stage-5 baselines |
