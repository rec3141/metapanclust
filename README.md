# metapanclust

Experiments in **cluster-level pangenomics of metagenomic populations**: treat conserved gene
clusters (operons, cassettes, other co-located gene sets) as the units of genome content, learn how
often each cluster occurs within species pangenomes and along the GTDB tree, and ask of any genome or
population bin: *is the model surprised that this cluster is present, or absent?*

Status: early experiments. Nothing here is stable.

## Why clusters, why populations

- **Clusters are observed, pathways are inferred.** Whether a set of genes sits together on a contig
  depends on the assembly and annotation, not on a curated pathway model or inference rules. Pathways can
  be linked back afterwards.
- **The species population is a sensible unit for metagenomic ecology.** A co-assembly bin is a
  population. An accessory cluster carried by some of its cells shows up at a fraction of the core's
  coverage, and that ratio, per sample, estimates how many cells carry it.

## Surprise

`surprise = -log2 P(observed state | context)` for each cluster in a genome or population bin, reported
both ways (surprisingly present, surprisingly absent). It is computed at two levels of context:

1. **Taxonomy:** cluster prevalence within the species pangenome, shrunk up the GTDB tree
   (species → genus → family → …, hierarchical beta-binomial), so sparse clades borrow from their
   parents.
2. **Taxonomy and the rest of the genome:** predict each cluster's presence from the other clusters
   plus clade (one masked logistic model per cluster to start), so co-occurrence counts. Phylogeny-aware
   null models in the spirit of Goldfinder keep shared ancestry from looking like linkage.

Absence is discounted by completeness (CheckM2), by contig ends (a cluster that would sit across a
break is unknown, not absent) and by read depth (detection limit). Clusters on MGEs (geNomad) are scored
separately: horizontal transfer is supposed to be surprising.

## Plan

| Stage | What | Check against |
|---|---|---|
| 1 | **Cluster catalogue** over GTDB species representatives (Spacedust or an MMseqs2/Foldseek gene-neighbourhood graph on eggNOG OGs) | curated operons (RegulonDB/EcoCyc, SubtiWiki), antiSMASH/gutSMASH cluster classes, PADLOC defence systems |
| 2 | **Prevalence and taxonomy surprise**: species → tree shrinkage; per-cluster co-location propensity | held-out clades |
| 3 | **Context surprise**: masked per-cluster models | clusters hidden in held-out genomes |
| 4 | **Population frequencies**: per-sample cluster/core coverage ratios in co-assemblies | simulated strain mixtures; CAMI III strain data |
| 5 | **Binning by surprise**: score bin merges and contig moves by the change in surprise; accept accessory contigs whose coverage follows the core at a lower ratio | CAMI III gold standard (contig membership known) |

Inputs come from the [danaSeq](https://github.com/rec3141/danaSeq) `mag_analysis` pipeline (Bakta gene
coordinates, eggNOG OGs, KOs, CheckM2, GTDB-Tk, geNomad, per-sample depths, 95% ANI dereplication).

## Related work

- **AMPHORA** (Wu & Eisen 2008) → universal single-copy markers; **CheckM** (Parks et al. 2015)
  lineage-specific and co-located marker sets: completeness/contamination are a narrow form of
  surprise over single-copy genes.
- **Spacedust** (Zhang, Mirdita & Söding, Nat Methods 2025): de novo conserved gene clusters across
  genomes with Foldseek; 72,843 clusters covering 58% of genes in 1,308 genomes.
- **PPanGGOLiN** / panModule (conserved co-located accessory modules in pangenome graphs),
  **Panorama** (systems from conserved genomic context across species), synteny-aware pangenome graphs
  (Henoch et al. 2026, anvi'o).
- **Coinfinder**, **Goldfinder**: phylogeny-aware gene co-occurrence and avoidance.
- **gLM** (Hwang et al. 2024), **Bacformer** (2025): language models of genes in genomic context.
- Metapangenomics: **anvi'o** (Delmont & Eren 2018), **MIDAS2**, **PanPhlAn**, **inStrain**.
- Pathway-level background: [docs/background.md](docs/background.md).

## Layout

```
docs/       design notes and results
scripts/    experiment scripts (small, one task each)
notebooks/  exploration
```

Data, databases and outputs stay outside the repo (see `.gitignore`). Licensed reference data
(MetaCyc/BioCyc) must never be committed.
