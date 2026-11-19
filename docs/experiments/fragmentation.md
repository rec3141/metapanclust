# Fragmentation and cluster detection probability (2026-11-18)

**Result.** Breaking a complete genome at random into contigs with Chesterfield bin-like lengths leaves
**83% of curated E. coli operons and 81% of six-genome Spacedust clusters intact** (mean over the 93
bins, before completeness). At the modal bin N50 of 10-20 kb only 74% / 70% stay intact; at N50 > 100 kb,
97-99%. Losses grow with span: at N50 10-20 kb, 87% of clusters under 2 kb stay intact but only 26% of
those spanning 10-20 kb. The probability of seeing a cluster intact has an exact closed form under random
breaks, and it matches the simulation to 0.5 percentage points (mean abs. error per bin x span class).
In the real Chesterfield bins, conserved clusters touch contig ends **less often than random breaks
predict**: the excess over the no-truncation null is 25% of the model's. Contig ends are enriched for
hypothetical, non-conserved genes, so assemblies break in variable DNA more than inside conserved clusters.
The random-break formula is therefore conservative (it overstates splitting), which is the safe direction
for scoring absence.

Code: `experiments/fragmentation/` (`detection.py` is the model; `simulate.py`, `analyze_sim.py`,
`real_bins.py`, `end_profile.py`). Data and outputs: `/data/scratch/metapanclust_data/fragmentation/`.

## Detection model (for the surprise model)

The model takes the bin's contig lengths L_i, its completeness c (CheckM2) and the cluster span s in bp
(first gene start to last gene end in a reference or a sister genome), and returns:

```
G        = sum(L_i) / c                                  # estimated true genome size
f(s)     = sum_i max(L_i - s + 1, 0) / sum(L_i)          # P(span inside one contig | on assembled DNA)
P_intact = c * (1 - beta * (1 - f(s)))                   # beta = 1: breaks independent of the cluster
P_lost   ~ (1 - c) * f(s) + (1 - c)^2 * (1 - f(s))       # approximation
P_partial = 1 - P_intact - P_lost
surprisingly absent = -log2(1 - p_prior * P_intact)
```

With beta = 1 this is exact if the cluster's position is independent of where the assembly broke, and
it holds whichever pieces were lost (short contigs first or at random), because only the binned contigs
and G enter. Code: `detection.py`: `p_intact`, `p_lost`, `detection`, `absent_logprob`, all vectorised
over spans. A cluster seen truncated at a contig end should be scored as present-but-unresolved, never
as absent.

**Beta.** Use 1 (conservative) for absence scores until a better calibration exists. The Chesterfield
end-touch excess suggests about 0.25 for conserved clusters. That estimate is a lower bound on true
splitting, because truncated fragments may also fail Spacedust's significance test and go unseen.

Rule of thumb, if contig lengths are unavailable: P_intact ~ c * exp(-s / (0.64 * N50)). Over the
simulation it fits only slightly worse than the exact form (per-cluster MAE 0.047 vs 0.046, mostly Monte
Carlo noise).

Example (P_intact, beta = 1, spans 2 / 5 / 10 / 20 kb):

| Bin | N50 | c | 2 kb | 5 kb | 10 kb | 20 kb | -log2 absent, 5 kb, p=0.95 |
|---|---|---|---|---|---|---|---|
| dastool-maxbin_115 | 10 kb | 0.79 | 0.58 | 0.34 | 0.23 | 0.12 | 0.57 |
| dastool-semibin_019 | 19 kb | 0.73 | 0.60 | 0.45 | 0.28 | 0.12 | 0.79 |
| dastool-vamb_623 | 67 kb | 0.97 | 0.91 | 0.84 | 0.73 | 0.57 | 2.30 |
| dastool-lorbin_029 | 118 kb | 0.81 | 0.80 | 0.77 | 0.73 | 0.66 | 1.91 |
| dastool-vamb_067 | 527 kb | 0.79 | 0.79 | 0.78 | 0.77 | 0.74 | 1.94 |

Without the fragmentation term, a 5 kb cluster missing from maxbin_115 would score -log2(1 - 0.95*0.79)
= 2.0 bits instead of 0.57.

## (1) Simulation

- **Inputs.** Six complete genomes (E. coli, B. subtilis, P. aeruginosa, M. tuberculosis, Synechocystis,
  M. jannaschii; Bakta CDS). 789 curated E. coli multi-gene transcription units (chains of same-TU pairs
  from `eco_pairs.tsv`; median 2 genes, 2.7 kb). 2,285 unique Spacedust cluster instances from the
  six-genome run (median 2 genes, 3.0 kb).
- **Fragmentation.** Each replicon is cut into pieces drawn with replacement from one bin's contig
  lengths (stationary start). Each piece is kept with probability equal to the bin's CheckM2
  completeness. 93 bins x 50 replicates.
- **Scoring.** A gene counts only if it lies wholly on a kept piece. *Intact* means all genes are on one
  kept contig; *lost* means no gene is kept; *partial* is the rest; *detected* means at least 2 of the
  cluster's genes are on one contig.
- **Bins.** Bin N50 ranges from 4.7 kb to 527 kb (median 34 kb; median completeness 0.66).

Fraction intact from fragmentation alone (both sets pooled), by bin N50 x cluster span:

| N50 \ span | <2 kb | 2-5 kb | 5-10 kb | 10-20 kb | >20 kb |
|---|---|---|---|---|---|
| <10 kb (9 bins) | 0.78 | 0.51 | 0.21 | 0.06 | 0.01 |
| 10-20 kb (25) | 0.87 | 0.72 | 0.49 | 0.26 | 0.08 |
| 20-50 kb (26) | 0.94 | 0.86 | 0.73 | 0.55 | 0.29 |
| 50-100 kb (14) | 0.96 | 0.92 | 0.84 | 0.73 | 0.51 |
| 100-250 kb (14) | 0.99 | 0.97 | 0.94 | 0.89 | 0.77 |
| >250 kb (5) | 1.00 | 0.99 | 0.98 | 0.97 | 0.93 |

By gene count at N50 10-20 kb: 2 genes 0.79, 3 genes 0.69, 4-5 genes 0.59, 6-9 genes 0.46, 10+ genes 0.27.
Curated operons and Spacedust clusters behave the same at equal span (2-5 kb: 0.82 vs 0.82; 10-20 kb:
0.53 vs 0.53). Span, not cluster type, is what matters.

With completeness (Spacedust set):

| N50 | mean c | intact | partial | lost | >= 2 genes detected |
|---|---|---|---|---|---|
| <10 kb | 0.45 | 0.23 | 0.26 | 0.52 | 0.34 |
| 10-20 kb | 0.55 | 0.38 | 0.19 | 0.43 | 0.47 |
| 20-50 kb | 0.45 | 0.39 | 0.08 | 0.54 | 0.43 |
| 50-100 kb | 0.59 | 0.53 | 0.06 | 0.41 | 0.56 |
| 100-250 kb | 0.81 | 0.78 | 0.03 | 0.19 | 0.80 |
| >250 kb | 0.87 | 0.86 | 0.01 | 0.13 | 0.87 |

## (2) Model vs simulation

Comparison over 465 bin x span-class cells:

| Quantity | MAE | Max abs. error | Bias | r |
|---|---|---|---|---|
| P_intact, fragmentation only | 0.005 | 0.06 | -0.0001 | 0.9995 |
| P_intact, with completeness | 0.006 | 0.05 | -0.0004 | 0.9995 |
| P_lost (approximation) | 0.023 | 0.18 | -0.019 | 0.994 |

P_lost errs where genes cut by a break are lost, which the approximation ignores. It matters only for
low-completeness bins.

## (3) Real MAGs: conserved clusters at contig ends

The test set is 139,041 Spacedust cluster instances in the 93 Chesterfield bins (median 3 genes,
3.6 kb). An instance "touches an end" if it includes the contig's first or last CDS. The null assumes
uniform placement within the contig and no truncation, 2/(n-w+1). The model is the simulation above,
using each bin's contigs.

| Bin N50 | Instances | Observed | Null (no truncation) | Model (random breaks) |
|---|---|---|---|---|
| <10 kb | 11,281 | 0.360 | 0.337 | 0.538 |
| 10-20 kb | 36,126 | 0.229 | 0.201 | 0.348 |
| 20-50 kb | 31,401 | 0.113 | 0.085 | 0.168 |
| 50-100 kb | 21,928 | 0.078 | 0.053 | 0.108 |
| 100-250 kb | 28,317 | 0.030 | 0.020 | 0.041 |
| >250 kb | 9,988 | 0.008 | 0.005 | 0.015 |

- **Overall.** 13.3% of instances touch an end, against 11.1% under the null. The excess over the null
  is 0.25 of what random breaks predict, and the observed/model ratio is 0.67. Per bin, observed
  correlates with the model at r = 0.89.
- **Why.** Contig ends are depleted of conserved genes. Over 905 contigs of 40+ CDS, the share of CDS in
  any Spacedust cluster is 0.44 at the terminal gene, 0.54-0.59 at genes 1-3, about 0.61 at genes 5-10
  and 0.67 in the interior. The share of hypothetical proteins runs 0.76, then 0.64-0.58, then about
  0.52, then 0.45. The gradient over genes 1-10 is not explained by partial terminal ORFs: breaks fall
  in variable (island/MGE/strain-variable) DNA.
- **Caveat.** The 0.25 mixes two effects: fewer real splits, and truncated fragments that Spacedust no
  longer calls significant. Separating them needs a planted test: cut bins' contigs further, rerun
  Spacedust, and see what fraction of truncated clusters is still recovered.

## Limitations and next steps

- Contamination inflates sum(L_i) slightly; for high-contamination bins, pass contig lengths minus
  flagged contaminant contigs, or G = genome size from the species pangenome.
- The span s must come from the cluster catalogue: the median span across genomes carrying it, or the
  catalogue representative.
- Next: calibrate beta per cluster class (core vs MGE-adjacent: MGE clusters sit nearer breaks, so their
  beta should be near or above 1) using the same end-touch test on catalogue loci.
