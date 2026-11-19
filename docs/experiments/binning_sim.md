# Binning by surprise: simulated contamination and incompleteness (2026-11-18)

Question (README stage 5): can a species pangenome tell a bin that has picked up another species
from one that has picked up more of its own species, and how does that compare with a CheckM-style
single-copy-marker estimate?

**Short answer.** Yes, for contaminants from outside the species, with one caveat about contig length.
Accept a merge when the added contigs cost no more bits per Mb than the bin already does. That rule
separated foreign-species additions from same-species additions with AUC 1.0 at every contamination
level (5-40%), in both species and for related and unrelated contaminants. Marker duplication ranks
the two the wrong way round (AUC 0.0-0.4): a same-species merge duplicates markers and a contaminant
mostly does not. Short contigs (< 10 kb) are where the rule fails, because a genome's own short
contigs carry its unique genes.

## Setup

Code: `experiments/binning_sim/` (`01_build.sh` pangenomes, projections and ANI; `02_simulate.py`
bins and scores; `03_analyze.py` tables and figures; `04_run.sh` runs 02 and 03). Data and full
tables: `/data/scratch/metapanclust_data/binning_sim/` (`results/summary.md`, `results/sim_*.tsv`,
figures `results/*.png`). The whole run takes about 25 CPU-minutes, mostly PPanGGOLiN.

- **Held out:** the 6 most complete public genomes of each species, chosen by PPanGGOLiN persistent
  completeness. Pairwise skani ANI among them is 96.4-99.8% (UBA9145) and 98.6-99.6% (MAG-120531).
  Each pangenome was rebuilt without them (`ppanggolin all`, defaults):

  | Pangenome | Training genomes | Families | Persistent | Single-copy markers* | Modules |
  |---|---|---|---|---|---|
  | UBA9145 sp001438145 | 48 | 5,496 | 1,813 | 1,732 | 90 |
  | MAG-120531 sp000173115 | 63 | 3,853 | 1,459 | 1,448 | 38 |

  *Markers are persistent families that are multicopy in < 5% of training genomes (PPanGGOLiN
  `dup_margin`).
- **Gene assignment:** each held-out genome and contaminant was annotated against the reduced
  pangenome with `ppanggolin projection`, in two modes:
  - strict: 80% identity, 80% coverage (the PPanGGOLiN default);
  - relaxed: 30% identity, 50% coverage. This is a rough stand-in for CheckM's HMMs, which find
    markers across phyla.

  A gene that maps to no family is *foreign*.
- **Contaminants (B):**
  - unrelated: the other species' 6 held-out genomes (Gammaproteobacteria vs Bacteroidota);
  - related ("near"): Chesterfield DAS Tool bins.
    - For UBA9145: Marimicrobium, Halieaceae, same order Pseudomonadales; CheckM2 92% complete.
    - For MAG-120531: MAG-121220-bin8 and Hel1-33-131, same family Flavobacteriaceae; CheckM2 61-62%
      complete.

  No congeneric or same-genus contaminant was available.
- **Simulated bins:** 50-kb pseudo-contigs (real contigs cut into windows, genes kept with their
  window); 5 random draws per condition and genome.

  | Kind | Contents |
  |---|---|
  | clean | the held-out genome alone |
  | incomplete | 10-50% of the genome's length removed |
  | contam_B | B contigs added at 5/10/20/40% of the base's length |
  | strain | contigs of another held-out genome of the same species, at the same fractions |
  | contam_B_inc / strain_inc | the same two additions to a base reduced to 70% |
  | single-contig moves | one contig offered to the 70% base: from its own removed part, from another strain, or from B |

  That makes about 3,600 bins over 4 pangenome × contaminant pairs.
- **Scores:**
  - comp: completeness, the fraction of markers present.
  - dup / dup_rel: duplication, extra marker copies per marker, from strict or relaxed mapping.
  - foreign: the fraction of genes that are foreign.
  - S, total surprise in bits over all pangenome families. Prevalence is p = (k+1)/(n+2) over the n
    training genomes, and the terms are:
    - present family: -log2 p;
    - absent family: -log2(1 - p·comp);
    - each foreign gene: -log2 1/(n+2), about 5.6-6.0 bits.

    Optional extras: a multicopy term, and the same total over modules (S_mod).
  - Merge rule: **excess = (S(bin+piece) - S(bin)) / Mb added - S(bin)/Mb**. Reject the merge if
    excess > 0.

## Results

The full per-condition tables are in `results/summary.md`. Figures:
`experiments/binning_sim/fig/merge_uba9145_mag120531.png`,
`fig/merge_mag120531_near_mag120531.png` and `fig/incomplete_uba9145.png`.

### 1. Strict marker duplication cannot see contamination from outside the species

At 80% identity, a contaminant's genes almost never land in the pangenome's families: 2 of about
1,900 genes for each Bacteroidota genome against UBA9145. Strict duplication therefore stays near the
clean baseline while the foreign fraction climbs. Contamination at 40% of length:

| Pangenome / contaminant | dup | dup_rel | foreign | dup after a 40% same-species merge |
|---|---|---|---|---|
| UBA9145 / Bacteroidota | 0.004 (clean 0.004) | 0.136 | 0.32 (clean 0.03) | 0.41 |
| UBA9145 / Halieaceae | 0.012 | 0.224 | 0.33 | 0.41 |
| MAG-120531 / Gammaproteobacteria | 0.003 (clean 0.003) | 0.106 | 0.32 (clean 0.05) | 0.42 |
| MAG-120531 / Flavobacteriaceae | 0.043 | 0.371 | 0.33 | 0.42 |

Relaxed (CheckM-like) duplication does detect contamination: AUC 0.97-1.0 against clean and
incomplete bins, already at 5%. It underestimates contamination from another phylum, though: about
a third of the true value, because only about 23% of the persistent families have a 30%-identity
homolog there. It is closer for the same-family contaminant.

Neither duplication measure can tell contamination from a same-species merge. Against strain merges
at equal added length the AUC is 0.0 for strict and 0.0-0.39 for relaxed duplication: a strain merge
duplicates more markers than a contaminant does.

### 2. Surprise separates foreign from same-species content; the merge rule works

Total surprise per Mb at 0 / 5 / 10 / 20 / 40% added (UBA9145, unrelated B):

| Bin | Bits/Mb |
|---|---|
| contam_B | 727 → 971 → 1154 → 1480 → 1990 |
| strain | 727 → 700 → 678 → 651 → 594 |
| strain_inc (70% base, 1144 bits/Mb) | 1099 → 1042 → 965 → 838 |

So adding B contigs raises surprise per Mb and adding same-species contigs lowers it, for both a
complete and an incomplete base. The other three pangenome × contaminant pairs look the same (see
`summary.md`).

Merge rule (excess > 0 means reject):

| Pangenome / contaminant | AUC B vs same-species, all levels | B merges rejected | Same-species merges rejected |
|---|---|---|---|
| UBA9145 / Bacteroidota | 1.0 | 100% | 0-7% |
| UBA9145 / Halieaceae | 1.0 | 100% | 0-10% |
| MAG-120531 / Gammaproteobacteria | 1.0 | 100% | 7-33% |
| MAG-120531 / Flavobacteriaceae | 1.0 | 100% | 0-17% |

The B pieces cost about 5,000-5,900 bits/Mb. The bins themselves cost 430-1,450 bits/Mb.

Adding the multicopy term to S makes separation worse (AUC 0.74-1.0), because it penalizes the
duplication that same-species merges create.

The higher rejection rate for MAG-120531 same-species merges comes from one held-out genome (see 5).

### 3. Single-contig moves: reliable at ≥ 10 kb, unreliable below

| Pangenome / contaminant | AUC, < 10 kb | AUC, 10-30 kb | AUC, 30-50 kb |
|---|---|---|---|
| UBA9145 / Bacteroidota | 0.94 | 1.0 | 1.0 |
| UBA9145 / Halieaceae | 0.95 | 0.99 | 1.0 |
| MAG-120531 / Flavobacteriaceae | 0.80 | 0.96 | 1.0 |
| MAG-120531 / Gammaproteobacteria | 0.52 | 0.98 | 1.0 |

The bin's *own* contigs under 10 kb were rejected 26-35% of the time (UBA9145) and 69-73%
(MAG-120531). Short contigs in MAGs are where the genome's unique content sits: islands and mobile
elements, which also fragment assemblies. A surprise-only rule would strip that accessory content
from bins. This supports the README's plan to accept such contigs on coverage (tracking the core at a
lower ratio) rather than on surprise.

### 4. Incompleteness: surprise adds nothing beyond marker completeness

Over clean and incomplete bins, marker completeness tracks true completeness (Spearman ρ 0.98).
Absence surprise without the completeness discount (c = 1) per Mb tracks it equally well (ρ -0.98).
Both are driven by the same high-p families, so surprise gives no new information on completeness.

With the discount, total surprise per Mb still rises as bins get less complete (727 → 1,312 bits/Mb
at 50% removed, UBA9145). As an *absolute* bin-quality score it therefore confounds incompleteness
with contamination: AUC 0.25-0.34 for 5% contamination against clean and incomplete bins. Content
surprise per gene (absence term left out) avoids this: AUC 0.80-0.91 at 5% and 0.92-1.0 at 20%.

The merge rule is not affected, because it compares a bin with itself.

Module surprise is uninformative here (AUC 0.12-0.37 for contamination). There are few modules
(90 and 38), and foreign genes never complete one, so module surprise per Mb *falls* with
contamination.

### 5. By-product: a public MAG-120531 genome with 18% foreign genes

GCA_963976715.1 (held out) has 17.8% of its genes outside the reduced pangenome; the other 11
held-out genomes have 0.7-6.4%. Some details:

- **Contigs:** 51 small contigs (201 kb, 249 genes) carry no pangenome gene at all.
- **GC:** these contigs have 37.5% GC, against 38.1% for the rest.
- **Homologs:** only 25 of the 249 genes have even a 30%-identity homolog in the pangenome.
- **Size:** the genome is 2.40 Mb, against 1.94-2.11 Mb for the other five held-out genomes.

PPanGGOLiN's persistent-family contamination for this genome is 0.2%. This is either contamination
from a GC-matched organism or a block of mobile elements. geNomad and per-contig coverage would
decide; the surprise score cannot. Its contigs are what get same-species merges rejected in (2) and
(3), and they lower the foreign-fraction AUC for MAG-120531 at 5-10% contamination to 0.86.

## Limits

- **Few species.** Only two species pangenomes; the contaminants are the other species and three
  related Chesterfield bins (same order, or same family).
  - The hard case is untested: a congeneric or same-genus contaminant at 80-95% ANI, whose genes
    partly map at 80% identity. Strict mapping makes anything outside the species look foreign, so
    the AUC of 1.0 here is close to the easy limit.
  - The near contaminants are themselves 61-92% complete MAGs.
- **Pangenomes built from MAGs.** All training genomes are MAGs, so contamination in them raises p for
  foreign families. GCA_963976715.1 shows such genomes exist in the GTDB species set.
- **Held-out set.** The held-out genomes are the most complete ones and often come from the same
  studies. Strain diversity is limited (ANI ≥ 96.4%; ≥ 98.6% for MAG-120531).
- **Simulated contigs.** Pseudo-contigs are 50-kb windows, not real MAG contig-length distributions.
  Fragmentation enters only through the single-contig test.
- **Model simplifications.** Families are scored independently (no context model, stage 3) and with
  no coverage. Every foreign gene costs the same, so "foreign fraction" and S carry much the same
  signal for contamination. S adds a common scale (bits) that also covers absences and rare families.
- **Projection.** Projection assigns each gene its best family only. The relaxed mode is a crude
  substitute for CheckM marker HMMs, not CheckM itself.

## For grex

Nothing needed: everything ran locally under `nice` with 1-2 threads. The obvious next step would be a
cluster job: repeat this over many GTDB species with congeneric contaminants (same genus,
85-95% ANI), using Chesterfield/CAMI bins with known contig membership.
