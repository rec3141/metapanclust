# Level-1 surprise over catalogue families, Chesterfield 93 bins (2026-11-18)

This is the first cross-species surprise over the cluster catalogue (`catalogue.md`). It uses only the
Chesterfield co-assembly: 93 bins x 7,178 families = 667,554 bin-family cells, 18,748 of them present
(14,764 intact, 3,984 seen only at contig ends). Each bin is scored against a prior estimated from the
**other 92 bins**: a hierarchical beta-binomial along the GTDB-Tk lineage, observed through the
fragmentation detection model.

## Result

- **Taxonomy shrinkage cuts held-out log loss by 16% against a flat prior**: 0.1864 to 0.1561 bits per
  cell, a total of 124,444 to 104,181 bits. On families seen in at least 2 bins the cut is 25% (0.297 to
  0.223). On present cells it is 5.36 to 4.36 bits.
  - Cross-fitting the concentrations (tuned on even-numbered families, scored on odd and vice versa)
    gives the same 0.1561, so the 7 tuned κ values do not overfit.
  - Fitted κ: phylum, class, order and family 3; genus 30; species and domain uninformative. A clade
    contributes about as much as 3 bins from its parent, which is strong shrinkage. The genus adds
    little beyond the family, because most genera here have 1-2 bins.
  - The fully pooled ("hard") prior, which takes the nearest non-empty clade frequency (κ = 0.001), scores
    0.1648. It is better than flat but 6% worse than shrinkage.
- **The detection model is worth about as much as the taxonomy.** Scored as two states (seen / not
  seen):

  | Prior | With detection | Without detection (κ refitted) |
  |---|---|---|
  | flat | 0.167 | 0.180 |
  | taxonomy | **0.137** | 0.155 |

  Taxonomy plus detection is 24% below flat without detection.
- **Leakage warning.** Running the EM detection correction of the training counts with the hierarchical
  prior itself (`model_spec.md` §3.4 A as written) gives 0.110, which looks 30% better than the honest
  number. It is a leak. A held-out bin's own presence enters its clade-mates' posteriors and comes back
  through them: a family seen in TMED14 alone got p̂ = 0.14 for TMED14. In the main fit the training
  posteriors use the flat prior instead, where the held-out bin weighs 1/92.
  **§3.4 A needs a per-bin refit (or a leave-two-out posterior) for its leave-one-out to hold.**
- **Calibration** (predicted P(seen) against observed) is good below 0.05 and somewhat underconfident
  above 0.2:

  | Predicted P(seen) | Observed rate |
  |---|---|
  | 0.32 | 0.41 |
  | 0.61 | 0.70 |
  | 0.84 | 0.87 |

- **Core operons are unsurprising where present.** The S10-spc-alpha superoperon (CF00000) costs a median
  0.06 bits in its 85 bins. When it is absent it costs 2.0 bits (max 2.85). Absence surprise is low across the board. Only 39
  absent states in the whole dataset exceed 2 bits (max 2.89). The ceiling comes from detection: the
  median P_intact of a missing 11-17 kb operon is 0.07-0.21.
- **Contamination shows up in absences, not in presences.**
  - Absence bits per expected Mb correlate with CheckM2 contamination (Spearman 0.66, partial 0.65 given
    completeness). Among bins at least 70% complete, absence separates ≥ 10% contaminated bins from the
    rest with AUC 0.85, against 0.77 for presence. 7 of the 10 most absence-surprising bins are ≥ 10%
    contaminated (base rate 18/93).
  - The likely mechanism: a contaminant's markers inflate CheckM2 completeness, so expected content that
    is really missing gets scored as if it should have been detected.
  - Presence bits per Mb do not track contamination over all bins (ρ = 0.12, p = 0.25).
- **Misbinned contigs carry more presence surprise.**
  - Loci on coverage_test's off-population contigs score a median **3.6 bits**, against 1.6 on
    core-like contigs (MWU p = 1e-22). Without TMED14 the medians are 4.1 vs 1.8.
  - 42% of off-population loci exceed 5 bits, against 24% of core-like loci. On no-reads contigs it is 51%.
  - **But** accessory-like contigs score almost as high (38% > 5 bits), and the within-bin paired test is
    weak: 14 of 23 bins, median +0.36 bits, Wilcoxon p = 0.06. Presence surprise points to "not this
    clade's usual content". It cannot by itself tell contamination from accessory DNA.
- **TMED14 (maxbin_102) ranks 11th of 93 on presence bits per Mb.** Half of its 429 loci sit on
  off-population contigs, scoring 3.7 bits on average (89 of 212 above 5 bits). Its core-like loci score
  2.7 (39 of 194 above 5 bits).
- **UBA9145 (maxbin_040) ranks only 57th.** Within the bin, its accessory-like (4.0 bits) and
  off-population (3.7) loci still score above its core-like ones (2.5).
- **Half the presence surprise is singleton families.** Families found in one bin hold 28,333 of 58,573
  presence bits (48%) and 4,121 of the 5,366 present states above 5 bits. 2,967 of those 5,366 are
  all-"hypothetical protein" loci. A family seen in no other bin costs about 5-10 bits by construction.

## Method

Code: `experiments/catalogue_surprise/` (`run_all.sh`, about 6 min at 2 threads). Data:
`/data/scratch/metapanclust_data/catalogue_surprise/`.

1. **Observations** (`build.py`). The state of family k in bin j is:
   - **I**: any of its loci in j does not touch a contig end (`at_contig_end` in `chesterfield_members.tsv`);
   - **P**: present, but every locus touches a contig end;
   - **A**: absent.

   Span s_k is the median bp span (first gene start to last gene end, from `fam30/genes.tsv`) of the
   family's loci that are not at a contig end. If all are at an end, the longest is used (median 2.0 kb).
   Contig lengths come from the bin FASTAs; completeness and contamination from CheckM2.
2. **Detection** (`model_spec.md` §2, `fragmentation/detection.py`, β = 1):
   - Given z = 1: P(I) = P_intact, P(A) = P_lost, P(P) = the rest.
   - Given z = 0: false presence ε_jk = 1e-4 + 1 - exp(-r_j G_j ρ_k), where r_j is CheckM2
     contamination, G_j = Σ L / c, and ρ_k is the loci of k per bp in the other bins. Of ε,
     P(I) = ε f_j(s_k) and P(P) = ε (1 - f).
   - No composition-loss λ_k and no assembly term a_j.
3. **Prior, leave-one-bin-out** (`surprise.py`).
   - Ranks: all → domain → phylum → class → order → family → genus → species. An unnamed rank (`g__`,
     `s__`) is a node of its own for that bin, so its prior falls back to the parent.
   - Root: Beta(a0, b0) by method of moments over families.
   - Each clade: p̂_t = (X_t − x_j + κ_ℓ p̂_parent) / (N_t − 1 + κ_ℓ), where x_j is the held-out bin's
     count.
   - Training counts x_g = P(z_g = 1 | O_g) under the flat leave-one-out prior, iterated 3 times (see the
     leakage warning).
   - κ_ℓ is chosen per rank on a log grid (0.03 to 1e4), by coordinate descent on held-out log loss.
4. **Surprise** (`analyze.py`).
   - Absent: s = −log2 P(A | p̂). P(A | z=0) is about 1, so this is already content-only.
   - Present: s = −log2[P(O | p̂) / P(O | z=1)] = −log2(p̂ + (1 − p̂) b/a).
     - This is 0 when p̂ = 1 and −log2 p̂ when false presence is negligible.
     - It removes the detection cost that the full −log2 P(O) puts on every observation. With that
       cost left in, presence bits per Mb track low completeness (ρ = −0.66 in the first version).
   - Per bin: presence bits per Mb of bin, and absence bits per expected Mb (Σ_k p̂ s_k), as in §4.2.
5. **Variants.**
   - r_j = 0 (contamination ignored in ε) changes per-bin presence bits per Mb by a median 3.7% (max
     14%, `S_pres_r0_per_Mb`).
   - Counts not corrected for detection: taxonomy 0.172 vs flat 0.197 (`logloss_none.tsv`).

## Log loss (bits per bin-family cell, held out)

| Model | All 667,554 cells | Families in ≥ 2 bins | Present cells | Absent cells |
|---|---|---|---|---|
| flat (root only) | 0.1864 | 0.2967 | 5.36 | 0.0371 |
| taxonomy, fully pooled (κ = 0.001) | 0.1648 | 0.2387 | 4.59 | 0.0370 |
| **taxonomy, shrunk (fitted κ)** | **0.1561** | **0.2227** | **4.36** | **0.0346** |
| taxonomy, shrunk, cross-fitted κ | 0.1561 | 0.2227 | 4.36 | 0.0346 |
| *leaky: EM with hierarchical prior* | *0.1096* | *0.1608* | *2.87* | *0.0300* |

## Core operons

Detection-free "bits present" per family (medians over bins):

| Family | Present in | p̂ where present | Bits present | Absent in | p̂ where absent | P_intact where absent | Bits absent (max) |
|---|---|---|---|---|---|---|---|
| CF00000 S10-spc-alpha | 85 | 0.95 | 0.06 | 8 | 0.93 | 0.07 | 2.01 (2.85) |
| CF00009 L11-L10-L12-rpoBC | 48 | 0.70 | 0.48 | 45 | 0.62 | 0.21 | 0.60 (2.30) |
| CF00011 uL13-uS9 | 52 | 0.70 | 0.43 | 41 | 0.74 | 0.45 | 0.67 (2.34) |
| CF00030 alpha operon | 37 | 0.74 | 0.43 | 56 | 0.42 | 0.43 | 0.46 (1.84) |
| CF00033 atp F0F1 | 32 | 0.60 | 0.68 | 61 | 0.28 | 0.30 | 0.26 (2.39) |
| CF00005 nuo | 41 | 0.73 | 0.45 | 52 | 0.37 | 0.20 | 0.41 (2.45) |
| CF00017 groESL | 38 | 0.61 | 0.56 | 55 | 0.47 | 0.35 | 0.38 (1.43) |
| CF00014 dcw | 42 | 0.68 | 0.51 | 51 | 0.51 | 0.13 | 0.46 (2.17) |

- **Only CF00000 behaves like a universal family.** For the others p̂ stays at 0.6-0.75 even where present,
  for three reasons:
  - the catalogue splits these systems across families (3 families hold 80% of atp and of nuo genes);
  - operons break at contig ends;
  - 56% of bins are < 70% complete.

  For atp, 22 of the 61 bins without CF00033 have ≥ 3 ATP synthase genes elsewhere. For nuo it is 10 of
  52. A system-level score (union of the system's families) is the next step for these. Per-family
  absence understates surprise.
- **The strongest core-operon absences are a mix of biology, catalogue splits and binning loss:**
  - **UBA845 vamb_049** (97% complete) lacks nuo: 2.45 bits, P_intact 0.63. It has no nuo genes at all
    and carries the Na⁺-translocating Nqr instead, a real biological absence.
  - **Peredibacter vamb_623** (97%) lacks CF00033: 2.39 bits. It has ATP synthase α, β, γ and δ in
    CF02329, a 2-bin family. This is a catalogue split.
  - **Nanopelagicus lorbin_029** (81%, P_intact 0.69) lacks CF00000: 2.45 bits. None of the S10-spc core
    genes (uL2-uL6, uS3, uS10, uS19) are in the bin, though other ribosomal proteins are. This is the
    composition-loss case (`model_spec.md` λ_k: the rp-operon contig was not binned).
- **Strongest absences overall** (`tops.log`):
  - **Thioglobus_A lorbin_001** (99% complete) lacks the PntAB transhydrogenase CF00023: 2.89 bits.
  - **Pseudothioglobus sp913057955, vamb_711_sub** (CheckM2 88% complete, 21.6% contaminated) accounts
    for 5 of the top 15 absences: NusA/IF-2/RbfA, DnaA-dnaN-recF-gyrB, uL13-uS9, leu, YchF. Its
    completeness is overstated, and `pangenome_survey.md` already found a contaminant nuo module in this
    species.

## Most surprising states per bin

Full tables: `states.tsv` (every present state, plus absences > 1 bit, with Bakta medoid products and
p̂ under both priors) and `top_states_per_bin.tsv` (top 5 present and top 5 absent per bin).

- **Present.** The top of every bin is singleton families, 5-10 bits, mostly hypothetical proteins. Among
  families in ≥ 3 bins with an annotated product, the highest are in phylogenetically isolated bins
  (Peredibacter vamb_623, UBA845 vamb_049; 6.2-6.6 bits, p̂ ≈ 0.01). Examples:
  - Peredibacter: cytochrome bc1 b subunit, protoheme IX farnesyltransferase, UDP-GlcNAc oxidation;
  - UBA845: UbiB, MsrB;
  - Pseudothioglobus lorbin_018: CdsA / pseudouridine synthase.

  These are families shared with 2-3 distant bins. Level 1 can only call them "not this clade's usual
  content".
- **TMED14 maxbin_102.**
  - Absent: Lon/PBP/cytidylate kinase (CF00629), the histidine-utilization cluster (CF00179) and
    KdsB/RfbD (CF00409), each 2.8-3.0 bits at p̂ ≈ 0.93-0.95.
  - Present: TonB-dependent receptor and SLT transglycosylase families at 2.9 bits each.
- **UBA9145 maxbin_040.**
  - Absent: DnaX/RecR (CF00056) and an ABC/transaldolase cluster (CF00546), 2.5-3.0 bits.
  - Present: an RsmI/MraZ locus (CF04531, singleton), error-prone polymerase / PurC / AroC (CF01023) and
    a MacB efflux locus, 5.1-5.2 bits.

## Per bin

`per_bin.tsv`. Medians over the 93 bins:

| Measure | Median |
|---|---|
| Families present | 184 |
| Present families in no other bin (singletons) | 29 |
| Presence bits | 476 (374 per Mb; 229 per Mb excluding singleton families) |
| Absence bits | 241 (235 per expected Mb) |

- The top presence bins per Mb are mostly < 35% complete. With little binned DNA, a few singleton
  families dominate the total. Presence bits per Mb still correlate weakly with completeness (ρ = −0.26).
- The top absence bins are complete and contaminated:

  | Bin | CheckM2 compl./contam. |
  |---|---|
  | semibin_272_sub | 0.69 / 0.45 |
  | semibin_223 | 0.74 / 0.25 |
  | maxbin_069_sub | 0.71 / 0.22 |
  | maxbin_043 | 0.94 / 0.32 |
  | semibin_267_sub | 0.82 / 0.27 |
  | maxbin_115 | 0.79 / 0.08 |
  | maxbin_050 | 0.97 / 0.04 |
  | maxbin_102 | 0.84 / 0.18 |

## Circularity and what the numbers do not show

- **The families are defined from these same 93 bins.**
  - A locus exists only if it had Spacedust support from another genus.
  - Family edges were cut by similarity among these loci.
  - The family universe therefore includes only content seen in at least one bin.

  Singleton families exist only because the held-out bin carries them, so their presence surprise is
  partly a catalogue artefact: their partners matched below the cut. A family that every bin lacks
  cannot be scored at all.
- **Leave-one-out is not leave-one-population-out.** The catalogue's best-match and cross-genus rules
  saw the held-out bin. ρ_k, a0/b0 and κ are global (7 + 2 numbers), as are Leiden and MMseqs2 families.
  Species-level leakage is minimal: only UBA10364 has two named-species bins.
- **The held-out log loss compares priors fairly with each other, but not with an external truth.** A
  prior estimated from GTDB pangenomes (PanGBank) or from the marine and freshwater co-assemblies
  (overlapping populations: not independent, see `coverage_test.md`) is the needed out-of-sample test.
- **The detection model omits composition loss (λ_k) and depth (a_j).** Both would lower absence bits.
  lorbin_029's missing rp operon is the λ_k case.
- **Contamination enters ε through CheckM2 r_j.** Setting r_j = 0 changes presence totals by a median 3.7%,
  because ρ_k is small for rare families, so the contamination checks are not circular through ε. The
  absence-contamination link does go through CheckM2 completeness, which is itself inflated by
  contamination.
- **Presence surprise at level 1 does not separate contamination from genuine accessory content**
  (accessory-like contigs score like off-population ones). Coverage and anchoring (`coverage_test.md`,
  `pangenome_survey.md`) or level 2 (context) are needed for that.

## Outputs (`/data/scratch/metapanclust_data/catalogue_surprise/`)

- `inputs.npz`, `bins.tsv`, `loci_span.tsv`: observation matrix, spans, contig lengths, bin metadata.
- `fit_{flat,none,full}.npz`, `logloss_*.tsv`, `kappa_*.json`: fits per training-count mode. `flat` is
  the main fit; `full` is the leaky contrast.
- `per_bin.tsv`, `states.tsv`, `top_states_per_bin.tsv`, `core_operons.tsv`, `contig_class.tsv`,
  `calibration.tsv`, `*.log`.

## Next

1. Fix the leave-one-out in `model_spec.md` §3.4 A: per-bin refit, or leave-two-out posteriors.
2. Score systems (union of a system's families) as well as families. Match truncated loci by
   containment, so that split operons stop looking absent.
3. Use an external prior: GTDB species representatives through the same catalogue, then score
   Chesterfield bins against it.
4. Add absence bits per expected Mb as a contamination flag next to CheckM2, and test it on bins with
   known chimeras (TMED14, vamb_711_sub).
