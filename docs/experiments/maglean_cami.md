# MAGLean on CAMI II marine (2026-11-24)

[MAGLean](https://github.com/rec3141/maglean) (Lean 4 formalisation and proof-carrying QC auditor; native `maglean`
CLI, commit 0871d7a, built with Lean v4.36.0) run on 300 CAMI II marine bins: 20 sampled from each of the 15
submitted binnings of the pooled gold-standard assembly. Scripts: `experiments/cami_marine/maglean/`
(`build_evidence.py`, `compare.py`); full output: [maglean_cami_results.txt](maglean_cami_results.txt).

## Evidence given to MAGLean

| field | source |
|---|---|
| marker sets | CheckM v1 domain-level collocated sets (Bacteria 58 sets / 104 markers; Archaea 108 / 150), domain by Kraken2 bp-majority |
| contig markers | Pyrodigal (meta) + pyhmmer against those 207 CheckM HMMs at their trusted cutoffs |
| contig taxonomy | Kraken2, standard 16 GB DB (2025), lineage domain..species from the DB's own taxonomy |
| RNA | not available: all false, so the MIMAG tier cannot reach HQ |

The Kraken2 database contains RefSeq relatives of most CAMI marine genomes (601 of 777 are "known strains"), so
taxonomic purity is optimistic here compared with novel environmental MAGs.

## Results (300 bins; 139 with ≥ 5% true contamination)

**Completeness and contamination estimates:** CheckM2 is better on both. Completeness error is 7.0 points for
CheckM2 against 9.5 for MAGLean's domain-level markers (Spearman 0.89 vs 0.85). Contamination error is 14.3 against
20.8 points (Spearman 0.64 vs 0.46). MAGLean is limited by using domain-level sets rather than lineage-specific ones.

**Flagging contaminated bins (AUC):**

| score | true contamination ≥ 5% | ≥ 10% |
|---|---|---|
| CheckM2 contamination | 0.842 | 0.871 |
| MAGLean marker contamination | 0.749 | 0.784 |
| MAGLean conservative impurity, genus rank | 0.804 | 0.778 |
| MAGLean conservative impurity, species rank | **0.863** | 0.856 |
| MAGLean max chimerism over ranks | 0.848 | 0.842 |
| **CheckM2 + species impurity** (logistic, CV grouped by genome) | **0.904** | **0.922** |

Taxonomic impurity at species rank is as good as CheckM2 on its own, and the two are complementary.
metapanclust's novelty signal added only +0.007 to CheckM2 ([cami_marine.md](cami_marine.md)); this adds +0.06.

**Alerts:**
- The `Contaminated` tier (marker contamination ≥ 10%) catches 71 of 116 bins with ≥ 10% true contamination, and
  only 1 of 106 clean bins.
- `STEALTH_CHIMERA_ALERT` fires on 22% of bins with ≥ 5% true contamination but also on 14% of clean bins. The
  false alarms come from contigs Kraken2 leaves unclassified or misplaces: their median unclassified fraction at
  species rank is 8%, against 0% for clean bins without the alert. Conservative purity counts unclassified bp
  against the bin.
- `UNCLASSIFIED_DILUTION_ALERT` fires on 13-52% of bins in every class, so on these data it is not discriminating.

## Notes for MAGLean

- The stealth-chimera alert would gain from a minimum classified fraction, or from using the
  `[purityCons, purityOpt]` interval rather than `purityCons` alone.
- Lineage-specific marker sets (CheckM's lineage workflow, or GTDB-Tk placement) would close most of the gap to
  CheckM2 on completeness and contamination.
- An RNA field from barrnap and tRNAscan-SE (both in the danaSeq container) would make the HQ tier reachable.

## Rerun with PR #1 (`ea0c4f9`) and issues filed

PR #1 changes `STEALTH_CHIMERA_ALERT` to optimistic purity. Rerun on the same 300 bins, the alert still fires on 13.2%
of clean bins (14.2% before). 9 of the 14 clean-bin alerts fail at species rank only, where Kraken2 splits a CAMI
new species or new strain across congeners. AUCs and tiers are unchanged.

The CLI's ≥95% Hoeffding interval contains the true bp-completeness for 88% of bins: 92% of clean bins and 80% of
bins with ≥ 10% contamination.

Filed: maglean #8 (stealth-chimera false positives) and #9 (dilution alert not discriminating). Commented with
real-data numbers on #4 (interval coverage), #5 (domain-level marker family against CheckM2) and #6 (first real-MAG
run; converter location).

## Rerun with maglean `0e5ca7c` (fixes for #8 and #9)

What changed: the stealth-chimera alert is evaluated at genus rank and above and needs ≥ 50% classified bp; the
dilution alert fires only when an impure coarser rank (< 90%) is masked at a finer rank; species impurity is exported
as `species_conservative_impurity_exact`. Rerun on the same 300 bins
([maglean_cami_results_0e5ca7c.txt](maglean_cami_results_0e5ca7c.txt)):

| true contamination | bins | stealth alert (before → after) | dilution alert (before → after) |
|---|---|---|---|
| < 1% | 106 | 13.2% → **1.9%** | 13.2% → **0%** |
| 1–5% | 55 | 7.3% → 7.3% | 45.5% → 0% |
| 5–10% | 23 | 17.4% → 8.7% | 52.2% → 0% |
| ≥ 10% | 116 | 21.6% → 18.1% | 37.9% → 3.4% |

- **Stealth alert:** precision is now 79% (23 of 29 alerts on bins with ≥ 5% true contamination).
- **Dilution alert:** now specific, but it rarely fires on these data.
- **Unchanged:** AUCs, estimates and tiers. The exported species impurity reproduces AUC 0.863.
- **New contig-based interval:** median radius ±30.9%; it contains the true completeness for 93% of bins (99% of
  clean bins, 85% of bins with ≥ 10% contamination). The marker-set interval stays at 88%.
- **Provenance:** the binary reported `git_commit` `0b946d1` while built from `0e5ca7c`'s tree. Filed as maglean #10.
