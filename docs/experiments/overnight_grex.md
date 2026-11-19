# Overnight grex runs (2026-11-18/19)

Raw outputs: workstation `/data/scratch/metapanclust_data/grex_overnight/`; grex `~/scratch/metapanclust/{chf_3di,pangbank_scale}/`.

## Chesterfield, Spacedust structure mode (job 8187201, one L40S)

ProstT5 3Di for 173,368 binned proteins in 4,759 s (36 proteins/s); all-vs-all clustersearch over 93 bins in 1,623 s.

| | MMseqs2 mode | 3Di mode |
|---|---|---|
| pairwise conserved clusters | 264,176 | 378,224 (+43%) |
| genes in >= 1 cluster | 101,134 (58.3%) | 105,873 (61.1%) |

Much smaller coverage gain than on the six distant benchmark genomes (18.1% -> 26.5%): within one co-assembly
most bins have a relative close enough for sequence search, so structure search mostly adds more cluster pairs
for genes already covered. Not yet fed into the catalogue.

## PanGBank scale (job 8187309, 2 genoa nodes; 455 pangenomes, 291 species / 52 phyla for A, 40 genera / 203 focal species for B)

Full tables: `pangbank_scale_results.txt`.

**A, context model vs genomes per species** (median % of prevalence log-loss removed, logistic C = 0.1):

| N genomes | 15 | 25 | 50 | 100 | 200 | 400 | 800 |
|---|---|---|---|---|---|---|---|
| accessory families | 20.0 | 26.6 | 36.9 | 51.9 | 63.4 | 69.8 | 78.9 |
| modules | 24.7 | 33.7 | 42.5 | 52.0 | 62.0 | 68.1 | 75.9 |
| kNN-8 (clade proxy) | 8.0 | 10.6 | 15.7 | 22.4 | 28.3 | 33.3 | 38.7 |

The gain grows by about 9-10 points per doubling of N (median per species) and has not levelled off at 800
genomes; it is mostly context, not clade structure (kNN reaches half as much). Best regularisation moves from
C = 0.03-0.1 at small N to 0.3 at 800.

**B, genus shrinkage vs species-only empirical Bayes** (held-out bits per genome, m training genomes):

- All clusters: shrinkage loses at every m (median -10%, -9%, -7% at m = 5, 10, 20; wins for 12%, 5%, 4% of focal species).
- Only on clusters shared with congeners does it help, and only at m = 5 (+2.4%, wins 64%); by m = 10-20 it is neutral.
- Fitted genus kappa is ~1 (IQR 1-2).

This confirms phylo_prior.md and phylo_prior_check.md at scale: once a species has a handful of genomes, its own
data screen off the genus; tree information belongs in the prior only for species with ~1-2 genomes or for shared
clusters, ideally via the ancestral-state predictor.
