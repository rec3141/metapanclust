# Background: pathway-level benchmark (November 2026)

These experiments started from benchmarking danaSeq's open MetaCyc pathway scorer against PathoLogic
(Pathway Tools, via MetaPathways 4.0.3), with each organism's family (order for *E. coli*) removed from
both tools' reference sequences. Truth is EcoCyc for *E. coli* K-12 and, for the other organisms,
MetaCyc pathways with experimental evidence in that organism (recall only: these sets are incomplete).

| Organism (truth size) | PathoLogic calls / recall | Scorer all tiers | Scorer high tier |
|---|---|---|---|
| *E. coli* K-12 (346, EcoCyc) | 474 / 0.98 (precision 0.71) | 745 / 0.90 (precision 0.42) | 242 / 0.55 (precision 0.78) |
| *B. subtilis* 168 (121) | 321 / 0.84 | 584 / 0.83 | 169 / 0.40 |
| *P. aeruginosa* PAO1 (94) | 417 / 0.59 | 761 / 0.77 | 214 / 0.31 |
| *M. tuberculosis* H37Rv (78) | 253 / 0.91 | 603 / 0.81 | 140 / 0.27 |
| *Synechocystis* PCC 6803 (44) | 240 / 1.00 | 501 / 0.84 | 119 / 0.27 |
| *M. jannaschii* (46) | 123 / 0.78 | 222 / 0.89 | 51 / 0.39 |

Findings that motivated this repo:

- Recall is roughly tied, but the scorer needs about twice as many calls. Its calls that PathoLogic does
  not make are almost all wrong (3% in *E. coli* against EcoCyc).
- The wrong extra calls come from complete pathways outside the organism's MetaCyc taxonomic range
  (about a quarter), weak evidence classes (EC numbers without sequence support), close variants of
  real pathways, and small or incomplete pathways. MetaCyc taxonomic ranges are a coarse, curated
  yes/no filter: a learned, clade-resolved prior should do better.
- Gene clustering separated correct from wrong calls only weakly when pooled (*E. coli*: share of
  supporting genes within 3 genes of each other, 0.33 vs 0.26), because clustering only carries signal
  for pathways that are normally clustered. Learning per-unit clustering expectations needs many
  genomes, which led to taking gene clusters themselves as the unit.
