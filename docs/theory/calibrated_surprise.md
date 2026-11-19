# Calibrated surprise: statistical and information-theoretic footing (draft, 2026-11-18)

**What this is.** A theory note for the surprise score of `model_spec.md` §4 and the binning rule of §7. It
answers four questions:

1. **Calibration.** What should a bin's total surprise be compared with, so that the number means the same
   thing across species and across reference sets of different size n?
2. **Dependence.** How much do co-occurring clusters overcount, and what replaces the independent sum?
3. **Multiple testing.** How to call "surprisingly present/absent" over ~10^5-10^6 bin × cluster cells.
4. **MDL.** Binning as coding: where the merge/move rule comes from, what the right normaliser is, and why
   absence surprise lines up with contamination.

Equations are LaTeX (`$…$`). Logs are base 2 unless written ln. A light numerical check on the
`catalogue_surprise` outputs is in §5; its script and logs are in
`/data/scratch/metapanclust_data/theory_surprise/` (`check.py`, `check2.py`, `check3.py`, `*.log`,
`per_bin_calibrated.tsv`; `check4.py` has the final E± definitions).

## Summary

- **Excess surprise per cell is a weighted residual.** For a two-state cell with predicted P(seen) = q and
  outcome y ∈ {0,1}, the excess over the model's expectation is
  $e = s - h = (y-q)\,\log_2\frac{1-q}{q}$. It has mean 0 and variance $q(1-q)\log_2^2\frac{1-q}{q}$ under the
  model. Bin totals of e are then centred and can be standardised (z-score) or given an exact saddlepoint
  tail. This is the comparable-across-species quantity, not raw bits.
- **Direction belongs to the prediction, not the outcome.** Split cells by whether q < 1/2 (a presence
  would be surprising) or q > 1/2 (an absence would be). Splitting by the *observed* state, as in
  `model_spec.md` §4.2, gives sums that have no zero mean: every observed absence of a rare cluster adds a
  small negative excess, and thousands of them dominate the absence total.
- **Bits per cell are already log p-values; n enters correctly through the posterior predictive.** A
  module absent from all n genomes costs $\log_2(n+2)$ bits under Laplace's Beta(1,1) prior, and
  $2^{-s} = 1/(n+2)$ *is* its p-value. Bits are comparable as evidence; what is not comparable is the bin
  *total*, because the number of cells that can be surprising, and the expected total, both grow with n.
  Replace Beta(1,1) by the fitted hierarchical prior: $s = \log_2\frac{n+\kappa}{\kappa\,\hat p_{\mathrm{pa}}}$.
  That is close to $-\log_2 \hat p_{\mathrm{pa}}$ while n ≪ κ, and grows as $\log_2 n$ only once n ≫ κ.
- **Dependence.** The chain rule turns the independent sum into a joint code length:
  $S_{\text{joint}} = S_{\text{indep}} - \sum_{\text{edges}} \mathrm{pmi}$ on a Chow–Liu tree. For the 36
  UBA9145 modules that always co-occur this gives 5.75 bits (one event) instead of 207. The level-2
  pseudo-likelihood (summing full conditionals) errs the other way: it scores that block near 0 bits.
  For z-scores the fix is the variance: $\mathrm{Var}(E) = w^\top \Sigma w$, and the design effect
  $\mathrm{VIF} = w^\top\Sigma w / \sum_k w_k^2\sigma_k^2$ deflates both the excess and the z.
- **Multiple testing.** The null distribution of every cell's surprise is known exactly from the model, so
  the expected number of null calls above any threshold is computable (Efron-style FDR without fitting a
  null). Raw surprise also gives an e-value for free: $2^{s}/|\mathcal O|$ has expectation 1 under the
  model. Conditional e-values multiply along the chain rule, so the joint code length gives a valid
  e-value under dependence. Use e-BH, or local fdr with the exact null.
- **MDL.** With mixture (Bayes) codes, total description length is $-\log_2$ of the marginal likelihood,
  so ΔL between two binnings is exactly $\log_2$ of the Bayes factor in `model_spec.md` §7 (plus the
  partition prior). The bits-per-Mb merge rule of `binning_sim.md` is an approximation: it codes the
  candidate contig's alternative at the bin's own average rate. The right normaliser is the code length of
  the contig under the best *alternative* home (another bin, or the domain-level prior), or, for a
  calibrated score, the model's expected cost of that contig (excess, not raw rate).
- **Check on the data (§5):**
  - Absence bits per expected Mb separate ≥ 10%-contaminated bins with AUC 0.845, but the model's
    *expected* absence bits alone give 0.80. With contamination removed from ε (r_j = 0), the observed
    absence score falls to AUC 0.54. Most of the published absence signal is therefore CheckM2's r_j fed
    back through the false-presence term, at $\approx \varepsilon/\ln 2$ bits on each of ~7,000 absent rare
    families. It is circular, not a discovery. The calibrated absence z-score keeps a weaker, honest
    signal (AUC 0.61-0.71).
  - The calibrated **presence-direction** excess with r_j = 0 (no CheckM2 input) separates the same
    contaminated bins with **AUC 0.96** (Z+; 0.84 when r_j is left in ε, because ε then explains the
    foreign content away). The contamination signal is in unexpected presence once it is measured against
    its expectation, not in absence.
  - Singleton families: 4,121 present cells observed against 4,433 predicted. In aggregate the model
    already expects them, so their 48% share of presence bits is expected surprise, not excess.
  - Cells whose observed binary state had probability < 0.01: 4,175 observed against 1,923 expected under
    the model. A 1% per-cell cut therefore has an FDR near 46%.

---

## 1. Calibration

### 1.1 Surprise, entropy and excess for one cell

For bin j and cluster k the model gives a predictive distribution $q_{jk}(o)$ over observed states
$o \in \mathcal O = \{I, P, A\}$ (`model_spec.md` (D1)). The surprise of the observed state $O_{jk}$ is

$$ s_{jk} = -\log_2 q_{jk}(O_{jk}). $$

If the model is calibrated, i.e. $O_{jk} \sim q_{jk}$, then

$$ \mathbb E[s_{jk}] = H(q_{jk}) = -\sum_o q_{jk}(o)\log_2 q_{jk}(o) \equiv h_{jk}, \qquad
   \mathrm{Var}[s_{jk}] = \sum_o q_{jk}(o)\log_2^2 q_{jk}(o) - h_{jk}^2 \equiv v_{jk}. $$

If instead the truth is $p^*$, $\mathbb E_{p^*}[s] = H(p^*) + D_{\mathrm{KL}}(p^*\|q)$: the cross-entropy. The
excess

$$ e_{jk} = s_{jk} - h_{jk} = \sum_o \big(\mathbf 1[O_{jk}=o] - q_{jk}(o)\big)\,\big(-\log_2 q_{jk}(o)\big) $$

is linear in the outcome indicators, with mean $\sum_o (p^*(o)-q(o))(-\log_2 q(o))$. That is 0 under the
model, and positive when the truth puts more mass on states the model thought unlikely.

**Two-state form.** Collapse to seen/not seen, with $q$ = P(seen) and $y$ = 1 if seen:

$$ e = (y - q)\,\log_2\frac{1-q}{q}, \qquad v = q(1-q)\,\log_2^2\frac{1-q}{q}. $$

So the excess is the Pearson residual $y-q$ weighted by the logit of the prediction in bits. Two points
follow:

- **Sign is set by the prediction.** A cell with q < 1/2 has positive weight: a presence scores +, an
  absence scores a little −. A cell with q > 1/2 has negative weight: an absence scores +.
- **Most cells contribute small negative excess.** An absent cluster with q = 0.01 scores s = 0.0145 and
  h = 0.081 bits, so e = −0.066. Summed over thousands of absent rare families, these terms swamp
  everything else in a total split by observed state. This is why §5 finds that the "absence excess bits"
  of `model_spec.md` §4.2 anti-correlate with contamination (ρ = −0.69): contaminated bins have larger ε,
  hence larger h on every rare cell.

**Directional sums (replacing §4.2's S+ and S−).**

$$ E_j^{+} = \sum_{k:\,q_{jk}<1/2} e_{jk}, \qquad E_j^{-} = \sum_{k:\,q_{jk}>1/2} e_{jk}, \qquad
   V_j^{\pm} = \sum_{k \in \pm} v_{jk}. $$

Each has mean 0 under the model. $E^+$ is "unexpected content", $E^-$ "missing expected content". For three
states, use the same split on $q(\text{seen}) = q(I) + q(P)$, and keep the three-state $e_{jk}$ in the
sum.

### 1.2 Standardisation and tail probabilities per bin

Under independence (§2 removes this),

$$ Z_j^{\pm} = E_j^{\pm} \big/ \sqrt{V_j^{\pm}} . $$

The normal approximation is poor for $E^+$. It is a sum of a few highly skewed terms: a q = 0.005 cell
contributes $(1-q)\cdot 7.6$ bits with probability 0.005, and −0.04 bits otherwise. Use the exact
cumulant generating function of a sum of independent discrete variables:

$$ K_j(t) = \sum_k \ln \sum_o q_{jk}(o)\, e^{t\,e_{jk}(o)}, $$

where $e_{jk}(o)$ is the excess that state o would give. The Lugannani–Rice saddlepoint tail for the
observed total $E$ is

$$ \Pr(E_j \ge E) \approx 1 - \Phi(\hat w) + \phi(\hat w)\Big(\frac{1}{\hat u} - \frac{1}{\hat w}\Big), \quad
   K_j'(\hat t) = E,\; \hat w = \mathrm{sgn}(\hat t)\sqrt{2(\hat t E - K_j(\hat t))},\; \hat u = \hat t\sqrt{K_j''(\hat t)}. $$

It is one 1-D root-find per bin, O(K) per Newton step. Restricting to the K' cells with
min(q, 1−q) > 1e-6 loses nothing measurable. Report $-\log_2$ of this tail as **calibrated bin surprise**
(bits of evidence against "this bin is a typical draw from its clade"). It is on the same scale for any
species and any n.

**Per cluster.** A single cell needs no approximation. The exact (mid-)p-value of observing a state at
least as unlikely is

$$ \pi_{jk} = \sum_{o:\,q(o) < q(O)} q(o) + \tfrac12 \sum_{o:\,q(o) = q(O)} q(o). $$

For a two-state cell whose observed state is the minority one, $\pi = q(O)$ (or q/2 mid-p), i.e.
$\pi = 2^{-s}$. **Per-cluster surprise in bits already is a log p-value**, and needs no recalibration
beyond the model being right.

### 1.3 Where n enters: posterior predictive and the prior

Level 1 with species data only and prior Beta(a, b) on θ, x carriers among n genomes:

$$ q = \Pr(z=1 \mid x, n) = \int \theta\, \mathrm{Beta}(\theta;\, a+x,\, b+n-x)\, d\theta = \frac{x+a}{n+a+b}. $$

This is the posterior predictive. It integrates over θ, so the estimation uncertainty in θ is already in q,
and $-\log_2 q$ is the correct prequential code length for the next genome (Dawid 1984). With
$a = b = 1$ (the prototype), x = 0 gives $s_{\text{present}} = \log_2(n+2)$: 5.75 bits for UBA9145 (n = 52),
3.81 for n = 12, 2.32 for n = 3 (`pangenome_survey.md`).

**What is right and what is not about $\log_2(n+2)$:**

- **Right:** with fewer references, a novel module really is less surprising. Under the model the chance
  that a fresh genome carries a module unseen in 3 genomes is 1/5, against 1/54 after 52. Bits are
  comparable *as evidence*: 2.32 bits for n = 3 is p = 0.2, which no one would flag.
- **Wrong 1: the prior.** Beta(1,1) puts a uniform prior on prevalence. Pangenome prevalences are
  U-shaped (core near 1, cloud near 0; `catalogue_surprise.md` finds a0, b0 < 1 at the root). The prior
  pseudo-count should come from the clade above, as in `model_spec.md` §3:

  $$ q = \frac{x + \kappa\,\hat p_{\mathrm{pa}}}{n + \kappa}, \qquad
     s_{\text{present}}(x{=}0) = \log_2\frac{n+\kappa}{\kappa\,\hat p_{\mathrm{pa}}}
     \;\approx\; \begin{cases} -\log_2 \hat p_{\mathrm{pa}} & n \ll \kappa \\
     \log_2 n - \log_2(\kappa \hat p_{\mathrm{pa}}) & n \gg \kappa \end{cases} $$

  With κ_genus = 30 (fitted in `catalogue_surprise.md`), a species with n = 3 is scored almost entirely by
  its genus. Only beyond n ≈ 30 does the species' own absence start to add $\log_2 n$. A module that is
  also rare in the genus ($\hat p_{\mathrm{pa}} = 0.01$) costs 6.7 bits at n = 3, not 2.3. A module common
  in the genus ($\hat p_{\mathrm{pa}} = 0.5$) costs 1.1 bits at n = 3 and 2.5 at n = 52, not 5.75. The "+2" of
  Laplace should be "+κ", and the "1" in the numerator should be "κ p̂_pa".
- **Wrong 2: totals.** The expected total $H_j = \sum_k h_{jk}$ also depends on n. With small n every q
  is pulled towards the parent, h per cell is larger, and the maximum attainable s is lower. A raw bit
  total therefore mixes the reference set's size with the bin's atypicality. The excess $E_j$ and its
  tail (§1.2) remove both effects to first order.
- **Effective n.** Genomes enter with weights $w_g = m_g^{-\gamma}$ (`model_spec.md` §3.2). For the
  posterior predictive, n in the formula is $\sum_g w_g$. When weights are uneven, the effective sample
  size for the Beta update is closer to Kish's $(\sum w_g)^2/\sum w_g^2$. Use the smaller of the two when
  the reference set is dominated by one study.

**Single-cell numbers (two-state, presence of a module absent from all n references, Laplace prior):**

| n | q = 1/(n+2) | s (bits) | h (bits) | e = s − h | z = e/√v = √((1−q)/q) | mid-p |
|---|---|---|---|---|---|---|
| 3 | 0.200 | 2.32 | 0.72 | 1.60 | 2.0 | 0.10 |
| 12 | 0.071 | 3.81 | 0.37 | 3.44 | 3.6 | 0.036 |
| 52 | 0.0185 | 5.75 | 0.13 | 5.62 | 7.3 | 0.0093 |

The single-cell z grows as $\sqrt n$ while the p-value falls as 1/n: for one cell, z is the wrong scale
(the variable is far from normal), and the p-value (equivalently s) is the right one. For bin totals,
where many cells add up, the saddlepoint tail is the right one.

**How surprise should behave as n grows.** Two different quantities, both useful:

- **Predictive surprise** $s = -\log_2 q$ (what we use). For a cluster the clade really lacks
  (θ = 0), it grows without bound, as $\log_2 n$. For a cluster at true prevalence θ > 0 it converges to
  $-\log_2\theta$. It does not shrink: more references make an atypical genome *more* surprising, which is
  correct.
- **Bayesian surprise** $D_{\mathrm{KL}}(\text{post after } j \,\|\, \text{post before } j)$ about θ
  (Itti & Baldi 2009): how much bin j should change the species model. For a first presence after x = 0 of
  n (Beta(1, n+1) → Beta(2, n+1)):

  $$ D_{\mathrm{KL}} = \ln(n+2) + \psi(2) - \psi(n+3) \;\approx\; \psi(2) + \ln\frac{n+2}{n+2.5} \;\to\; 1-\gamma = 0.42\text{ nats} = 0.61\text{ bits}, $$

  independent of n (the posterior mean doubles on a log scale whatever n is). For a bin that agrees with
  the majority it falls as O(1/n). This is the quantity for "should this bin enter the reference set",
  not for flagging content.

### 1.4 Calibration diagnostics to run routinely

- **Reliability table** of predicted P(seen) against observed rate. `catalogue_surprise/calibration.tsv`
  shows the shrunk prior under-predicts presence at the bottom (predicted 0.0053, observed 0.0115 in the
  0-0.01 bin) and in the middle (0.31 → 0.41). The bottom bin holds 364k cells, so excess presence
  surprise is inflated systematically. Recalibrate q with isotonic regression or a two-parameter logistic
  map, $\mathrm{logit}\,q' = \alpha + \beta\,\mathrm{logit}\,q$, fitted on held-out cells, before computing e.
- **PIT / null-count check:** the number of cells with π < t observed, against $N_0(t)$ of §3.1.
- **E[Z] ≈ 0 and Var[Z] ≈ 1 on clean bins.** Model_spec §9 E4 (clean complete genomes cut into MAG-like
  contigs). Over the 93 Chesterfield bins (not all clean), the median Z is 0.06 and the SD 5.7: the
  centring works, the spread does not. That is the dependence of §2 (VIF ≈ 30 if all of it is
  dependence) together with the miscalibration above.

---

## 2. Dependence

### 2.1 How independent summation overcounts

The total $S_j = \sum_k s_{jk}$ is the code length of the bin's content under a code that treats clusters
as independent. If the true joint is $P(O_{j1},\dots,O_{jK})$, the excess code length of the independent
code over the joint code is, in expectation, the **total correlation**:

$$ \mathbb E\Big[\sum_k s_k - s_{\text{joint}}\Big] = \sum_k H(O_k) - H(O_1,\dots,O_K) = C(O_1;\dots;O_K) \ge 0. $$

**UBA9145.** 36 flagellar/chemotaxis modules, each in 0 of 52 public genomes, all present in the bin.
Independent sum: 36 × 5.75 = 207 bits. If the 36 are one locus (or one co-binned genome), the event is
"the locus is present". With a block indicator B and $\Pr(B=1) = 1/54$ by the same Laplace argument,

$$ s_{\text{joint}} = -\log_2\Pr(B=1) - \sum_{m=1}^{36}\log_2\Pr(\text{module } m \mid B=1) \approx 5.75 + 36\,\log_2\frac{1}{\rho}, $$

with ρ the probability that a module is seen given that its locus is present. ρ = 1 gives 5.75 bits;
ρ = 0.95 (detection loss) gives 8.4 bits. Between 6 and 9 bits, against 207.

### 2.2 Chain rule: surprise under a joint model

Any ordering $k_1, k_2, \dots$ gives an exact decomposition

$$ s_{\text{joint}}(O) = \sum_{i} -\log_2 \Pr(O_{k_i} \mid O_{k_1},\dots,O_{k_{i-1}},\text{clade}). $$

Each term is a conditional surprise: the "new information" in cluster $k_i$ given what has been scored
already. This is a proper code length and a proper log-likelihood, unlike the level-2 sum.

**Why the level-2 pseudo-likelihood is wrong in the opposite direction.** `model_spec.md` §4.2 option 2
sums full conditionals $-\log_2\Pr(O_k \mid O_{-k})$. For a block of 36 perfectly linked modules, each one
is predicted with certainty by the other 35, so the block scores ≈ 0 bits. The pseudo-likelihood
undercounts linked blocks exactly as much as the independent sum overcounts them. Only an ordered
(autoregressive) conditional is a joint code.

**Tractable joint: Chow–Liu tree.** Approximate the joint residual dependence by a maximum-weight spanning
tree T on the clusters, with edge weights the mutual information of their states in the reference set
(after removing the clade effect; §2.4). The tree code length is

$$ s_T(O) = \sum_k s_k(O_k) - \sum_{(k,l) \in T} \mathrm{pmi}(O_k; O_l), \qquad
   \mathrm{pmi}(o_k; o_l) = \log_2\frac{\Pr(o_k, o_l)}{\Pr(o_k)\Pr(o_l)}. $$

For the 36-module block on a chain with every pair present together whenever either is,
$\mathrm{pmi} = \log_2(1/q) = 5.75$ on each of 35 edges, and $s_T = 207 - 35 \times 5.75 = 5.75$ bits.
pmi is estimated from the pairwise 2 × 2 tables with the same shrinkage as the marginals (a pair never seen
together in the species borrows its association from the genus or family, where the locus does occur).
The tree is fitted once per clade. Scoring is O(K).

**When the conditionals have no data.** A module never present in the species has no within-species
co-occurrence. The association must then come from:

1. the catalogue (physical adjacency in other clades: same locus);
2. higher clades where both occur;
3. as a fallback, contig co-location in the bin itself. Clusters on one contig whose loci are adjacent
   within 5 kb are one event unless the reference set shows them separating.

This is the `model_spec.md` §5.1 rule (merge same-locus clusters in the catalogue), turned into a scoring
rule.

### 2.3 Effective number of independent events (z-scores)

When only the standardised total is needed, the variance carries the dependence. With
$e = (e_1,\dots,e_K)$ the cell excesses, $\sigma_k = \sqrt{v_k}$ and R the correlation matrix of the
residuals,

$$ \mathrm{Var}(E_j) = \sum_{k,l} \sigma_k\sigma_l R_{kl}, \qquad
   Z_j^{\text{dep}} = \frac{E_j}{\sqrt{\sigma^\top R\,\sigma}} = \frac{Z_j}{\sqrt{\mathrm{VIF}_j}}, \qquad
   \mathrm{VIF}_j = \frac{\sigma^\top R\,\sigma}{\sigma^\top\sigma}. $$

The effective number of independent events (Kish/Bartlett form) is

$$ K_{\text{eff}} = \frac{\big(\sum_k v_k\big)^2}{\sum_{k,l} v_k v_l R_{kl}^2}. $$

For a block of m identical cells, VIF = m and $K_{\text{eff}}$ = 1. A deflated excess $E_j/\mathrm{VIF}_j$
reproduces "one event, ~5.6 bits" for UBA9145.

**Estimator.**

1. **Residuals.** In the reference genomes of the clade (leave-one-out), $r_{gk} = (\hat z_{gk} - \hat p_{k})
   / \sqrt{\hat p_k(1-\hat p_k)}$, with $\hat z$ the detection-corrected soft labels of §3.4.
2. **Sparsity.** Estimate $R_{kl}$ only for pairs that (a) share a locus in the catalogue or sit within
   5 kb in some genome, or (b) are among the top-M residual correlates of `model_spec.md` §5.1. Set all
   others to 0.
3. **Shrinkage.** Ledoit–Wolf towards 0 (n genomes is often < K). Phylogenetic correction: residuals
   within species after the clade offset, so shared ancestry is not counted as linkage (the Goldfinder
   point).
4. **In the bin.** Pairs present on the same contig in the bin, with no reference data, get $R_{kl} = 1$
   if adjacent (one locus) and 0 otherwise.

`check2.py` (§5) does step 1 on the ten families with mean q ≥ 0.6 across the 93 bins: $K_{\text{eff}}$ =
8.8 of 10. Core operons behave nearly independently once completeness and fragmentation are in the
model. The dependence that matters is among rare accessory modules, the UBA9145 case, where the
reference data are thinnest; there the catalogue adjacency rule (step 4) does most of the work.

### 2.4 Dependence that is not linkage

Two sources of correlated cells are not biology:

- **Bin-level nuisance shared by all cells**: completeness c, contamination r, fragmentation f. Errors in
  c shift every absence cell of a bin together, which is a rank-1 component of R ($R_{kl} = \lambda_k
  \lambda_l$). It explains part of the observed SD of Z (5.7 instead of 1). The fix is to integrate c over
  its uncertainty (CheckM2 error, about ±5%) inside q, or to estimate c from the content itself (§4.5),
  not to deflate.
- **Shared ancestry.** It is handled by the clade offset (level 1) and the within-species residuals above.

---

## 3. Multiple testing across bins × clusters

The catalogue run has 667,554 cells. Pangenome surveys at GTDB scale reach 10^7-10^8.

### 3.1 The null is known: exact expected false calls

Under the model, cell (j,k) takes surprise $s_{jk}(o)$ with probability $q_{jk}(o)$. The expected number of
cells at or above threshold t when every cell is null is

$$ N_0(t) = \sum_{j,k}\ \sum_{o:\, s_{jk}(o) \ge t} q_{jk}(o), $$

and an empirical-Bayes FDR estimate for the call set {s ≥ t} is

$$ \widehat{\mathrm{FDR}}(t) = \hat\pi_0\,\frac{N_0(t)}{\#\{(j,k): s_{jk} \ge t\}}, \qquad \hat\pi_0 \le 1 \text{ (take 1 to be conservative)}. $$

No null needs to be fitted: the model is the null. For the 93 bins (§5), the threshold "observed binary
state had probability < 0.01" gives 4,175 calls against $N_0$ = 1,923, an estimated FDR of 0.46. Half of
the strongly surprising cells are what the model itself expects to see once in a hundred times over
667,554 tries.

**Local fdr** for ranking: bin the s values and compute

$$ \mathrm{lfdr}(t) = \hat\pi_0\,\frac{f_0(t)}{f(t)}, \qquad f_0(t) \propto \frac{d}{dt}\big(-N_0(t)\big), $$

with f the observed histogram (Efron 2004, two-group model, theoretical null). If the model is
miscalibrated (§1.4), $N_0$ is off too. Recalibrate q first, or use Efron's empirical null on the bulk of
cells.

### 3.2 Testable cells (discreteness)

A cell whose minority state has probability $m_{jk} = \min_o q_{jk}(o)$ can never reach p < m. Tarone's
procedure: test only cells with $m_{jk} < \alpha/K_{\text{test}}$, choosing $K_{\text{test}}$ self-
consistently. In practice:

- **Absence calls** come only from cells with q(seen) high enough that $1 - q$ is small, i.e. clusters
  the clade expects *and* that detection could have seen. With completeness 0.7 and P_lost ≈ 0.3, the
  smallest attainable P(A) is about 0.3, so no single absence in such a bin can be significant. This is the
  2.0-2.9-bit ceiling of `catalogue_surprise.md`, and it is correct: one missing operon in a 70% bin is
  not evidence. Absence evidence must be aggregated (§1.2) or come from complete bins.
- **Presence calls** come from cells with small q. Every rare cluster is testable for presence, which is
  why presence dominates the call counts.

### 3.3 E-values

**Surprise gives an e-value directly.** For any cell,

$$ \mathbb E_q\big[2^{s}\big] = \sum_o q(o)\,\frac{1}{q(o)} = |\mathcal O|, \qquad\text{so}\qquad
   \varepsilon_{jk} = 2^{s_{jk}}/|\mathcal O| $$

has expectation 1 under the model: a valid e-value for "cell (j,k) follows the model". It is weak (it
spreads power over every alternative). The better e-value is a likelihood ratio against a named
alternative, which also has expectation 1 under the null:

$$ \varepsilon_{jk}^{\text{LR}} = \frac{\Pr_{\text{alt}}(O_{jk})}{q_{jk}(O_{jk})}. $$

Useful alternatives:

- **Foreign presence:** state drawn from the domain-level or donor-clade prior.
- **Contamination at fraction π:** ε raised by π G ρ_k.
- **Detection failure:** c lowered.

**Combining.**

- **Within a bin, under dependence:** the product of *chain-rule conditional* e-values,
  $\prod_i \Pr_{\text{alt}}(O_{k_i}\mid\cdot)/\Pr_{\text{null}}(O_{k_i}\mid\cdot)$, is a test martingale and
  stays valid under any dependence. In bits, $\log_2 \varepsilon_j = S_{\text{joint}}^{\text{null}} -
  S_{\text{joint}}^{\text{alt}}$, the code-length difference of §4. Without a joint model, the *average*
  of cell e-values is always valid, but conservative.
- **Across cells or bins:** e-BH (Wang & Ramdas 2022) controls FDR under arbitrary dependence. Reject
  the cells with the $\hat k$ largest e-values, $\hat k = \max\{k: \varepsilon_{(k)} \ge K/(\alpha k)\}$.
  This is the right tool for "surprisingly present/absent" calls, because cells within a bin are
  dependent in ways §2 only approximates.
- **Mixture alternatives** for bin-level contamination: $\varepsilon_j = \int \prod_k \frac{\Pr_\pi(O_{jk})}
  {q_{jk}(O_{jk})}\, d\Pi(\pi)$ with a prior Π on contamination fraction is a valid e-value. The
  maximised version (the GLR) is not.

### 3.4 Recommended reporting

| Level | Statistic | Threshold |
|---|---|---|
| Cell | s (bits) = −log2 p; e-value $2^{s}/|\mathcal O|$ or LR | e-BH at 5% over testable cells (§3.2), or lfdr < 0.2 |
| Bin, direction ± | $E^\pm$, $Z^{\pm,\text{dep}}$, saddlepoint tail in bits | BH on bin tail p-values (93-10^4 bins; mildly dependent) |
| Locus | chain-rule s of the block (§2.2) | as a cell |

---

## 4. Binning as coding (MDL)

### 4.1 The code

The data are the contigs of an assembly: per contig x, its sequence composition $t_x$, its per-sample
coverage $d_x = (d_{x1},\dots,d_{xI})$, and the clusters it carries. A binning $\mathcal B$ is a partition
of the contigs. A two-part (or, better, mixture) code describes:

$$ L(\mathcal B, \text{data}) = \underbrace{L(\mathcal B)}_{\text{partition}} + \sum_{b \in \mathcal B}
   \Big[\underbrace{L(\Theta_b)}_{\text{model}} + \underbrace{L(O_b \mid \text{clade}_b, c_b, r_b)}_{\text{content: } S_b}
   + \underbrace{L(d_{X_b} \mid \Theta_b)}_{\text{coverage}} + \underbrace{L(t_{X_b} \mid \Theta_b)}_{\text{composition}}\Big]. $$

- **Partition** $L(\mathcal B)$: $-\log_2$ of a prior on partitions, for example a Chinese restaurant
  process or the binner-agreement prior of `model_spec.md` §7.3.
- **Model** $\Theta_b$: the bin's population parameters: core abundance per sample, clade placement,
  completeness, φ of accessory clusters.
- **Content** $S_b = -\log_2 P(O_b \mid \cdot)$: exactly the bin's total surprise over all K clusters,
  under the *joint* code of §2.2. Absences are part of the bin's code, not of any contig: they cost bits
  for the bin as a whole.
- **Coverage:** the §6/§7.1 likelihood of each contig's read counts given the core profile and its φ.
- **Composition:** a TNF/GC likelihood (optional, `model_spec.md` O10).

**Mixture codes make MDL and Bayes identical.** If each term is a Bayes mixture code (parameters
integrated against their prior rather than transmitted at finite precision), then
$L = -\log_2 P(\text{data}, \mathcal B)$ exactly. A two-part code with ML parameters gives the same to
$O(1)$ per parameter ($\tfrac12\log_2 N$ per free parameter, the BIC/Laplace term; Rissanen 1996 for the
refinement with Fisher information).

### 4.2 Merge, split and move criteria

**Merge A and B** (H1: one population) against keeping them apart (H2):

$$ \Delta L_{\text{merge}} = L(\mathcal B_{\text{merged}}) - L(\mathcal B)
   = \underbrace{\big(S_{A\cup B} - S_A - S_B\big)}_{\text{content}}
   + \underbrace{\big(L^{\text{cov}}_{A\cup B} - L^{\text{cov}}_A - L^{\text{cov}}_B\big)}_{\text{coverage}}
   + \Delta L_{\text{comp}} + \Delta L(\mathcal B). $$

Merge iff $\Delta L_{\text{merge}} < 0$. Since $L = -\log_2 P$,

$$ -\Delta L_{\text{merge}} = \log_2\frac{P(\text{data}\mid H_1)}{P(\text{data}\mid H_2)} + \log_2\frac{P(H_1)}{P(H_2)}
   = \frac{\ln \mathrm{BF} + \ln \text{prior odds}}{\ln 2}, $$

which is the `model_spec.md` §7 rule, with $\ell_{\text{content}} = \ln 2\,([S_A + S_B] - S_{A\cup B})$,
term for term, and τ = 0. A threshold τ > 0 is a description-length margin: insist that the merge saves at
least τ/ln 2 bits. That is the MDL counterpart of controlling false merges when the code is misspecified
(§4.4).

The **coverage Occam term** is automatic. H2 has its own per-sample abundance profile for B, I free
parameters integrated out, about $\tfrac I2 \log_2(\text{reads}/I)$ bits. Merging saves that cost when the
two profiles are proportional, and pays $\sum_i$ of the misfit when they are not.

**Split** is the reverse of merge: the same ΔL with the sign changed, evaluated on a proposed partition of
b (for example by coverage clustering of b's contigs).

**Move X from A to B** (B may be "unbinned"):

$$ \Delta L_{\text{move}} = \big[S_{A\setminus X} + S_{B\cup X}\big] - \big[S_A + S_B\big]
   + \Delta L^{\text{cov}} + \Delta L^{\text{comp}}. $$

Note that $S_{A\setminus X} - S_A$ includes the absences that open up in A when X leaves: removing a contig
that carries expected clusters makes them absent. This is how a core contig defends its place. An
"unbinned" destination is a bin with domain-level prior and its own free coverage profile. Its code cost
is the reference point for "X belongs nowhere here".

### 4.3 Is per-Mb normalisation right?

The `binning_sim.md` rule accepts a piece X into bin A when

$$ \frac{S_{A+X} - S_A}{\mathrm{Mb}_X} \le \frac{S_A}{\mathrm{Mb}_A}. $$

In MDL terms this is $\Delta L_{\text{move}} \le 0$ with one specific, implicit alternative: *X would cost
$\mathrm{Mb}_X \cdot S_A/\mathrm{Mb}_A$ bits anywhere else*, the bin's own average rate. That is why it
works when the alternative is roughly constant per Mb, and why it fails in identifiable ways:

1. **The yardstick inherits the bin's own defects.** $S_A/\mathrm{Mb}_A$ contains A's absence bits. An
   incomplete or contaminated A has a high rate and accepts more (strain_inc base at 1,144 against 727
   bits/Mb for a complete base, `binning_sim.md` §2). The acceptance threshold should not depend on how
   much of A is missing.
2. **Short contigs.** A genome's own short contigs carry its unique (island) genes. Their cost in A is high
   *and* their cost anywhere else is about as high. The per-Mb rule compares against A's average instead
   of the alternative's cost and rejects them (26-73% of own < 10 kb contigs rejected). The MDL rule
   compares $S_{A+X} - S_A$ with $L_{\text{alt}}(X)$, the cost of X in its best other home. For a
   genuine island, both are about $\log_2(n+2)$ per novel gene, so the content term is near 0 and
   coverage decides. That is the README's intended behaviour.
3. **Length scaling is right only for the per-gene part.** Presence costs of foreign genes do scale with
   length (≈ 5.6-6 bits per foreign gene, about 5,000-5,900 bits/Mb). Absence costs do not scale with X's
   length at all: they belong to the bin. Dividing everything by Mb mixes a per-Mb quantity with a
   per-bin one.

**The right normaliser** is the competing code length, not length:

$$ \text{accept } X \to A \iff \big[S_{A+X} - S_A\big] + \Delta L^{\text{cov}} < L_{\text{alt}}(X), \qquad
   L_{\text{alt}}(X) = \min\Big(\min_{B \ne A} [S_{B+X} - S_B],\; S_{\text{domain}}(X)\Big) + \dots $$

When no alternative bin is modelled, the calibrated substitute is the **expected** cost of X under A, i.e.
the excess: accept if

$$ \Delta E = \big(S_{A+X} - S_A\big) - \mathbb E_A\big[S_{A+X} - S_A\big] \le \tau \sqrt{\mathrm{Var}_A(\cdot)}. $$

$\mathbb E_A$ is the expected content cost of a random stretch of A's genome of X's length. It is
estimable from the reference genomes of A's clade (expected rare-cluster density per Mb, by contig-length
class). Per-Mb normalisation is then legitimate, because the comparison is against an expected *rate*
that is specific to the kind of DNA (short contigs are island-rich, so their expected rate is higher).

### 4.4 Misspecification and why τ > 0

MDL selects the binning with the shortest code under *the model's* code. If cells are scored independently
when they are linked (§2), the content code length of a co-binned foreign genome is overstated by its total
correlation, and the content term outweighs coverage. Hence `model_spec.md`'s "do not use τ = 0 until
calibrated". With the chain-rule content code, the over-statement goes away and τ can come down towards
the prior-odds term alone.

### 4.5 Why absence surprise lines up with contamination

**The intended mechanism.** CheckM2's completeness c is part of the bin's code. A contaminant brings
markers the host lacks, so c is overstated. With c overstated, P(A | z=1) ≈ 1 − c (times fragmentation) is
understated, and every truly missing expected cluster costs more bits. In code terms: the absence code uses
the wrong parameter, and the penalty is the KL divergence between the code's c and the c that best
describes the content. That suggests a direct statistic: estimate completeness from the content,

$$ \hat c^{\text{content}}_j = \arg\min_c L(O_j \mid c), \qquad
   D_j = L(O_j \mid c^{\text{CheckM2}}_j) - L(O_j \mid \hat c^{\text{content}}_j) \;\ge 0, $$

a likelihood-ratio statistic in bits for "the markers and the clade's content disagree about how much of
the genome is here". It is large when markers say complete and the content says not: one signature of
contamination by a genome that supplies markers but not the host's accessory content. Model_spec's
duplication term gives the other signature: markers present twice.

**What the data say (§5).** In the current implementation this is *not* where the AUC 0.85 comes from.
ε_jk = 1e-4 + 1 − exp(−r_j G_j ρ_k) includes CheckM2's contamination r_j. Every absent rare family then
costs $-\log_2(1-\varepsilon) \approx 1.44\,\varepsilon$ bits, and there are about 7,000 of them per bin.
The absence total thus contains a term that is a deterministic function of r_j. Recomputing with r_j = 0
drops the AUC from 0.845 to 0.54. The model's own expected absence bits give 0.80, with no observation
involved. The calibrated excess (§1.1), which subtracts that expectation, keeps AUC 0.61-0.71. Only that
part is evidence from the content. The `catalogue_surprise.md` remark that "the contamination checks are
not circular through ε" was tested on presence bits (3.7% change) and does not hold for absence bits.

**Consequences:**

- Never put CheckM2 contamination into a score that is then validated against CheckM2 contamination.
- Use the presence-direction excess $Z^+_j$ with r_j = 0 as the contamination flag (AUC 0.96 against
  CheckM2 ≥ 10% here, §5), with $D_j$ as a complement, and test both against CAMI ground truth rather than
  CheckM2. If contamination must enter ε (for scoring rather than flagging), estimate it from the content
  or from other bins, not from CheckM2.
- In binning, the duplication term and the coverage term carry the contamination evidence. The absence
  term's job is to reward merges of complementary bins.

---

## 5. Numerical check on the catalogue run (93 bins × 7,178 families)

Inputs: `catalogue_surprise/fit_flat.npz` (the main fit) and `inputs.npz`, read-only. Three-state q as in
`analyze.py`. Scripts and logs: `/data/scratch/metapanclust_data/theory_surprise/`. Run time a few
seconds, 1 thread.

**Totals (per bin, three-state, independence):**

| | median | range |
|---|---|---|
| S (bits) | 934 | 158-3,631 |
| H = expected S | 1,045 | 239-2,358 |
| E = S − H | +3 | −630 to +1,303 |
| Z = E/√V | 0.06 | −8.4 to +16.9 (SD 5.7) |
| Absent families, observed vs expected | 6,994 vs 6,980 | |

The model is centred (median Z ≈ 0) but overdispersed (SD 5.7, not 1): dependence and bin-level nuisance
(§2.4) plus the under-prediction of rare presences (§1.4).

**Contamination flag (bins ≥ 70% complete, n = 41, of which 10 are ≥ 10% contaminated):**

| Score | AUC | Spearman ρ with contamination (93 bins) |
|---|---|---|
| Absence bits / expected Mb (as published) | 0.845 | 0.66 |
| … model-*expected* absence bits / expected Mb (no data) | 0.800 | 0.63 |
| … recomputed with r_j = 0 in ε | 0.542 | 0.47 |
| Absence excess, score form, cells with p̂ ≥ 0.5 (r_j in ε) | 0.713 | 0.00 |
| same with r_j = 0 | 0.613 | −0.11 |
| Missing-count O − E, p̂ ≥ 0.5, r_j = 0 | 0.635 | −0.15 |
| $Z^+$ (q-split, §1.1), r_j in ε | 0.842 | 0.23 |
| **$Z^+$, r_j = 0** | **0.955** | 0.39 |
| $E^+$ (bits), r_j = 0 | 0.965 | 0.40 |
| $Z^-$ (q-split), r_j in ε / r_j = 0 | 0.416 / 0.616 | 0.15 / 0.26 |
| Total Z, r_j in ε / r_j = 0 | 0.839 / 0.955 | 0.23 / 0.40 |
| Absence excess bits split by *observed* state (§4.2 as written) | 0.126 | −0.69 |
| Mean ε per bin (pure input) | 1.000 | 0.97 |

Readings:

- The published absence AUC is mostly an input echo (rows 2-3 and the last row).
- The calibrated presence-direction excess is the real contamination signal. With r_j = 0 (no
  CheckM2 input at all), $Z^+$ reaches AUC 0.955, against 0.77 for raw presence bits per Mb in
  `catalogue_surprise.md`. Putting r_j into ε *lowers* it to 0.84, because a large ε explains foreign
  presences away. Ten positives among 41 bins: treat it as a strong hint, to be confirmed on CAMI.
- The absence direction has few cells: a median of 13-18 families per bin have q(seen) > 1/2, so $E^-$
  is a sum of a dozen terms and carries little power. That fits the 2-3-bit absence ceiling.
- Splitting excess by observed state gives an anti-signal (AUC 0.13), as §1.1 predicts.

**Singletons.** Families present in exactly one bin: 4,121 present cells observed, 4,433 predicted by the
leave-one-out model in those columns. In aggregate their presence is fully expected, so their 48% share of
presence *bits* is not excess. Their median share of the presence-direction excess is 40% per bin. They
are worth flagging individually only through the FDR machinery of §3.

**Null calls.** Cells whose observed binary state had probability < 0.01: 4,175 observed, 1,923 expected
under the model (FDR ≈ 0.46 at that cut).

**Dependence among core families.** Families with mean q ≥ 0.6 (10 of them): $K_{\text{eff}}$ = 8.8,
largest eigenvalue 18% of the trace. Their residuals are close to independent.

---

## 6. Implications for the code and the spec

1. **`model_spec.md` §4.2:** replace S+/S− split by observed state with the q-split excess $E^\pm$,
   their variances, and the saddlepoint tail in bits. Report Z only alongside a VIF.
2. **§3 priors:** the prototype's Beta(1,1) is the source of $\log_2(n+2)$. Use the hierarchical
   pseudo-count κ p̂_pa everywhere, including the pangenome survey, and use Kish's effective n when
   weights are uneven.
3. **§4.2 dependence option 2 (pseudo-likelihood):** drop it for totals; it undercounts linked blocks.
   Use a chain-rule/Chow–Liu joint code, with catalogue adjacency as the default dependence structure for
   clusters without reference data. Count the UBA9145 flagellar set as one locus (6-9 bits).
4. **Calls:** cell-level e-values ($2^{s}/|\mathcal O|$ or LR) with e-BH, or FDR from the exact $N_0(t)$.
   Restrict absence calls to testable cells.
5. **Contamination:** remove r_j from ε in any score used to flag contamination. Flag with $Z^+$ (AUC
   0.96 at r_j = 0 here), add the content-completeness statistic $D_j$, and validate on CAMI ground truth.
   Retract the "absence flags contamination" headline of `catalogue_surprise.md`: at r_j = 0 it is
   AUC 0.54.
6. **Binning (§7):** keep the Bayes-factor rule (it is the MDL rule). Replace the bits-per-Mb merge
   criterion of `binning_sim.md` with the alternative-code comparison, or with excess against a
   contig-length-specific expected rate.
7. **Recalibration** (logistic on logit q, held-out cells) before any of the above. The shrunk prior
   under-predicts presence by a factor of 2 in the 0-0.01 bin.

## References

- Dawid 1984, prequential analysis (posterior predictive as code length).
- Rissanen 1996, Fisher information and stochastic complexity; Grünwald 2007, *The MDL Principle*.
- Efron 2004, large-scale simultaneous hypothesis testing (local fdr, theoretical vs empirical null).
- Lugannani & Rice 1980, saddlepoint tail approximation.
- Chow & Liu 1968, tree approximation of discrete joint distributions.
- Kish 1965, design effect / effective sample size.
- Tarone 1990, multiple testing with discrete statistics.
- Vovk & Wang 2021, e-values; Wang & Ramdas 2022, e-BH.
- Itti & Baldi 2009, Bayesian surprise.
