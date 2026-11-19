# Where the surprise-per-Mb merge rule fails, and the optimal test (2026-11-18)

Theory note for the same-genus binning simulation (`docs/experiments/binning_sim.md`,
`docs/experiments/binning_sim_genus.md`, grex run in
`/data/scratch/metapanclust_data/binning_sim_genus/grex_results/results/summary_grex.md`: 200 focal
species, 1,002 congeneric pairs, family cuts c = 0.80 / 0.90 / 0.95). It also bears on the content
term of the binning Bayes factor (`docs/experiments/model_spec.md` §7.2).

**Short answer.**

1. **At the 80% cut, strict_frac only measures gene content.** A two-parameter model of per-gene
   divergence fits all three family cuts across the whole 80-97% ANI range. Its parameters are
   gamma-distributed rates across genes (shape α ≈ 2) and amino-acid divergence ≈ 0.65 × nucleotide
   divergence. Under that model, essentially every shared orthologue clears the 80% cut once ANI is
   above about 88%. strict_frac at c = 0.8 is then just the shared-gene fraction
   g(A) ≈ 0.94 − 1.6(1 − A), which is why it predicts failure better than ANI.
2. **The rule fails where the gap in "new genes per added gene" falls to about one noise unit.**
   That gap is Δq = s_strain(c) − s_cont(A, c). Merge-to-merge noise is about 200-250 bits/Mb,
   roughly constant across ANI and cut, and it comes from island-clustered accessory content. With
   these two terms, AUC ≈ Φ(Δ/σ) predicts the observed per-bin AUC to within 0.01-0.08. The theory
   puts AUC = 0.9 at 92.3 / 94.0 / 95.3-96.4% ANI for c = 0.80 / 0.90 / 0.95. The observed
   breakpoints are 92-94, 94-95 and ≥ 96.
3. **The optimal test uses the identity values, not a hard cut.** The Neyman-Pearson statistic for
   "same species vs congener" scores each gene by its identity to the nearest focal family. Under
   the gamma model it is linear: a count of family-less genes plus the summed amino-acid divergence
   of the shared genes. The divergence part is uniformly most powerful over all congener ANIs.
   Against a contaminant at 95% ANI it needs about 7-26 genes for 5% false accepts at 5% false
   rejects, depending on strain diversity. The 80% hard cut needs about 145 and the 95% cut 23-41.
4. **Implications.**
   - Replace hard membership with an identity-weighted score: the divergence sum, ideally
     normalised by a per-family rate.
   - If a hard cut is kept, set it where the within- and between-species identity densities cross.
     Near the boundary that is c* ≈ 0.97-0.98. Above it, true strains start to look foreign.
   - Expect real power to be lower than the independent-gene numbers by an overdispersion factor:
     accessory content arrives in islands, and divergence in recombination blocks.

Computation: small pure-Python/numpy scripts on `pairs_grex.tsv` and the per-species `sim_*.tsv`, run
under `nice -n 15` on 1 thread. Scripts and outputs are in `/data/scratch/metapanclust_data/theory_merge/`
(`a.py` medians, `fit.py` model fit, `auc.py`/`ds.py` AUC check, `fail.py` failure ANI, `lr2.py` ROC
simulation, `lr2.out`).

---

## 0. Notation

| Symbol | Meaning |
|---|---|
| $A$ | genome-wide ANI (skani, over aligned fraction) between contaminant and focal species; $D = 1-A$ |
| $c$ | family identity cut (amino-acid identity to the family's representative, ≥ 80% coverage) |
| $x_i$ | amino-acid divergence of gene $i$ to its best focal family, $y_i = 1 - x_i$ its identity |
| $g(A)$ | fraction of the genome's genes that have a homologue (orthologue) in the focal pangenome at all |
| $s(A,c)$ | strict fraction: genes assigned to a focal family at cut $c$ (`strict_frac`) |
| $q = 1-s$ | foreign fraction |
| $\mu$ | mean amino-acid divergence of shared genes; $\mu_s$ for a same-species strain, $\mu_A$ for a congener |
| $n$ | training genomes per pangenome (15); foreign-gene cost $B_F = \log_2(n+2) = 4.09$ bits |
| $\rho$ | gene density, ≈ 900-950 genes/Mb |

## 1. Per-gene identity given ANI

### 1.1 Model

A contaminant gene either has a homologue in the focal pangenome, with probability $g(A)$, or is
private to the contaminant's lineage (gain/loss, HGT), with probability $1-g(A)$.

**Rate variation.** For a shared orthologue, nucleotide divergence is $d_i = D\,r_i$ with
$r_i \sim \mathrm{Gamma}(\alpha, \alpha)$ (mean 1, CV $1/\sqrt{\alpha}$). This is the standard
rates-across-loci model: purifying selection varies by gene, and so do recombination and HGT
history.

**Amino-acid divergence.** It is proportional to nucleotide divergence, $x_i = \kappa d_i$. Here
$\kappa < 1$ because most substitutions at these distances are synonymous. Saturation can be ignored
at $D \le 0.2$.

So $x_i \sim \mathrm{Gamma}(\alpha, \alpha/\mu_A)$ with $\mu_A = \kappa(1-A)$, and

$$
F(c;\mu) \equiv \Pr(y_i \ge c) = \Pr(x_i \le 1-c) = P\!\left(\alpha,\ \frac{\alpha(1-c)}{\mu}\right),
$$

where $P$ is the regularised lower incomplete gamma function. The expected strict fraction is

$$
\boxed{\; s(A,c) = g(A)\; P\!\left(\alpha,\ \frac{\alpha(1-c)}{\kappa(1-A)}\right) \;}
\tag{1}
$$

The ratio $s(A,c)/s(A,c')$ does not depend on $g$. So $\alpha$ and $\kappa$ can be fitted from the
ratios $s_{0.90}/s_{0.80}$ and $s_{0.95}/s_{0.80}$ alone (653 pairs with skani ANI, all three cuts).
$g(A)$ then follows as $s_{0.80}/F(0.8;\mu_A)$.

### 1.2 Fit

The least-squares grid fit over pairs gives **α = 2.0, κ = 0.65**. An additive offset for
divergence to the family representative was allowed and fitted to 0. The fitted content term is
close to linear:

$$
g(A) \approx g_s - \beta(1-A),\qquad g_s \approx 0.94,\ \beta \approx 1.6
\tag{2}
$$

That is about 1.6% of genes lost or gained per 1% ANI. skani's aligned fraction rises in step with
ANI, from 19% to 83% over the same bins.

**Two consequences of (1)-(2).**

- **The 80% cut.** At c = 0.8, $F(0.8;\mu_A) \ge 0.976$ for A ≥ 88%, so $s_{0.8}(A) \approx g(A)$.
  The 80% cut sees only gene content, not sequence divergence. That is why the shared-gene fraction
  predicts failure better than ANI: it *is* the quantity the rule measures, while ANI predicts it
  only through the noisy relation (2).
- **Stricter cuts.** These move part of the sequence divergence into the score. Each gene then
  contributes $1 - F(c;\mu_A)$ extra foreign probability. The gain is largest where
  $(1-c) \approx \mu_A$, i.e. where c is close to the congener's mean amino-acid identity.

### 1.3 Predicted vs observed strict_frac (medians per ANI bin)

| ANI bin | median ANI | n | obs s₀.₈₀ | obs s₀.₉₀ | obs s₀.₉₅ | F(0.80) | F(0.90) | F(0.95) | g(A) | pred s₀.₉₀ | pred s₀.₉₅ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 80-85 | 83.97 | 100 | 0.583 | 0.330 | 0.144 | 0.896 | 0.572 | 0.250 | 0.650 | 0.372 | 0.162 |
| 85-88 | 86.86 | 105 | 0.709 | 0.486 | 0.239 | 0.947 | 0.679 | 0.327 | 0.748 | 0.508 | 0.245 |
| 88-90 | 89.01 | 122 | 0.761 | 0.591 | 0.352 | 0.976 | 0.769 | 0.408 | 0.780 | 0.600 | 0.318 |
| 90-92 | 91.10 | 107 | 0.800 | 0.690 | 0.478 | 0.992 | 0.860 | 0.516 | 0.807 | 0.693 | 0.416 |
| 92-94 | 93.02 | 77 | 0.845 | 0.770 | 0.590 | 0.999 | 0.934 | 0.646 | 0.846 | 0.790 | 0.547 |
| 94-95 | 94.56 | 61 | 0.851 | 0.810 | 0.680 | 1.000 | 0.977 | 0.774 | 0.851 | 0.831 | 0.659 |
| 95-96 | 95.37 | 60 | 0.863 | 0.820 | 0.715 | 1.000 | 0.990 | 0.844 | 0.863 | 0.854 | 0.728 |
| ≥ 96 | 96.52 | 20 | 0.881 | 0.837 | 0.743 | 1.000 | 0.999 | 0.935 | 0.881 | 0.880 | 0.823 |
| strain (within species) | ~97-100 | — | 0.940 | 0.921 | 0.892 | | | | 0.94 | | |

Pred s₀.₈₀ equals obs s₀.₈₀ by construction, since g is calibrated on it.

- **Fit quality.** Errors in s₀.₉₀ and s₀.₉₅ are at most 0.06, across a 5-fold range of s₀.₉₅.
- **Residual pattern.**
  - Above 92% ANI, (1) over-predicts the 95%-cut fraction by 0.02-0.08.
  - Same-species strains also lose 5% of their genes between c = 0.8 and c = 0.95. Under α = 2
    that requires $\mu_s \approx 0.025$, an amino-acid divergence to the *representative* of
    about 2.5%, while a strain at 98-99% ANI has $\kappa D \approx 0.01$.
  - Both point to the same missing term: identity is measured to one representative sequence. That
    adds the within-family divergence $\delta_w$ of the focal species,
    $\mu \to \kappa(1-A) + \delta_w$. It matters only when $(1-c)$ is comparable to $\delta_w$,
    i.e. at c = 0.95.
  - A gamma with α = 2 also has a lighter left tail than the data. Some congener genes are much
    more conserved than average: recent recombination within the genus. A two-component rate
    mixture would fit better.

## 2. Surprise per Mb of a strain vs a congener

### 2.1 Expected cost of an added gene

Consider a merge into a near-complete base (`base = full`). Each added gene changes $S$ only if it
brings something the bin does not yet have:

| Added gene lands in | ΔS (bits) |
|---|---|
| family already present in the bin | 0 (duplication is not in $S_{tot}$) |
| new accessory family k (prevalence $p_k = (k+1)/(n+2)$) | $-\log_2 p_k + \log_2(1-p_k\,\mathrm{comp})$ (absence refunded) |
| no family at cut c (foreign) | $B_F = \log_2(n+2)$ |

With $q$ the foreign fraction of the added genes and $a$ the fraction landing in families new to
the bin, at mean cost $\bar B_{acc} < B_F$:

$$
\frac{\Delta S}{\mathrm{Mb}_{add}} \approx \rho\left[q\,B_F + a\,\bar B_{acc}\right].
\tag{3}
$$

The merge rule compares this with $S(\mathrm{bin})/\mathrm{Mb}$, a constant for a given base. The
calibrated rule compares it with a strain quantile. Either way, separation depends only on the
difference between contaminant and strain in (3).

**The gap.**

- **Foreign term.** For a congener, $q_c = 1 - s(A,c)$. For a strain, $q_s(c) = 1 - g_s F(c;\mu_s)$.
- **Accessory term.** Both strains and congeners bring accessory families new to the bin, so $a$
  partly cancels. Congener genes that fall into focal *rare* families pay $\bar B_{acc}$ instead
  of $B_F$, so the effective gap is a fraction $\eta < 1$ of the foreign gap:

$$
\Delta(A,c) \equiv \mathbb{E}\!\left[\tfrac{\Delta S}{\mathrm{Mb}}\right]_{cont}-\mathbb{E}\!\left[\tfrac{\Delta S}{\mathrm{Mb}}\right]_{strain}
\approx \eta\,\rho\,B_F\,\big[s_s(c) - s(A,c)\big].
\tag{4}
$$

With ρ ≈ 940 genes/Mb, B_F = 4.09 bits and η ≈ 0.8 (fitted, see below), $\eta\rho B_F \approx 3{,}100$
bits/Mb per unit of strict-fraction gap.

### 2.2 Why the distributions overlap: island-clustered noise

The per-merge spread of $\Delta S/\mathrm{Mb}$ comes from the grex `sim_*.tsv` files, full base,
mean over focal species. It is

$$
\sigma_s \approx 180\text{-}280,\qquad \sigma_c \approx 160\text{-}260 \quad \text{bits/Mb},
$$

nearly independent of ANI and cut.

**This is far above binomial noise.** A 20% addition has about 700 genes. At q ≈ 0.06 the binomial
SD of q is 0.009, or about 30 bits/Mb. The ratio implies an overdispersion factor

$$
\varphi = \sigma^2_{obs}/\sigma^2_{binom} \approx 40.
$$

Foreign and accessory genes come in islands of tens of genes, so the number of *independent*
units in a merge is roughly the number of islands, not the number of genes. This is the dominant
noise source and it is the same for strains and congeners.

**Gaussian approximation.** With overdispersed noise,

$$
\mathrm{AUC}(A,c) \approx \Phi\!\left(\frac{\Delta(A,c)}{\sqrt{\sigma_s^2+\sigma_c^2}}\right),
\qquad
\Pr(\text{reject}\mid\text{cont, calibrated}) \approx \Phi\!\left(\frac{\Delta - z_{0.95}\,\sigma_s}{\sigma_c}\right).
\tag{5}
$$

**Check against the simulation.** Δ and σ are measured per pair from `sim_*.tsv`; AUC is the
median per ANI bin. The table shows predicted / observed AUC:

| ANI (2-pt bin) | 82 | 84 | 86 | 88 | 90 | 92 | 94 | 96 |
|---|---|---|---|---|---|---|---|---|
| Δ, c = 0.80 (bits/Mb) | 1326 | 1096 | 745 | 646 | 483 | 361 | 290 | 291 |
| AUC c = 0.80 | 1.00 / 1.00 | 1.00 / 1.00 | 0.99 / 0.98 | 0.99 / 0.97 | 0.94 / 0.92 | 0.91 / 0.89 | 0.81 / 0.80 | 0.77 / 0.75 |
| Δ, c = 0.90 | 2310 | 1981 | 1456 | 1225 | 832 | 605 | 370 | 348 |
| AUC c = 0.90 | 1.00 / 1.00 | 1.00 / 1.00 | 1.00 / 1.00 | 1.00 / 1.00 | 0.99 / 0.98 | 0.96 / 0.95 | 0.87 / 0.85 | 0.79 / 0.75 |
| Δ, c = 0.95 | 2874 | 2673 | 2328 | 1995 | 1523 | 1130 | 633 | 475 |
| AUC c = 0.95 | 1.00 / 1.00 | 1.00 / 1.00 | 1.00 / 1.00 | 1.00 / 1.00 | 1.00 / 1.00 | 1.00 / 0.99 | 0.97 / 0.94 | 0.85 / 0.82 |

- **AUC.** The Gaussian form is slightly optimistic, by 0.01-0.03, because σ has heavy tails from
  big islands, but it reproduces the whole curve.
- **Calibrated rejection at 94-95% ANI.** From (5): c = 0.80 gives
  $\Phi(1.26-1.645)=0.35$ (observed 0.35); c = 0.90 gives 0.45 (observed 0.54); c = 0.95 gives 0.80
  (observed 0.94).

**Is the gap theory-predictable?** The theoretical gap (4) from the fitted strict fractions is
$3{,}100 \times (0.940-0.851) = 276$ bits/Mb at c = 0.8 and 94.6% ANI, against 290 observed. At
84% ANI it is 1,110, against 1,096 observed. This match is what fixes η ≈ 0.8. So the gap *is*
predictable from strict_frac alone.

### 2.3 Failure ANI as a function of c and focal strain diversity

**Condition.** Set $\Delta(A^\dagger,c)/\sigma = z^*$ with $\sigma = \sqrt{\sigma_s^2+\sigma_c^2} \approx 300$
bits/Mb. Using (1), (2) and (4):

$$
g_s F(c;\mu_s) - \big[g_s - \beta(1-A^\dagger)\big]\,P\!\left(\alpha,\frac{\alpha(1-c)}{\kappa(1-A^\dagger)}\right)
= \frac{z^*\sigma}{\eta\rho B_F}\ \approx\ 0.097\,z^* .
\tag{6}
$$

**Two regimes.**

- **Loose cut** ($1-c \gg \mu_A$, here c ≤ 0.85). Both $F$ terms are about 1, and
  $\beta(1-A^\dagger) = 0.097 z^*$ gives

  $$1-A^\dagger \approx 0.061\,z^*.$$

  AUC = 0.9 ($z^* = 1.28$) fails at $A^\dagger \approx 92\%$, whatever the cut. Only gene content
  separates the species, and it differs by only ~8% at 95% ANI. The rule cannot reach the species
  boundary.
- **Strict cut** ($1-c \sim \mu_A$). The divergence term adds
  $g\,[1-F(c;\mu_A)] - g_s[1-F(c;\mu_s)]$ to the left side. This pushes $A^\dagger$ up until
  $F(c;\mu_s)$, the strains' own loss at the cut, starts to cancel it.

Numerical solution of (6), with α = 2, κ = 0.65, β = 1.6 and g_s = 0.94:

| c | μ_s (strain aa divergence to family) | s_strain(c) | A† at AUC = 0.9 | A† at AUC = 0.75 | A† at calibrated rejection 50% |
|---|---|---|---|---|---|
| 0.80 | any ≤ 0.025 | 0.94 | 92.3 | 95.9 | 92.4 |
| 0.85 | any ≤ 0.025 | 0.94 | 92.8 | 95.9 | 92.9 |
| 0.90 | 0.005-0.025 | 0.94 | 94.0 | 96.0 | 94.0 |
| 0.95 | 0.005 | 0.94 | 96.4 | 97.3 | 96.4 |
| 0.95 | 0.015 | 0.93 | 96.3 | 97.1 | 96.3 |
| 0.95 | 0.025 | 0.85 | 95.3 | 96.0 | 95.4 |
| 0.97 | 0.005 | 0.94 | 97.7 | 98.2 | 97.7 |
| 0.97 | 0.015 | 0.85 | 97.0 | 97.5 | 97.0 |
| 0.97 | 0.025 | 0.65 | 95.5 | 96.0 | 95.5 |

**Observed breakpoints** (median AUC < 0.9) are 92-94% at c = 0.80, 94-95% at c = 0.90 and ≥ 96% at
c = 0.95. These match the theory, with the c = 0.95 data consistent with μ_s ≈ 0.02-0.025, the
value implied by the strain ratio in §1.3.

**Strain diversity.** It enters in two ways.

1. **Through $F(c;\mu_s)$.** A diverse focal species (large μ_s, or identity measured to a single
   representative) loses its own strains' genes at strict cuts. Once $(1-c) \lesssim \mu_s$, a
   stricter cut no longer helps: at c = 0.97 and μ_s = 0.025, $A^\dagger$ falls back to 95.5%.
2. **Through σ_s and g_s.** More accessory diversity means larger per-merge noise and lower
   prevalence for rare families. Larger pangenomes raise $g_s$ (fewer private genes per strain) and
   lower σ_s. That improves the loose-cut regime only as $\sqrt{\cdot}$ in (6).

**Best hard cut.** Within the gamma model, the cut that best separates one shared gene's identity is
where the two densities cross:

$$
x^* = \frac{\mu_s\mu_A\ln(\mu_A/\mu_s)}{\mu_A-\mu_s},\qquad c^* = 1-x^*.
\tag{7}
$$

At the boundary (A = 95%, μ_A = 0.0325) this gives $c^* = 0.978$ for μ_s = 0.015 and
$c^* = 0.972$ for μ_s = 0.025. For A = 90%, $c^* \approx 0.96$-0.97. A single cut of about 0.95 is a
reasonable compromise across 88-96% ANI: it also keeps strains' family loss small, at
$s_s(0.95) = 0.89$ observed.

## 3. The likelihood-ratio optimal statistic

### 3.1 Per-gene likelihoods

For a gene on the candidate contigs, record its best hit to the focal pangenome as either "no
family" (∅) or an identity $y_i$ to the nearest family member, without a cut. The two hypotheses:

- **H₀ (same species):** a strain with content $g_s$ and divergence $\mu_s$.
- **H₁ (congener at ANI A):** content $g(A)$ and divergence $\mu_A = \kappa(1-A)$.

$$
f_j(\varnothing) = 1-g_j,\qquad
f_j(y) = g_j\,\mathrm{Gamma}\!\left(1-y;\ \alpha,\ \alpha/\mu_j\right),\qquad j\in\{0,1\}.
$$

**Neyman-Pearson statistic.** For independent genes the NP statistic is the log LR:

$$
\Lambda = \sum_i \ln\frac{f_1(y_i)}{f_0(y_i)}
= m_\varnothing \ln\frac{1-g_A}{1-g_s}
+ m\left[\ln\frac{g_A}{g_s} + \alpha\ln\frac{\mu_s}{\mu_A}\right]
+ \alpha\left(\frac{1}{\mu_s}-\frac{1}{\mu_A}\right)\sum_{i\,\text{shared}} x_i
\tag{8}
$$

Here $m_\varnothing$ is the number of family-less genes and $m$ the number of shared genes. The
$(\alpha-1)\ln x$ terms cancel because both hypotheses share α.

**Three properties follow.**

1. **Sufficient statistic.** It is $(m_\varnothing,\ T = \sum_{shared} x_i)$: the count of
   foreign genes and the *summed divergence* of the shared genes. The hard-cut rule replaces $T$
   by a count $\sum 1[x_i > 1-c]$ and merges it into the foreign count with the same weight. That
   discards information whenever $(1-c)$ is not tuned to both μ's.
2. **UMP for the divergence part.** For fixed $m$, the coefficient of $T$ is positive for every
   $\mu_A > \mu_s$, so the family has a monotone likelihood ratio in $T$. A test that rejects for
   large $\bar x = T/m$ is uniformly most powerful against *all* congener ANIs at once. The only
   unknown-A nuisance is the relative weight of the foreign count. A GLRT, or a prior over
   congener ANI (e.g. the genus's representative-ANI distribution), handles that.
3. **Interpretation.** $\bar x/\kappa$ is an estimate of "ANI of the candidate contigs to the
   pangenome", computed only from genes. The optimal test is essentially "estimate ANI to the
   nearest focal family members and compare with the species' own within-species ANI
   distribution", plus a gene-content term that matters below ~90% ANI.

### 3.2 Information per gene and hard-cut efficiency

With Gaussian approximations, the per-gene discrimination of a statistic $t$ is
$d' = (\mathbb{E}_1 t - \mathbb{E}_0 t)/\sqrt{\mathrm{Var}_0 t + \mathrm{Var}_1 t}$. Then
$\mathrm{AUC}(m) \approx \Phi(\sqrt{m}\,d')$, and the gene count for false-accept rate $a$ and
false-reject rate $b$ is

$$
m(a,b) = \left(\frac{z_{1-a}\sqrt{V_0}+z_{1-b}\sqrt{V_1}}{\mathbb{E}_1t-\mathbb{E}_0t}\right)^2 .
\tag{9}
$$

**Hard cut.** $t = 1[\text{foreign at }c]$ is Bernoulli with $p_0 = 1-g_sF(c;\mu_s)$ and
$p_1 = 1-g_AF(c;\mu_A)$. At A = 94.6% and c = 0.8: $p_0 = 0.06$, $p_1 = 0.149$, so

$$
m = \left(\frac{1.645(0.237+0.356)}{0.089}\right)^2 \approx 120
$$

genes for 5% / 5%.

**Divergence mean**, shared genes only. Per-gene variance is gamma plus sampling:
$V_j = \mu_j^2/\alpha + \mu_j/L$, with L ≈ 300 aa. At μ_s = 0.015 and μ_A = 0.0325 (95% ANI):

$$
m = \left(\frac{1.645(0.0127+0.0252)}{0.0175}\right)^2 \approx 13 .
$$

If the per-family rate $r_k$ is known from the pangenome, divide by it: $x_i/r_k$ removes the
gamma term, $V_j \to \mu_j/L$, and $m \approx 3$.

### 3.3 Predicted ROC

Monte Carlo with independent genes from the fitted model (L = 300 aa binomial identity noise, LR
from binned densities at the true A). Each cell gives AUC and the miss rate at 5% false rejects of
strains (`lr2.out`):

| μ_s | contaminant ANI | genes | count at c = 0.80 | count at c = 0.90 | count at c = 0.95 | identity LR |
|---|---|---|---|---|---|---|
| 0.01 | 93% | 10 | 0.78; 0.66 | 0.86; 0.48 | 0.985; 0.07 | 1.00; 0.00 |
| 0.01 | 93% | 50 | 0.96; 0.17 | 0.99; 0.03 | 1.00; 0.00 | 1.00; 0.00 |
| 0.01 | 95% | 10 | 0.71; 0.76 | 0.75; 0.71 | 0.93; 0.31 | 0.999; 0.01 |
| 0.01 | 95% | 50 | 0.91; 0.36 | 0.94; 0.25 | 1.00; 0.00 | 1.00; 0.00 |
| 0.02 | 95% | 10 | 0.71; 0.77 | 0.74; 0.72 | 0.86; 0.49 | 0.92; 0.33 |
| 0.02 | 95% | 20 | 0.79; 0.64 | 0.83; 0.57 | 0.94; 0.24 | 0.98; 0.10 |
| 0.02 | 95% | 50 | 0.91; 0.38 | 0.94; 0.28 | 0.99; 0.03 | 0.999; 0.00 |
| 0.02 | 96% | 20 | 0.75; 0.74 | 0.76; 0.72 | 0.83; 0.58 | 0.88; 0.46 |
| 0.02 | 96% | 100 | 0.95; 0.25 | 0.95; 0.21 | 0.99; 0.07 | 0.997; 0.01 |

Gene counts for 5% false accepts and 5% false rejects (Gaussian approximation (9), independent
genes):

| contaminant ANI | c = 0.80 | c = 0.90 | c = 0.95 | identity LR (μ_s = 0.01 / 0.02) |
|---|---|---|---|---|
| 90% | 39 | 15 | 4-6 | 4 / 6 |
| 93% | 77-80 | 39-41 | 9-15 | 5 / 11 |
| 95% | 145 | 105-109 | 23-41 | 7 / 26 |

**Dominance.**

- **The LR dominates every hard cut at every ANI.** That is guaranteed by the NP lemma for the
  simple alternative, and by property 2 for the composite one.
- **The gain is largest at the boundary,** where the cut count has almost no information
  ($d' \approx 0.19$ per gene at c = 0.8 and 95% ANI) while divergence still has $d' \approx 0.45$-0.8.
- **The c = 0.95 count is the best hard cut,** close to (7), and gets within a factor of 1.5-3 of
  the LR.

**Limits of the model.**

- **Independence.** These are independent-gene numbers, which the simulation shows is optimistic
  for content counts (φ ≈ 40, §2.2). Divergence is also correlated along the genome (recombination
  blocks), but less so than presence/absence of islands. A conservative working assumption is
  effective $m_{eff} = m/\varphi_x$ with $\varphi_x$ of order 3-10 for divergence. Not measured
  here.
- **Within-species tail.** H₀ is not a point. Strains of the focal species span 95-100% ANI, so the
  null distribution of $\bar x$ is a mixture over the species' strain divergences. The test
  threshold has to sit above the *most divergent* strain the species admits. Against a congener at
  96% ANI and a focal species with members at 95-96% from each other, no content or identity test
  can separate them. That is the species-concept limit, not a statistical one. Identity to the
  *nearest member*, not the representative, minimises μ_s and pushes this limit as far as it can
  go.

## 4. Implications

### 4.1 Which family cut

- **c = 0.80 (PPanGGOLiN default) is the worst choice for contamination tests.** Above ~88% ANI it
  measures only gene content, which differs by ~1.6% per 1% ANI. By (6) it fails at ~92% ANI
  whatever the strain diversity.
- **c ≈ 0.95 is the best single hard cut over 88-96% ANI** (eq. 7, table §2.3). Above it (0.97),
  the strains' own genes start falling out ($s_s$ = 0.85 at μ_s = 0.015, 0.65 at μ_s = 0.025). The
  failure point then moves *down* again, along with the bin's completeness and absence terms.
- **The family-building cut and the membership cut should be separate.** Keep the 80% cut for
  *building* the pangenome (families, prevalence, completeness) and use a strict or
  identity-weighted membership only in the *merge test*. Rebuilding families at 95% would shatter
  them and inflate the foreign term for strains as well.

### 4.2 Replace hard membership by identity-weighted scores

- **Score.** Per candidate contig set, compute

  $$\Lambda = m_\varnothing w_\varnothing + \gamma\sum_{shared}\frac{x_i}{\hat r_k},$$

  where $x_i$ is the divergence to the nearest member of the best family and $\hat r_k$ is the
  family's relative rate. $\hat r_k$ can be estimated from the within-pangenome identity
  distribution, or set to 1 if unavailable. Calibrate the threshold on held-out strains, as the
  calibrated rule does now.
- **In the Bayes factor** (model_spec §7.2), this adds an identity term to $\ell_{content}$:

  $$\ell_{id} = \sum_i \ln f_{0}(y_i) - \ln f_{1}(y_i).$$

  $f_0$ is fitted from the pangenome's own leave-one-out identities. $f_1$ is either a congener
  model at the genus's typical ANI or a GLRT over A. It is cheap, because mmseqs already reports
  identity.
- **Keep the foreign count.** It carries most of the information below ~88% ANI, where gene
  content differs by more than 20%.

### 4.3 Small contigs: sample-size limits

Gene density is ≈ 1 gene/kb, so a contig of $\ell$ kb carries about $\ell$ genes, of which ~$0.9\ell$
are shared at the boundary. From (9), assuming independent genes and 5% / 5% error against a 95%-ANI
congener:

| Statistic | Genes needed | Contig length (kb) |
|---|---|---|
| hard cut 0.80 | ~145 | ~150; unreachable with single contigs |
| hard cut 0.95 | 23-41 | 25-45 |
| divergence LR | 7-26 (μ_s 0.01-0.02) | 8-30 |
| divergence LR, rate-normalised | ~3-5 | ~5 |

- **Overdispersion.** These sizes multiply by the overdispersion factor for real genomes. For
  content counts that factor is ~40, from the merge noise in §2.2. That explains why single
  20-50 kb contigs were rejected at only 6-15% at 93-95% ANI with the 80% cut.
- **Contig-level decisions.** Below ~10 kb even the identity LR cannot decide on content alone at
  the species boundary. Such contigs need coverage covariance (model_spec §7.1) or should be left
  to the bin they were assigned to.
- **Error targets.** For error rates $a = b$, $m \propto (2z_{1-a})^2$. Going from 5% to 1% errors
  costs a factor of $(2.33/1.645)^2 = 2.0$ in genes.

### 4.4 What to test next

1. **Recompute the merge AUC with an identity score.** Add per-gene identity and nearest-member
   identity to the mmseqs output in `03_simulate.py`, without changing the search. Then recompute
   AUC with $\bar x$ and with (8) on the existing grex genomes. Prediction: AUC ≥ 0.95 up to 95%
   ANI for 20%-genome merges, and ≥ 0.9 for 50 kb contigs.
2. **Measure φ for divergence.** The per-contig variance of $\bar x$ against its binomial
   expectation sets the real gene counts in §4.3.
3. **Fit the within-family term.** Fit $\delta_w$ (identity to representative vs nearest member)
   and a two-component rate mixture. These should remove the 0.02-0.08 over-prediction of
   s₀.₉₅ above 92% ANI.
