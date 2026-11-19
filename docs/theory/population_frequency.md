# Population frequency of accessory clusters from multi-sample coverage (theory, 2026-11-18)

What can be inferred about the fraction $f$ of a population's cells that carry an accessory gene cluster,
from read depth across samples, and how to tell "member at frequency $f$" from "contaminant from another
population". Builds on [model_spec §6-7](../experiments/model_spec.md) and the numbers in
[coverage_test](../experiments/coverage_test.md).

## Result

- **With nanopore reads the unit of information is reads, not bases.** A contig of length $L$ at summed
  depth $D_+$ carries about $D_+ L / \ell_\text{eff}$ independent counts, where $\ell_\text{eff} =
  \mathrm{E}[\ell^2]/\mathrm{E}[\ell]$ is the length-weighted aligned read length. The sampling variance of a
  depth ratio is then $\mathrm{Var}(\ln \hat f) \approx \ell_\text{eff} / (f D_+ L) + \sigma_b^2$.
- **This reproduces the 10 kb rule.** The spread of ratio_adj on 312 known-member (persistent-family)
  contigs of the two focal bins gives $\ell_\text{eff} \approx 2.8$-$3.1$ kb, consistently in the <5 kb,
  10-20 kb and >20 kb classes. The 5-10 kb class is noisier than counting allows (an apparent 9 kb). At
  $D_+ \approx 4$ the SD of $\ln \hat f$ for a core contig is 0.41 at 3 kb, 0.26 at 10 kb and 0.19 at 20 kb.
  The smallest frequency that can be told apart from 1 (one-sided 5%, power 0.8) is 0.34 at 5 kb, 0.47 at
  10 kb and 0.59 at 20 kb. The observed accessory median, 0.44, is therefore only resolvable on contigs of
  about 10 kb and longer (T2). The contig-specific bias floor $\sigma_b$ is small: ≤ 0.04 in log units,
  from contigs over 50 kb.
- **Membership and frequency separate cleanly.** Given the contig's total count $X_+$, the split of its reads
  across samples does not depend on $f$. So "follows this core at some $f$" vs "follows another profile" is a
  test on that split alone. Its power grows as $X_+ \cdot 2\phi^2$, where $\phi$ is the angle between the two
  square-root abundance profiles. At $\phi = 10°$ a contig needs about 50-90 reads, at $5°$ 200-360, and at
  $2°$ 1,300-2,200 (T3). With 81 of 93 bins at $r \geq 0.8$ to another bin, and all mass in 1-4 samples,
  most pairs sit at small $\phi$. Correlation is weak there for a structural reason. The test's power comes
  from samples where one profile is ~0 and the other is not: the "off" statistic is the limit of this LRT.
- **Composition and gene content should enter as calibrated likelihood ratios,** estimated on the known
  members vs the co-occurring-bin contigs, not as raw distances. They are not independent of each other given
  "contaminant": both track phylogenetic distance. Calibrate them jointly or discount them.
- **Frequency that varies across samples (beta hierarchy) is easier to detect than absolute $f$.** Between-sample
  contrasts cancel the contig's length and GC bias. With L = 20 kb and no biological noise, a change from
  0.3 to 0.6 between two groups of 5 samples needs about 2.5x core depth per sample (T4). Replicate noise
  $\tau$ quickly makes depth useless: adding samples is what helps then.
- **Junction-spanning reads give membership directly, and $f$ without any normalisation.** $\hat f =
  J_\text{ins}/(J_\text{ins}+J_\text{empty})$ at an insertion site is binomial. Its variance $f(1-f)/J$
  depends on depth but not on element length. It is free of abundance, GC, length-class and core-definition
  bias. Correction to the premise: it does **not** need far less depth than coverage in general. Per site
  $J \approx 2 D\,\mathrm{E}[(\ell-2w)_+]/\mathrm{E}[\ell]$, which at the current read lengths is about the
  information of a 3-5 kb contig's coverage (T6). Junctions win for short elements (< ~$\ell_\text{eff}$),
  near $f = 0$ or 1, and above all because a single well-anchored spanning read is linkage evidence that
  coverage can never give. At $\phi \to 0$ that is the only coverage-free route to membership besides
  composition and content. They need anchors longer than the boundary repeats: with IS-flanked islands
  ($w \approx 1.5$ kb), reads must average ≥ 5 kb (T5).

## 1. Generative model, MLE, posterior, detection limit

### 1.1 Model

Sample $i = 1..n$, population $b$ with core coverage (cells × genome copies per base) $\theta_i$. Contig $k$ of
mappable length $L_k$ is carried at $f_k$ copies per cell, so $f_k \in [0,1]$ for single-copy accessory
content, $f = 1$ for core, and $f > 1$ for multi-copy. Reads of length $\ell \sim F_\ell$ start uniformly
at rate $\theta_i / \mathrm{E}[\ell]$ per base. The read count overlapping $k$ is

$$
X_{ik} \sim \mathrm{Poisson}(\lambda_{ik}), \qquad
\lambda_{ik} = \frac{\theta_i\, f_k\, b_k\, (L_k + \bar\ell - 2m)}{\mathrm{E}[\ell]},
$$

where $b_k$ is the contig's multiplicative bias (GC, mappability, end effects, MAPQ filtering) and $m$ is
the minimum overlap for a read to be assigned. Overdispersion from library and mapping noise is a
per-sample gamma factor, so $X_{ik}$ is negative binomial. Conditioned on the total, it is beta-binomial
(model_spec §6.2). The core $C$ (length $L_C \gg L_k$) gives $Y_i$ with $f = 1$ and $b_C$ averaged to ~1.

**Depth vs counts.** jgi depth is the mean per-base depth $\bar d_{ik}$. With reads of fixed length
$\ell$ the per-base depths at $x, y$ have covariance $(\theta/\ell)(\ell - |x-y|)_+$, so

$$
\mathrm{Var}(\bar d) = \frac{\theta}{\ell L^2}\iint (\ell-|x-y|)_+\,dx\,dy = \theta\,g(L/\ell),\qquad
g(x) = \begin{cases} 1 - x/3 & x \le 1 \\ 1/x - 1/(3x^2) & x > 1\end{cases}
$$

$$
\mathrm{CV}^2(\bar d) = \frac{g(L/\ell)}{\theta} \;\xrightarrow{L \gg \ell}\; \frac{\ell_\text{eff}}{\theta L},
\qquad \ell_\text{eff} = \frac{\mathrm{E}[\ell^2]}{\mathrm{E}[\ell]} \text{ for a read-length distribution.}
$$

A contig shorter than a read behaves like a single point: $\mathrm{CV}^2 \to 1/\theta$, however short it is.
Counting read midpoints instead of averaging depth gives $\mathrm{CV}^2 = \mathrm{E}[\ell]/(\theta L)$, smaller by
$1 + \mathrm{CV}_\ell^2$, about 2x for nanopore length distributions. This is a free gain (§Implications).

### 1.2 MLE

Under $H_\text{const}$ ($f$ the same in all samples), condition on $N_i = X_i + Y_i$:

$$
X_i \mid N_i \sim \mathrm{Bin}(N_i, \pi), \quad \pi = \frac{f\rho}{1+f\rho}, \quad \rho = \frac{b_k L_k}{L_C}.
$$

The binomials share $\pi$, so the sufficient statistic is $(X_+, N_+)$ and

$$
\hat\pi = X_+/N_+,\qquad \hat f = \frac{\hat\pi}{\rho(1-\hat\pi)} = \frac{X_+}{\rho\, Y_+}.
$$

The MLE is the **ratio of sums** over samples, not a mean of per-sample ratios. Samples are weighted by
their counts automatically, and samples where the population is absent contribute nothing (this is
what restricting to "informative samples" approximates). $\rho$ needs $b_k$: ratio_adj's division by the
median of trusted contigs of the same length class is a plug-in estimate of $b_k$ as a function of length.
That calibration is exactly what makes $\hat f$ comparable across length classes.

Fisher information (Poisson limit, $L_C \gg L_k$, $Y_+$ effectively known):

$$
I(f) = \frac{\mu_1}{f}, \quad \mu_1 = \rho Y_+ \approx \frac{D_+ L_k}{\ell_\text{eff}}, \qquad
\mathrm{Var}(\hat f) \approx \frac{f}{\mu_1}, \qquad
\mathrm{Var}(\ln \hat f) \approx \frac{1}{f\mu_1} + \sigma_b^2 .
$$

$\mu_1$ is the number of reads the contig would get at $f = 1$. $D_+ = \sum_i \theta_i$ is summed core depth.
$\sigma_b^2$ is the residual variance of $\ln b_k$ after length-class calibration. It does **not** shrink
with depth or with samples, because $b_k$ is a property of the contig shared by all samples. With
beta-binomial overdispersion $\omega$, per-sample variance is inflated by $1 + (N_i - 1)/(\omega + 1)$, which
caps the information a single sample can give at about $\omega$ reads.

### 1.3 Posterior

Prior from the pangenome (model_spec §6.3): $f \sim \pi_0\delta_0 + \pi_1\delta_1 + (1-\pi_0-\pi_1)\,
\mathrm{Beta}(p\eta, (1-p)\eta)$, with $p$ the species prevalence and $\pi_0 = (1-p)^{N_\text{eff}}$,
$\pi_1 = p^{N_\text{eff}}$. In the Poisson limit the likelihood is $f^{X_+}e^{-f\mu_1}$, so the continuous part of
the posterior is

$$
p(f \mid X) \propto f^{X_+ + p\eta - 1}(1-f)^{(1-p)\eta - 1} e^{-\mu_1 f}, \quad f \in (0,1),
$$

with normaliser $B(a', \beta)\,{}_1F_1(a'; a'+\beta; -\mu_1)$ (Kummer), $a' = X_+ + p\eta$, $\beta = (1-p)\eta$.
The point masses get weights $\pi_0\,\mathbb 1[X_+ = 0]$ (more precisely $\pi_0 \cdot P(X_+ \mid \text{background}\ \nu)$)
and $\pi_1 \mu_1^{X_+} e^{-\mu_1}$. The grid in model_spec §6.4 computes this directly. For intuition: when the
prior is weak, the posterior mean is $\approx (X_+ + p\eta)/(\mu_1 + \eta)$. The prior counts as $\eta$
pseudo-reads at $f = 1$-equivalent. At $\mu_1 \approx 7$ (a 5 kb contig at $D_+ = 4$) a prior of $\eta = 2$ is
22% of the information, so the prior matters on exactly the contigs where coverage is weak.

### 1.4 Detection limits

**From 0.** $P(X_+ = 0 \mid f) = e^{-f\mu_1} \le 0.05$ gives $f_\text{lo} = 3/\mu_1 = 3\ell_\text{eff}/(D_+L)$. This is
model_spec's $\varphi_\text{min}$ with $R_C$ counted in reads, not bases (a 20x difference for 3 kb reads). With a
cross-mapping background of $\nu$ reads, this becomes the Poisson test of $\mu_1 f + \nu$ vs $\nu$:
$f_\text{lo} \approx (3 + 2\sqrt{\nu})/\mu_1$.

**From 1.** One-sided test of $f = 1$ at level $\alpha$ with power $1-\beta$:

$$
-\ln f \;\ge\; z_\alpha\sqrt{\tfrac{1}{\mu_1}+\sigma_b^2} + z_\beta\sqrt{\tfrac{1}{f\mu_1}+\sigma_b^2}.
$$

**Why < 10 kb is noise.** Take $f = 0.5$ and $\sigma_b \to 0$: $\sqrt{\mu_1} \ge (1.645 + 0.84\sqrt 2)/0.69 = 4.1$, so
$\mu_1 \ge 17$ reads, and $L \ge 17\,\ell_\text{eff}/D_+ \approx 13$ kb at $D_+ = 4$, $\ell_\text{eff} = 3$ kb.

### 1.5 Check against the coverage test (known members)

312 contigs carrying persistent families (UBA9145 + MAG-120531). SD of $\ln$ ratio_adj (robust, MAD) vs
contig length. $\ell_\text{eff}$ solves $\mathrm{SD}^2 = g(L/\ell)/D$ at the class's median length and depth.
Script: `/data/scratch/metapanclust_data/theory_popfreq/fit.py`.

| Length | Contigs | SD ln ratio_adj | ratio_adj < 0.6 | median depth | fitted $\ell_\text{eff}$ |
|---|---|---|---|---|---|
| < 5 kb | 81 | 0.45 | 15% | 3.1 | 3.1 kb |
| 5-10 kb | 84 | 0.45 | 20% | 3.7 | 9.1 kb (excess) |
| 10-20 kb | 88 | 0.21 | 7% | 4.6 | 2.9 kb |
| ≥ 20 kb | 59 | 0.14 | 0% | 5.2 | 2.8 kb |

(>50 kb: SD 0.03-0.04 on 5 contigs. This bounds $\sigma_b$.)

- Three of four classes give the same $\ell_\text{eff}$ of about 3 kb. That is a plausible
  length-weighted *aligned* read length for this environmental nanopore data (mean ≈ 1.6 kb for a log-normal
  with $\sigma = 0.8$). It is an upper bound, because overdispersion is folded in. Reading the actual
  read-length distribution from the BAMs would make it a test rather than a fit.
- The 5-10 kb class has about 3x the variance that counting predicts. Candidates are misassembled
  short contigs, repeat-collapsed ends, and the end-effect correction (raw ratio falls with length)
  being too coarse at the 5-7 / 7-10 kb boundary. Short contigs are not just noisy, they are also worse
  calibrated.
- Predicted false-accessory rate at the 0.6 cutoff for $f = 1$ (T1, $D_+ = 4$): 11% at 3 kb, 7% at 5 kb,
  2.6% at 10 kb, 0.4% at 20 kb. Observed: 15%, 20%, 7%, 0%. Same order, with heavier tails than Gaussian
  below 20 kb.

**T1. SD of $\ln\hat f$ at f = 1 / P($\hat f$ < 0.6 | f = 1)**, depth-averaging estimator, $\ell_\text{eff}$ = 3 kb,
$\sigma_b$ = 0.04.

| $D_+$ \ L | 3 kb | 5 kb | 10 kb | 20 kb | 50 kb |
|---|---|---|---|---|---|
| 4 | 0.41 / 0.107 | 0.35 / 0.071 | 0.26 / 0.026 | 0.19 / 0.004 | 0.13 / 0 |
| 10 | 0.26 / 0.025 | 0.22 / 0.011 | 0.17 / 0.001 | 0.13 / 0 | 0.09 / 0 |
| 20 | 0.19 / 0.003 | 0.16 / 0.001 | 0.12 / 0 | 0.09 / 0 | 0.07 / 0 |
| 50 | 0.12 / 0 | 0.11 / 0 | 0.08 / 0 | 0.07 / 0 | 0.05 / 0 |
| 100 | 0.09 / 0 | 0.08 / 0 | 0.07 / 0 | 0.06 / 0 | 0.05 / 0 |

**T2. Detection limits.** $f_\text{hi}$ is the largest f distinguishable from 1 (one-sided 5%, power 0.8).
$f_\text{lo} = 3/\mu_1$ is the smallest f distinguishable from 0 (no background).

| $D_+$ \ L | 3 kb | 5 kb | 10 kb | 20 kb | 50 kb |
|---|---|---|---|---|---|
| 4 | 0.26 / 0.75 | 0.34 / 0.45 | 0.47 / 0.23 | 0.59 / 0.11 | 0.71 / 0.045 |
| 10 | 0.47 / 0.30 | 0.54 / 0.18 | 0.63 / 0.09 | 0.72 / 0.045 | 0.80 / 0.018 |
| 20 | 0.60 / 0.15 | 0.65 / 0.09 | 0.72 / 0.045 | 0.79 / 0.022 | 0.84 / 0.009 |
| 50 | 0.73 / 0.06 | 0.76 / 0.036 | 0.81 / 0.018 | 0.84 / 0.009 | 0.88 / 0.004 |
| 100 | 0.79 / 0.03 | 0.82 / 0.018 | 0.85 / 0.009 | 0.87 / 0.005 | 0.89 / 0.002 |

At $D_+ = 4$, a 5 kb contig cannot be told from 0 below f ≈ 0.45 and cannot be told from 1 above 0.34. No
frequency in between is resolved: at that length and depth, coverage alone carries no information about $f$.
$f_\text{hi}$ saturates below 0.9 because of $\sigma_b$. Frequencies above ~0.85 cannot be told from core by coverage at
any depth unless $\sigma_b$ is calibrated better (junctions do not have this floor, §4).

## 2. Identifiability: member at f vs contaminant

### 2.1 Factorisation

Let the target core profile be $p_i = \theta_i/\theta_+$ and a candidate source profile be $q_i$ (another bin's core,
or free). Write $X_i \sim \mathrm{Poisson}(\mu\,r_i)$, with $r = p$ under $H_\text{same}$ and $r = q$ under $H_\text{cont}$. Then

$$
P(X \mid H) = \underbrace{\mathrm{Poisson}(X_+ \mid \mu)}_{\text{scale: } f \text{ or contaminant abundance}}
\times \underbrace{\mathrm{Mult}(X \mid X_+, r)}_{\text{shape: membership}} .
$$

The scale parameter is free under both hypotheses ($f$ under $H_\text{same}$, the contaminant's relative abundance
under $H_\text{cont}$), so the marginal likelihoods share the first factor up to priors. The Bayes factor for
membership comes only from how the reads are split across samples:

$$
\ln \Lambda = \sum_i X_i \ln\frac{p_i}{q_i}, \qquad
\mathrm{E}[\ln\Lambda \mid H_\text{same}] = X_+\,\mathrm{KL}(p\|q), \quad
\mathrm{E}[\ln\Lambda \mid H_\text{cont}] = -X_+\,\mathrm{KL}(q\|p).
$$

**The scale is not useless.** The prior on $f$ (pangenome prevalence) vs the prior on contaminant abundance
enters through the first factor. A contig at ratio 0.02 of a core fits "rare contaminant" better than
"accessory at f = 0.02" when the pangenome prior on f puts little mass there.

### 2.2 Power as a function of the angle

Map profiles to the unit sphere: $u = \sqrt p$, $v = \sqrt q$, $\cos\phi = \sum_i\sqrt{p_iq_i}$ (Bhattacharyya). For
nearby profiles $\mathrm{KL} \approx 2\|u-v\|^2 = 4(1-\cos\phi) \approx 2\phi^2$, and the per-read variance of the LLR is
$\approx 4\phi^2$. The two hypotheses' LLR distributions are then $\pm 2\phi^2 X_+$ apart with SD $2\phi\sqrt{X_+}$, so

$$
X_+ \;\ge\; \frac{(z_\alpha + z_\beta)^2}{4\phi^2} \quad\text{(known alternative)},
$$

When the alternative is free (G-test of homogeneity, df = $n_\text{inf}-1$), the non-centrality is
$\lambda \approx X_+\sum_i (q_i - p_i)^2/p_i \approx 4\phi^2 X_+$.

**T3. Reads on the contig needed** (5%, power 0.8).

| Angle $\phi$ | Known alternative | Free profile, df 1 | Free profile, df 3 |
|---|---|---|---|
| 30° | 6 | 7 | 10 |
| 20° | 13 | 16 | 22 |
| 10° | 51 | 64 | 89 |
| 5° | 203 | 258 | 358 |
| 2° | 1,269 | 1,610 | 2,237 |

**Why 81/93 is weak.** A correlation of 0.8 between depth profiles dominated by one sample does not
translate into a fixed $\phi$. What matters is how much of each profile's mass sits where the other's is
small. When both put > 90% of their mass on barcode04, $\phi$ is a few degrees, and a 5-10 kb contig
(7-15 reads at $D_+ \approx 4$) cannot reach the hundreds of reads that T3 requires.

**The useful limit is the "off" statistic.** If $p_i \approx 0$ in some samples where $q_i > 0$, then
$\ln(p_i/q_i) \to -\infty$ and a few reads there are decisive. With a background (cross-mapping, index hopping)
rate $\varepsilon$, set $p_i \leftarrow (1-\varepsilon)p_i + \varepsilon/n$. Each read in an absent sample then contributes
about $\ln(\varepsilon/(n q_i))$, typically -3 to -6 nats. The 17 barcode16 contaminant modules are this case:
almost all their reads sit where $p_i < 0.01$. In the opposite direction, a contig with no reads outside the
target's samples gains little ($\mathrm{KL}(p\|q)$ is small when $q$ also covers those samples). This is why
"follows the core" (r 0.98) is weak evidence for membership, while "has reads where the core is absent" is
strong evidence against it.

**Unidentifiable case.** At $\phi = 0$ (two populations with identical abundance across samples, as with
co-occurring partners in one sample), "member at f" and "contaminant at relative abundance f" have
identical coverage likelihoods for any number of reads. Only non-coverage evidence decides: composition,
content, linkage (§4). This is the regime of most of the 93 bins.

**Overdispersion.** With Dirichlet-multinomial concentration $\omega$, the effective read count is
$X_+(\omega+1)/(X_+ + \omega)$. That caps the attainable LLR at about $\omega\cdot 2\phi^2$ whatever the depth. Estimate $\omega$ from
split-core halves (model_spec §6.2).

### 2.3 Combining evidence

$$
\ln\mathrm{BF} = \ell_\text{shape} + \ell_\text{scale} + \ell_\text{comp} + \ell_\text{content} + \ell_\text{link},
\qquad \text{posterior odds} = \mathrm{BF}\times\text{prior odds}.
$$

- **Addition is valid only for conditionally independent evidence given each hypothesis.** Under
  $H_\text{same}$ the evidence lines are roughly independent: coverage noise is technical, composition and content are
  genome properties. Under $H_\text{cont}$ composition (TNF distance) and content (share of genes in the species'
  public families) both decrease with the contaminant's phylogenetic distance, so they are positively
  dependent. A closely related contaminant fools both. Fix: estimate the joint density
  $P(\text{tnf}, \text{content} \mid H)$ on the two calibration sets the coverage test already has: known members
  (persistent contigs) and contaminant-like contigs (co-occurring bins' contigs scored as members). Or, more
  simply, use $\max$ rather than sum for the two phylogenetic terms, plus a fixed discount.
- **Selection bias.** Binners chose these contigs partly because their TNF and coverage matched. Conditional
  on being binned, $H_\text{cont}$ contigs look member-like. Likelihoods for $H_\text{cont}$ therefore have to come from
  binned contaminants (the calibration set above), not from random genomes. The same applies to
  $\ell_\text{comp}$: model_spec's concern about double counting is really a concern about using an uncalibrated
  $P(\text{tnf} \mid H_\text{cont})$.
- **Content term, per gene, calibrated on the coverage test.** Per gene, the share in public UBA9145
  families is 216/290 = 0.745 on contigs that follow the core and 1/237 = 0.004 on off-population contigs.
  A gene in a public family contributes $\ln(0.745/0.004) \approx +5.2$ nats, and a gene not in one
  $\ln(0.255/0.996) \approx -1.4$. Genes on one contig are not independent (operons, shared ancestry of the island).
  Count per module or use an effective number of genes (e.g. divide by the mean module size). Novel
  accessory islands are exactly the members with low content scores, so this term must never veto alone:
  give it a floor (e.g. $\ge -3$ nats per contig).
- **Gene-content surprise from the main model.** $\ell_\text{content} = \ln 2\,(S_\text{rest}(X) + S_A - S_{A+X})$
  (model_spec §7.2) already is a Bayes factor. Use it as is, with the same dependence caveat against composition.
- **Prior odds.** The base rate of contamination among accessory-like contigs is estimable. On the
  coverage-test classes, off-population (4.2% of bp) vs accessory-like (8.1%) suggests prior odds of order
  1:2 among low-ratio contigs. Use a per-bin estimate where possible.

## 3. Frequency varying across samples

### 3.1 Hierarchical beta model

$$
f_i \sim \mathrm{Beta}(\bar f\psi, (1-\bar f)\psi), \qquad X_i \mid f_i \sim \mathrm{Poisson}(f_i\, m_i), \quad m_i = \rho Y_i ,
$$

(or the binomial form of §1.2). The marginal is Poisson-beta. $\mathrm{Var}(f_i) = \bar f(1-\bar f)/(\psi+1)$, and
$\psi\to\infty$ recovers $H_\text{const}$. Environmental structure enters as $\mathrm{logit}\, \bar f_i = \beta_0 + \beta^\top z_i$
(beta regression) with $\psi$ absorbing the residual.

**The bias floor cancels.** $b_k$ is multiplicative and shared by all samples, so it drops out of every
between-sample contrast $\ln(f_i/f_j)$. The same holds for length-class calibration and GC. Heterogeneity is
limited only by counting and by biological replicate noise. This does not hold for the absolute value
of $f$.

### 3.2 Heterogeneity tests

- **Global.** LRT of $\psi = \infty$ vs $\psi < \infty$, with a boundary null: $\tfrac12\chi^2_0 + \tfrac12\chi^2_1$. Or Tarone's
  score test of extra-binomial variation. In practice the Pearson dispersion $\chi^2 = \sum_i (X_i - \hat f m_i)^2/(\hat f m_i)$ on $n_\text{inf}-1$ df,
  which is fast and adequate at these counts. Its technical-overdispersion baseline comes from core
  split-halves: test against the core's own dispersion, not against Poisson.
- **Contrast** (two sets of samples $A$, $B$, e.g. two environments):
  $\ln(\hat f_B/\hat f_A) \sim N\big(\ln(f_B/f_A),\ \tfrac{1}{f_A\mu_A}+\tfrac{1}{f_B\mu_B}+\tfrac{2\tau^2}{n}\big)$,
  where $\mu_g = n\,d\,L/\ell_\text{eff}$ is the reads at $f = 1$ in group $g$ ($n$ samples at per-sample core depth $d$), and
  $\tau$ is the between-replicate SD of $\ln f$ within a group (biology plus unmodelled technical noise).
- **Heterogeneity confounds membership.** A contig whose $f_i$ varies may be a member under selection,
  or a contaminant with a profile at angle $\phi$. These are the same model when $f_i$ is free:
  $f_i \propto q_i/p_i$. Only the beta prior's smoothness (moderate $\psi$) and the other evidence lines separate them.
  So a significant heterogeneity test should send the contig to the §2 test, not straight to "selection".

### 3.3 Minimal design

**T4. Per-sample core depth d needed** to detect $f_A \to f_B$ (two groups of n samples each,
two-sided 5%, power 0.8, L = 20 kb, $\ell_\text{eff}$ = 3 kb). "-" means not reachable at any depth.

| change | $\tau$ | n = 2 | 3 | 5 | 10 | 20 |
|---|---|---|---|---|---|---|
| 0.3 → 0.6 | 0 | 6.1x | 4.1x | 2.5x | 1.2x | 0.6x |
| 0.3 → 0.6 | 0.2 | 18x | 7.2x | 3.3x | 1.4x | 0.7x |
| 0.3 → 0.6 | 0.4 | - | - | - | 2.6x | 0.8x |
| 0.5 → 0.7 | 0 | 18x | 12x | 7.1x | 3.6x | 1.8x |
| 0.5 → 0.7 | 0.2 | - | - | - | 8.0x | 2.5x |
| 0.1 → 0.3 | 0 | 6.5x | 4.3x | 2.6x | 1.3x | 0.7x |
| 0.1 → 0.3 | 0.4 | - | 14x | 4.5x | 1.6x | 0.7x |

For other lengths, scale $d$ by $20\,\text{kb}/L$ (for $L \gtrsim \ell_\text{eff}$, counting-limited rows only).

- **Total reads, not samples, set the counting term.** For $\tau = 0$ the required $n\cdot d$ is constant (about 12x
  per group for 0.3 → 0.6 at 20 kb). Once $\tau > 0$ there is a minimum number of samples:
  $n > 2\tau^2(z_{\alpha/2}+z_\beta)^2/\ln^2(f_B/f_A)$. That is 2 samples per group for 0.3 → 0.6 at $\tau = 0.2$, and 6 for 0.5 → 0.7.
  Beyond it, extra depth stops helping.
- **The current data cannot test heterogeneity.** The focal populations have core > 0.2x in 2-4 samples, and
  ~4x in one. Even $\tau = 0$ needs ≥ 2 samples per group at ≥ 6x on 20 kb contigs. A design that can test
  it: ≥ 5 samples per condition where the population is at ≥ 3x each, and contigs ≥ 20 kb (or junction counts, §4).

## 4. Long reads: junction-spanning evidence

### 4.1 Estimator

Consider an accessory element of length $A$ inserted at a site in a core region. The "empty" allele is the
same site without it, joining left flank $F_L$ directly to right flank $F_R$. Count, per sample, reads that
align across the left junction with at least $w$ bp on each side:

- $J^\text{ins}_i$: reads that run $F_L \to$ element. Plus the same at the right junction.
- $J^\text{emp}_i$: reads that run $F_L \to F_R$ directly, with ≥ w on each side.

All of these reads come from cells of the population that have the flank, so

$$
J^\text{ins}_i \mid J_i \sim \mathrm{Bin}(J_i, f_i), \qquad \hat f_i = \frac{J^\text{ins}_i}{J_i}, \qquad
\mathrm{Var}(\hat f_i) = \frac{f_i(1-f_i)}{J_i},
$$

Equal opportunity requires the same anchor $w$ on both alleles, and the element longer than $w$.
Requiring reads to span the whole element ($\ell > A + 2w$) instead would make the insertion allele rarer
in proportion to $\mathrm{E}[(\ell-A-2w)_+]$ and bias $\hat f$ low.

**Exactly what cancels:** the abundance $\theta_i$, the core definition and its length $L_C$, the contig's GC
and length bias $b_k$, the length-class calibration, and sample-specific mapping rates. These all affect both
alleles equally because both share the flank. What remains: allele-specific mappability (keep $w$ above the
boundary repeat length) and chimeric reads (rate $\varepsilon_c$, roughly $10^{-3}$-$10^{-2}$ for ligation libraries, which
adds false $J^\text{ins}$ at rate $\varepsilon_c$ times junction-adjacent reads).

### 4.2 Information vs coverage

Expected junction reads per site per unit depth: a read covers a point with probability proportional to its
length, and spans it with ≥ w each side if the point falls in its middle $\ell - 2w$. So

$$
\mathrm{E}[J_i] = \theta_i\,u(w), \qquad u(w) = \frac{\mathrm{E}[(\ell-2w)_+]}{\mathrm{E}[\ell]},
\qquad J = 2\,D_+\,u(w) \text{ for both junctions.}
$$

Fisher information for $f$: junctions give $J/(f(1-f))$, coverage gives $\mu_1/f = D_+L/(f\,\ell_\text{eff})$. The ratio is

$$
\frac{I_\text{junc}}{I_\text{cov}} = \frac{2u(w)\,\ell_\text{eff}}{(1-f)\,L}.
$$

- **Short elements ($L \lesssim 2u\,\ell_\text{eff}/(1-f)$, about 6 kb at f = 0.5 on current data) and high f: junctions win.** As $f\to 1$ the binomial variance vanishes,
  while coverage hits its $\sigma_b$ floor. Junctions are the only way to tell f = 0.95 from 1.
- **Long elements at mid f: coverage integrates more reads.** A 50 kb island at f = 0.5 gets ~9x more
  information from coverage than from its two junctions at 3 kb reads.
- **The claim "far less depth" is therefore true for what matters most here, not for precision in general.**
  Junctions are unbiased and need no calibration. They work for elements of any length, including the
  < 10 kb contigs where coverage is noise. Above all, they give linkage: a single spanning read with anchors
  ≥ w at > 95% identity is evidence of co-residence with likelihood ratio about $1/\varepsilon_c$, so 5-7 nats. That
  is more than the coverage shape test gives at $\phi < 10°$ with < 50 reads. The membership question at
  $\phi \approx 0$ (§2.2), unanswerable by coverage, becomes answerable with a handful of reads.

**T5. Usable fraction $u(w)$** of junction-covering reads, log-normal read lengths ($\sigma$ = 0.8).

| mean read | $\ell_\text{eff}$ | w = 250 | 500 | 1,000 | 2,000 |
|---|---|---|---|---|---|
| 2 kb | 3.8 kb | 0.76 | 0.56 | 0.31 | 0.11 |
| 3 kb | 5.7 kb | 0.83 | 0.68 | 0.45 | 0.22 |
| 5 kb | 9.5 kb | 0.90 | 0.80 | 0.63 | 0.39 |
| 10 kb | 19 kb | 0.95 | 0.90 | 0.80 | 0.63 |
| 20 kb | 38 kb | 0.98 | 0.95 | 0.90 | 0.80 |

**T6. SD of $\hat f$ at f = 0.5:** junctions (both, $J = 2D_+u$) vs coverage on an element of length L
(depth-averaging, $\sigma_b$ = 0.04).

| reads | w | $D_+$ | J | junction SD | coverage SD, L 3 / 10 / 20 / 50 kb |
|---|---|---|---|---|---|
| mean 1.6 kb ($\ell_\text{eff}$ 3.0 kb, current data) | 500 | 4 | 3.8 | 0.26 | 0.29 / 0.19 / 0.14 / 0.09 |
| | 500 | 20 | 19 | 0.11 | 0.13 / 0.08 / 0.06 / 0.04 |
| | 1,500 | 4 | 1.0 | 0.49 | (same) |
| | 1,500 | 20 | 5.1 | 0.22 | |
| mean 5 kb ($\ell_\text{eff}$ 9.5 kb) | 500 | 4 | 6.4 | 0.20 | 0.33 / 0.29 / 0.22 / 0.15 |
| | 500 | 20 | 32 | 0.09 | 0.15 / 0.13 / 0.10 / 0.07 |
| | 1,500 | 4 | 3.9 | 0.25 | (same) |
| | 1,500 | 20 | 20 | 0.11 | |

- **Longer reads make coverage worse and junctions better.** At fixed depth, longer reads mean fewer
  independent coverage counts, while junction reads gain anchor room. At mean 5 kb, junctions beat
  coverage for every element up to ~20 kb.
- **The anchor length is set by the boundary repeats.** Islands are commonly flanked by IS elements
  (0.8-2.5 kb), tRNA att sites, or rRNA operons. An anchor must extend past the repeat into unique sequence,
  so w ≈ repeat length + 200-500 bp. With w = 1.5-2 kb, reads with mean < 3 kb lose most of their junction
  counts (T5). Read lengths of ≥ 5 kb mean (N50 ≥ 10 kb) are where the method becomes routinely usable.
  Current data are about 1.6 kb mean aligned (if the §1.5 fit is right), which is enough only for
  repeat-free boundaries (w ≈ 500).
- **Detection from 0 with junctions:** $P(J^\text{ins} = 0) = (1-f)^J \le 0.05$ gives $f_\text{lo} \approx 3/J$, independent of L.

### 4.3 Where the empty allele comes from

- **Within-contig islands.** The assembly kept the insertion allele. The empty allele's reads have split
  alignments $F_L \to F_R$ with a deletion of size $\approx A$ (minimap2 reports this as a large D or as a
  supplementary alignment). These are counted directly from the BAM.
- **Accessory content on its own contig** (common: strain variation breaks contigs). The junction reads are
  reads whose primary/supplementary alignments link a core contig end to the accessory contig. These are the
  GFA links of the assembly graph, weighted per sample. The empty allele is a read linking the two core
  contigs directly. Graph alignment (GraphAligner to the Flye/metaMDBG GFA) counts both natively.
- **Element absent from the assembly** (f too low to assemble): the empty allele still gets counted, and
  $J^\text{ins} = 0$ bounds f from above. The coverage method cannot do this.

### 4.4 Joint model

Junction counts and coverage are close to independent: junction reads are a small subset of all reads and
sit at the boundaries. The joint log-likelihood in $f_i$ is

$$
\ell(f_i) = \ln\mathrm{Bin}(J^\text{ins}_i \mid J_i, f_i) + \ln\mathrm{BetaBin}(X_i \mid N_i, \pi(f_i), \omega),
$$

under the §3 hierarchy. The junction term also calibrates $b_k$. On contigs with both kinds of evidence,
$\ln(\hat f_\text{cov}/\hat f_\text{junc})$ estimates $\ln b_k$, and it can replace the length-class median.
A persistent difference (coverage f > junction f) means the element sits at more than one site (copy
number > 1, or a transposable cassette), not that it is at higher frequency.

## Implications for the pipeline

- **Report $\hat f$ as a ratio of summed reads, with a posterior.** Never use per-sample ratio means. Use
  read-midpoint counts per contig (from the BAM, MAPQ ≥ 20) instead of jgi mean depth: about 2x the
  information for free at nanopore read-length spreads.
- **Gate accessory calls by $\mu_1$, not by contig length.** Require $\mu_1 = D_+L/\ell_\text{eff} \ge 17$ for "f < 1 at
  f ≈ 0.5", and $\mu_1 \ge 3/f$ for "present". At current depth ($D_+ \approx 4$, $\ell_\text{eff} \approx 3$ kb) that is
  L ≥ 13 kb, which matches the empirical 10 kb rule. At 20x it would be 2.5 kb.
- **Membership test: use the shape LLR with an explicit background $\varepsilon$, not correlation.** Report $\phi$ (Bhattacharyya
  angle) to the nearest other bin for every bin. Where $\phi < 5°$, mark coverage-based membership as
  unidentifiable and rely on junctions, composition and content.
- **Combine evidence as calibrated Bayes factors** with the known-member and co-occurring-contaminant sets
  as training data, modelling composition and content jointly.
- **Add junction counting** (BAM split/supplementary alignments, or GraphAligner on the GFA) for every
  accessory-like contig with a core-contig neighbour in the graph. This is the cheapest upgrade: no new
  sequencing, and it answers the flagellar-locus question (f ≈ 0.5 and linked to the core?) directly,
  along with the chimeric-junction caveat in the coverage test.
- **For new sequencing,** longer reads (≥ 5 kb mean) help more than more depth for membership. For $f$
  heterogeneity, ≥ 5 samples per condition with the population at ≥ 3x each.

## Files

- Scripts and outputs: `/data/scratch/metapanclust_data/theory_popfreq/`. `fit.py` and `check.py` (spread vs
  length on known members), `tables.py`, `tables2.py` (T1-T6), `tables.txt`, `tables2.txt`. All
  ran in seconds, one thread.

## Open

- Read the aligned read-length distribution from the BAMs to test $\ell_\text{eff} \approx 3$ kb directly. If
  it is much longer, the excess variance is overdispersion and $\omega$ is the limiting parameter.
- Explain the 5-10 kb excess variance before trusting any contig < 10 kb, even at higher depth.
- Simulate (strain mixtures with known f, nanopore read-length model) to check T2-T6, including the
  chimera rate's effect on junction false positives.
