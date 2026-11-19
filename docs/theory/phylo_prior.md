# A phylogenetic gain/loss prior for cluster surprise (theory note, 2026-11-18)

**Summary.**

- Model each cluster's presence as a two-state gain/loss Markov chain on the tree. The covariance of
  presence between two genomes is then $\pi(1-\pi)e^{-\mu d}$, where $d$ is their patristic distance
  and $\mu$ is the cluster's turnover rate.
- The hierarchical beta-binomial (`model_spec.md` §3) matches this covariance exactly on an ultrametric
  "star of stars" tree. Its rank concentrations are then fixed by geometry:
  $\kappa_\ell \approx h_\ell/(h_{\mathrm{pa}}-h_\ell)$ for slow clusters and $\kappa_\ell \approx e^{2\mu h_\ell}-1$
  for fast ones, with $h$ the crown heights of the ranks.
- With GTDB RED depths, the slow-cluster values are $\kappa \approx 2$-$3.5$ for class to family and
  $\approx 0.5$ for genus. These are the values fitted in Chesterfield (3) and PanGBank B (0.5).
- The beta-binomial is only the *second-moment* (linear) approximation of the chain. It fails in three
  ways:
  - it cannot jump: a conjugate Beta adds counts linearly, while a binary ancestral state flips after
    one observation;
  - it uses one $\kappa$ per rank for clusters whose turnover rates span orders of magnitude;
  - it treats rank groups as exchangeable stars with fixed depth.
- Together these explain why genus shrinkage loses to species-only empirical Bayes in PanGBank once
  $m \ge 5$. Conspecific genomes *screen off* the genus, and the leftover pseudo-count $\kappa(1-\bar p)$
  costs bits on thousands of clusters.
- The right prior is the chain itself, run by Felsenstein pruning with detection-model tip vectors, per
  cluster rates pooled by nonparametric EB over rate classes, and an up-down pass that gives exact
  leave-one-out for free.
- A light check on the Chesterfield bins agrees: how informative a relative is depends on cluster
  frequency in the way rate heterogeneity predicts, and genus pairs are no more concordant than family
  pairs (§5).

Notation: genome $g$, cluster $k$, true state $z_{gk}\in\{0,1\}$, observation $O_{gk}\in\{I,P,A\}$ (intact,
partial, absent). The tree $\mathcal T$ has branch lengths in some time-like unit; §4 discusses the choice
of unit.

---

## 1. The gain/loss process

### 1.1 Assumptions

- **A1.** For each cluster $k$, presence evolves along every branch as an independent continuous-time
  Markov chain on $\{0,1\}$: gain rate $g_k$ ($0\to1$), loss rate $l_k$ ($1\to0$).
  - Gains by HGT are treated as "gain from anywhere", so the donor is not modelled.
  - Copies are ignored.
- **A2.** Clusters are conditionally independent given their rates. Level 2 (context) is what models
  their dependence.
- **A3.** Optional lineage heterogeneity: the rates on branch $b$ are $r_b\,g_k$ and $r_b\,l_k$. Here
  $r_b$ is a clade-level turnover multiplier shared by all clusters (open vs. streamlined lineages).
- **A4.** At the root, the state has distribution $\pi_k$ (stationary). §1.6 covers the non-stationary
  case.

Write $\mu=g+l$ for the turnover rate and $\pi=g/\mu$ for the stationary prevalence. The transition matrix
over a branch of length $t$ (with $e\equiv e^{-\mu t}$) is

$$
P(t)=\begin{pmatrix} 1-\pi(1-e) & \pi(1-e)\\ (1-\pi)(1-e) & \pi+(1-\pi)e\end{pmatrix},
\qquad
P_{11}(t)-\pi = (1-\pi)e^{-\mu t}.
$$

### 1.2 Induced covariance

At stationarity, take two tips $i,j$ whose path through their MRCA has length $d_{ij}$. By time
reversibility, the path is a single chain of length $d_{ij}$:

$$
\Pr(z_i=1,z_j=1)=\pi P_{11}(d_{ij}),\qquad
\operatorname{Cov}(z_i,z_j)=\pi(1-\pi)\,e^{-\mu d_{ij}},\qquad
\operatorname{Corr}(z_i,z_j)=e^{-\mu d_{ij}}.
$$

A relative at distance $d$ therefore carries information on a scale of $1/\mu$. Three regimes follow:

- **slow / core** ($\mu d \ll 1$): even phylum-level relatives inform;
- **fast / accessory** ($\mu d_{\text{genus}} \gtrsim 2$): only conspecifics, or only the same strain,
  inform;
- **mixtures of rates**: averaging over clusters of unequal rate gives
  $\bar\rho(d)=\int e^{-\mu d}\,dF(\mu)$. This curve is log-convex: it falls fast at short distances and
  then has a long tail. A single fitted exponential is therefore wrong at both ends.

With lineage multipliers (A3), $d_{ij}$ is replaced by $\sum_{b\in\text{path}} r_b t_b$.

### 1.3 Posterior predictive by pruning

Each tip carries a likelihood vector $L_g=(L_g(0),L_g(1))$. For an error-free genome this is $(1,0)$ or
$(0,1)$; for a MAG see §3. Upward (Felsenstein) pass:

$$
L_n(s)=\prod_{c\in \mathrm{ch}(n)} \sum_{s'} P_{ss'}(t_c)\,L_c(s'),\qquad
Z_k=\sum_s \pi_{\text{root}}(s)\,L_{\text{root}}(s).
$$

The downward (outside) pass gives, for each node $n$ with parent $p$ and sibling(s) $\sigma$:

$$
U_n(s)=\sum_{s_p} U_p(s_p)\,\Big[\prod_{\sigma}\sum_{s''}P_{s_p s''}(t_\sigma)L_\sigma(s'')\Big]\,P_{s_p s}(t_n),
\qquad U_{\text{root}}=\pi_{\text{root}}.
$$

- **New genome.** A query $x$ is placed (by GTDB-Tk / pplacer) on edge $(p,n)$, at distance $u$ below $p$,
  with pendant length $t_x$. Split the edge there, at a new node $a$. Then

$$
\Pr(z_x=1\mid \text{data})=\frac{\sum_{s} \alpha_a(s)\,P_{s1}(t_x)}{\sum_s \alpha_a(s)},\qquad
\alpha_a(s)=\Big[\textstyle\sum_{s_p}U'_p(s_p)P_{s_p s}(u)\Big]\Big[\textstyle\sum_{s'}P_{s s'}(t_n-u)L_n(s')\Big],
$$

  where $U'_p$ is the outside message at $p$ times the messages from $p$'s other children.
- **Exact leave-one-out for free.** For a tip $g$ already in the tree, its LOO predictive is
  $\Pr(z_g=1\mid \text{rest}) \propto \sum_s U_{\text{pa}(g)}'(s)\,P_{s1}(t_g)$. Its own vector never enters.
  This removes the leak found in `catalogue_surprise.md` with no refitting. To leave out a 99.5% ANI
  group, set the vectors of its members to $(1,1)$ and re-run only the path to the root, at
  $O(\text{depth})$ cost.
- **Surprise.** Surprise is $-\log_2 \Pr(O_x\mid \text{data})$, with
  $\Pr(O_x)=\sum_z \Pr(z_x=z\mid\text{data})\Pr(O_x\mid z)$. The detection-free "bits present" of §4 in
  `catalogue_surprise.md` carry over unchanged.

**Linear (kriging) view.** To second order, the predictive is the best linear predictor

$$
\Pr(z_x=1\mid z)\approx \pi + c_x^\top C^{-1}(z-\pi\mathbf 1),\qquad C_{ij}=\pi(1-\pi)e^{-\mu d_{ij}}.
$$

The weights $C^{-1}c_x$ concentrate on the closest relatives and decay once those relatives "screen" the
more distant ones. On a tree the screening is exact under the Markov property. Given the state at the
crown of $x$'s species, nothing outside the species matters. The beta-binomial (§1.4) reproduces this
linear predictor and nothing more.

### 1.4 The beta-binomial as a moment-matched approximation

Take an ultrametric "star of stars". Each rank-$\ell$ node $t$ sits at height $h_\ell$, all rank-$\ell$
nodes have equal height, and the children of a node attach directly to it. Two tips whose deepest shared
node is $t$ then have $\operatorname{Corr}=e^{-2\mu h_t}$.

In the HBB, $\theta_t\mid\theta_{\mathrm{pa}}\sim\mathrm{Beta}(\kappa\theta_{\mathrm{pa}},\kappa(1-\theta_{\mathrm{pa}}))$
gives $E[\theta_t(1-\theta_t)\mid\theta_{\mathrm{pa}}]=\theta_{\mathrm{pa}}(1-\theta_{\mathrm{pa}})\,\kappa/(\kappa+1)$.
So $H_t\equiv E[\theta_t(1-\theta_t)]=\pi(1-\pi)\prod_{\ell\le t} q_\ell$ with
$q_\ell=\kappa_\ell/(\kappa_\ell+1)$. Two tips with deepest shared node $t$ then have

$$
\operatorname{Cov}(z_i,z_j)=\operatorname{Var}(\theta_t)=\pi(1-\pi)-H_t
\;\Rightarrow\;
\operatorname{Corr}_{\mathrm{HBB}}(t)=1-\prod_{\ell\le t}q_\ell .
$$

Set this equal to $e^{-2\mu h_t}$ at every rank and write $x_\ell=2\mu h_\ell$:

$$
\boxed{\;\kappa_\ell=\frac{1-e^{-x_\ell}}{e^{-x_\ell}-e^{-x_{\mathrm{pa}}}}\;}
\qquad
\begin{cases}
\mu h_{\mathrm{pa}}\ll 1: & \kappa_\ell\to \dfrac{h_\ell}{h_{\mathrm{pa}}-h_\ell}\quad\text{(rate-free: pure geometry)}\\[2ex]
\mu h_{\mathrm{pa}}\gg 1: & \kappa_\ell\to e^{2\mu h_\ell}-1\quad\text{(grows exponentially with turnover)}
\end{cases}
$$

So the HBB is the exact second-moment match of the chain on a rank-ultrametric star-of-stars. It is
wrong whenever any of these conditions fails:

1. **Non-ultrametric or variable-depth ranks.**
   - GTDB RED normalises rank depth only approximately (interval about ±0.05-0.1 around the median).
     For slow clusters $\kappa$ depends on $h_\ell/h_{\mathrm{pa}}$, which then varies by a factor of
     2-4 between clades of the same rank. For fast clusters $\kappa\sim e^{2\mu h}$ is exponentially
     sensitive to $h$.
   - Inside a rank node the tree is not a star. Two species in a genus can be sisters or sit at the
     genus crown. The HBB treats them as exchangeable and so averages over the very structure that
     matters for fast clusters.
2. **Rate heterogeneity across clusters.** Per-cluster $\kappa_\ell$ ranges from the geometric floor
   $h_\ell/(h_{\mathrm{pa}}-h_\ell)$ for core clusters to $\gg 100$ for fast ones. One $\kappa_\ell$ per rank is
   a compromise that is too large for slow clusters and too small for fast ones. The lognormal $\kappa$
   hyperprior of §3.5 helps, but only if a cluster's rate can be identified from its own data.
3. **Discreteness (the non-linear part).**
   - The chain's ancestral states are binary. Under the chain, one present genome in a species with
     slow within-species turnover sets the species-crown state to 1. A second conspecific is then
     present with probability $\approx 1-l\,t$, which does not depend on the genus.
   - The conjugate Beta instead updates linearly. After $x=m$ presences its predictive is
     $1-\kappa(1-\bar p)/(m+\kappa)$, an error that shrinks as $1/m$ rather than with the branch length.
   - Matching moments cannot fix this, because the HBB's marginal at each rank is a Beta while the
     chain's is close to a two-point mixture.
4. **Non-stationarity and single origins.** A cluster born once (a lineage-specific island) has
   $\pi_{\mathrm{root}}=0$ and a single gain event: Dollo-like. Its prevalence pattern is a clade
   indicator, which the HBB can only approximate with $\kappa\to0$ at that rank.

### 1.5 Which ranks inform

A relative at shared-node height $h$ carries a share $e^{-2\mu h}$ of the maximum correlation.

| Turnover | Informative relatives | What the HBB should do |
|---|---|---|
| $\mu h_{\text{phylum}} < 0.3$ (core, slow shell) | all ranks | shrink at every rank, $\kappa_\ell \approx h_\ell/(h_{\mathrm{pa}}-h_\ell)$ |
| $\mu h_{\text{family}}\sim 1$ | order, family, genus | shrink at family and genus, deep ranks give only $\pi$ |
| $\mu h_{\text{genus}} > 2$ (accessory, MGE) | conspecifics, near-clones | species only, everything above sets $\pi$ |

### 1.6 Non-stationary root and rate classes

- **Root.** Replace $\pi_{\text{root}}$ by a free root prior $\rho_k$, pooled across clusters (§3.3).
- **Dollo-type clusters.** For clusters plausibly born once, use a two-component mixture over clusters:
  stationary, and "born on some branch $b$ with probability $\propto t_b$, then losses only". The
  second component's likelihood also follows from pruning, by summing over the branch of birth (an
  $O(N)$ extra pass).

---

## 2. Explaining the empirical findings

Approximate GTDB RED medians per rank: phylum 0.33, class 0.47, order 0.63, family 0.77, genus 0.92.
These are from memory of Parks et al. 2018 and should be checked against r226. With node height
$h=1-\mathrm{RED}$:

| Rank | $h_\ell$ | slow-limit $\kappa_\ell=h_\ell/(h_{\mathrm{pa}}-h_\ell)$ | fitted (Chesterfield) | fitted (PanGBank B) |
|---|---|---|---|---|
| class | 0.53 | 3.8 | 3 | – |
| order | 0.37 | 2.3 | 3 | – |
| family | 0.23 | 1.6 | 3 | – |
| genus | 0.08 | 0.53 | 30 | – |
| species (crown, given genus) | ~0.02-0.03 | ~0.3-0.6 | uninformative | **0.5** (median) |

### 2.1 Why phylum..family get $\kappa\approx3$

At these ranks most of the catalogue's informative families are slow ($\mu h_{\mathrm{pa}}\lesssim1$),
so the fitted $\kappa$ sits near the geometric value of 2-4. The chain predicts something specific: the
fitted value should depend on the depth ratio of adjacent ranks and **not** on the cluster mix. This can
be tested on a larger tree (§4.3).

### 2.2 Why genus adds little in Chesterfield ($\kappa_{\text{genus}}=30$)

There are three reasons, in decreasing order of importance.

1. **Construction filter.** A catalogue locus needs Spacedust support from *another genus*
   (`catalogue_surprise.md`, circularity). Under the chain, genus rank adds information only through
   events on the branches between the family crown and the genus crown. The cluster gained there and
   confined to one genus is exactly what the filter removes. What survives is either slow (family
   already tells you) or fast, with $\kappa_{\text{genus}}\approx e^{2\mu h_g}-1\gg1$.
2. **Sampling.** Most genera have 1-2 bins, so the $\kappa_{\text{genus}}$ likelihood is flat. The
   optimiser drifts to large $\kappa$, which is equivalent to "use the family".
3. **Variable genus depth.** Chesterfield genera are dominated by deep marine clades (SAR11-like,
   Pseudothioglobaceae). Their genus-pair distances overlap family-pair distances (check in §5).

### 2.3 Why genus shrinkage hurt in PanGBank ($m\ge5$ conspecifics)

- **Screening off.** On the tree, the new genome and its genus-mates are conditionally independent given
  the state at the species crown $c$:

  $$\Pr(z_x\mid z_{\text{species}},z_{\text{genus}})=\sum_s \Pr(z_c=s\mid z_{\text{species}},z_{\text{genus}})\,P_{s z_x}(h_s).$$

  Congeners enter only through $\Pr(z_c=s\mid\cdot)$. For slow-in-species clusters ($\mu h_s\ll1$), $m=5$
  concordant conspecifics already pin $z_c$: the posterior odds against the observed state shrink like
  $(\mu h_s)^m$. Congeners then change the prediction by $O((\mu h_s)^m)$, far below the 1-2% gains
  that were being sought. For fast-in-species clusters the congeners are themselves uninformative,
  because $e^{-2\mu h_g}\approx0$.
- **The pseudo-count cost.** The HBB posterior predictive $(x+\kappa\bar p_g)/(m+\kappa)$ with
  $\kappa=0.5$ assigns
  - $\kappa\bar p_g/(m+\kappa)$ to a cluster absent in all $m$ conspecifics (congener_only, with
    $\bar p_g\sim0.5$): 0.045 at $m=5$, against $\approx g\,h_s\approx 0$ under the chain;
  - $\kappa(1-\bar p_g)/(m+\kappa)$ to the absence of a cluster present in all $m$ (focal_only,
    $\bar p_g\approx0$): 0.09 at $m=5$, or 0.14 bits per test genome per cluster.

  Over thousands of genus clusters, these leftover pseudo-counts are the measured losses (focal_only
  +253 bits, congener_only +73 bits at $m=5$).
- **Why `sp_eb` wins.** `sp_eb` fits one Beta(a, b) over all clusters of the focal species. Its a and b
  are well below 1, so its effective pseudo-count is $b$ (or $a$) $\ll\kappa=0.5$. A Beta(a, b) with
  $a,b\to0$ tends to the two-point mixture that the chain implies for slow clusters. Of the two
  predictors, `sp_eb` is the better approximation of the discrete-state predictive. It is not a better
  use of the data.
- **What $\kappa_{\text{genus}}\approx0.5$ means.**
  - The slow-limit formula gives $\kappa=h_s/(h_g-h_s)=0.5$ at $h_s/h_g=1/3$. Most genus clusters
    therefore have $\mu h_g\lesssim1$: they are fixed (0 or 1) within each species and differ between
    species. A gain/loss event happened on the species stems, between genus crown and species crown.
  - The PanGBank text calls these clusters "species markers". That is correct, and the right reading is
    that the HBB at this rank is in the geometric regime: it says *nothing* about rate except that
    $\mu h_s\ll 1$.
  - Equivalently, in the fast-parent formula $e^{2\mu h_s}-1=0.5$ gives $\mu h_s\approx0.2$: at most
    about one turnover event per five species-crown depths.
  - Even the oracle $\kappa$ loses. No single $\kappa$ fixes the linear-update error, which is a
    discreteness error, not a variance error.
- **Where genus should help.** It should help only when the species state is uncertain:
  - $m\in\{0,1,2\}$;
  - all conspecifics are MAGs with low detection probability (§3);
  - a cluster's species-crown state is ambiguous (mixed $x$).

  This is the "shared" stratum, where the 50/50 mix did win slightly.

### 2.4 The 16% Chesterfield gain is real but deep

Most of the 16% taxonomic gain comes from slow and mid-rate families informed at order or family level
(§5 table). That is the regime where the HBB and the chain agree to second order. Fast families
(singletons, rare) contribute almost all the remaining loss, and the HBB cannot reduce it. The chain
cannot either without closer relatives. The chain predicts that the **gain grows with sampling density
near the query, not with the number of ranks**.

---

## 3. MAG detection inside the tree

### 3.1 Tip likelihood vectors

From `fragmentation.md` and `model_spec.md` §2, for genome $g$ (completeness $c_g$, contig lengths
$L_{g\cdot}$) and cluster $k$ (span $s_k$):

$$
L_{gk}(z)=\Pr(O_{gk}\mid z):\qquad
\begin{array}{c|cc}
 O & z=1 & z=0\\\hline
 I & c_g\,[1-\beta(1-f_g(s_k))] & \varepsilon_{gk} f_g(s_k)\\
 P & 1-P_I-P_A & \varepsilon_{gk}(1-f_g(s_k))\\
 A & \approx(1-c_g)f_g + (1-c_g)^2(1-f_g) & 1-\varepsilon_{gk}
\end{array}
$$

These vectors plug directly into §1.3.

- An absence in a 40%-complete MAG has $L\approx(0.6,1)$, a weak vote for 0.
- An intact observation has $L\approx(c,\,10^{-4})$, a near-certain 1.
- Pruning needs no other change. Uncertain tips are standard: they are how ambiguous characters are
  handled in phylogenetics.
- Contamination enters only through $\varepsilon$, as before.
- Optional composition loss $\lambda_k$: multiply the $z=1$ column for $O=A$ by $(1+\lambda_k/P_A)$, i.e.
  add $\lambda_k$ to $P_A$ and renormalise.

### 3.2 Fitting rates: EM over expected transitions

The complete data are the full histories on every branch. For a two-state chain the sufficient
statistics per branch are:

- the expected numbers of gains $\hat G_b$ and losses $\hat L_b$;
- the expected dwell times $\hat T_b^{0}, \hat T_b^{1}$.

Each comes in closed form from the up/down messages at the branch's ends: Holmes and Rubin 2002; Minin
and Suchard 2008. For one branch of length $t$, endpoint posterior $\Pr(s_p=a,s_n=b\mid\cdot)$, and the
conditional expectations $E[N_{01}\mid a,b,t]$ and $E[T^{0}\mid a,b,t]$, these are explicit functions of
$g,l,t$. The M-step is

$$
\hat g_k=\frac{\sum_b \hat G_{bk}}{\sum_b r_b \hat T^{0}_{bk}},\qquad
\hat l_k=\frac{\sum_b \hat L_{bk}}{\sum_b r_b \hat T^{1}_{bk}},\qquad
\hat r_c=\frac{\sum_{k}\sum_{b\in c}(\hat G_{bk}+\hat L_{bk})}{\sum_k\sum_{b\in c}(g_k\hat T^{0}_{bk}+l_k\hat T^{1}_{bk})}.
$$

Detection parameters ($\beta$, the $\varepsilon$ scale, $\lambda$) can be added to the M-step through the
expected tip states $\Pr(z_g=1\mid\text{all})$, which the same up-down pass provides.

### 3.3 Empirical-Bayes pooling of rates

Most clusters are seen in a handful of genomes. Their own $(g_k,l_k)$ MLE is unstable (often $\hat l=0$ or
$\hat g\to\infty$). Pool them as follows.

1. Use a fixed grid of $R$ turnover classes $\mu_r$ (log-spaced, e.g. $R=8$ spanning $10^{-2}$ to $10^{2}$
   per unit tree height) and $Q$ prevalence classes $\pi_q$ (logit-spaced, $Q=8$). This gives
   $RQ=64$ classes.
2. One up-pass per cluster evaluates $Z_k(r,q)$ for all classes at once (vectorised over the class
   axis).
3. Estimate the mixing weights $w_{rq}$ by nonparametric EB (NPMLE, Kiefer-Wolfowitz):

$$
\gamma_k(r,q)=\frac{w_{rq}Z_k(r,q)}{\sum_{r'q'}w_{r'q'}Z_k(r',q')},\qquad w_{rq}\leftarrow\frac1K\sum_k\gamma_k(r,q).
$$

   This iterates on a $K\times64$ matrix and needs no further tree passes. Covariates (MGE flag, span,
   module membership, persistent-partition share) enter as a multinomial-logit model for $w_{rq}(x_k)$,
   the analogue of §3.5's $\kappa$ regression.
4. The predictive is the class mixture
   $\Pr(z_x=1\mid\cdot)=\sum_{rq}\gamma_k(r,q)\Pr(z_x=1\mid\cdot,r,q)$. Rate uncertainty is integrated
   out, not plugged in.
5. Re-estimate the lineage multipliers $r_c$ by the M-step of §3.2 (clade-level, pooled over all
   clusters), then recompute $Z_k$. Two or three outer iterations are enough.

The marginal likelihood also compares models directly: $\sum_k\log\sum_{rq}w_{rq}Z_k(r,q)$ for chain
vs. HBB on the same tips.

---

## 4. Scalable algorithm

### 4.1 Tree and units

- Use the GTDB r226 bac120/ar53 species-representative trees (≈ 1.4×10⁵ leaves together).
- Inside each species, attach genomes as a subtree built from skani ANI: UPGMA on $1-\mathrm{ANI}$,
  scaled to the tree's units by the species-representative pendant lengths. A star at the species crown
  is the cheap fallback.
- Branch lengths are marker substitutions per site. Turnover is not clock-like with respect to them.
  The lineage multipliers $r_c$, at family level (≈ 10⁴ parameters estimated from ~10⁵ clusters), absorb
  this.
- Replicate weighting (§3.2 of the spec) becomes unnecessary for prediction, because clonal genomes
  sit on near-zero branches and are screened off automatically. It is still useful for runtime.

### 4.2 Per-cluster passes with sparsity and caching

**Absent-subtree cache.**

- A subtree in which cluster $k$ is unobserved everywhere has an upward message that depends only on
  the subtree, the class $(r,q)$, and the tips' absence vectors.
- Those vectors depend on $k$ only through span $s_k$. Bin spans into $S=4$ classes.
- Cache $M_n(r,q,S)$ for every node: $2N\times64\times4\times2$ doubles ≈ 0.4 GB for $N=10^5$.
- In one bottom-up pass also cache the outside messages for the all-absent configuration.

**Steiner pass.**

- For cluster $k$ with present set $V_k$ (size $n_k$), recompute messages only on the Steiner tree
  $\mathcal S_k$ = union of the root-paths of $V_k$. At each node of $\mathcal S_k$, children outside it
  contribute cached messages.
- For near-core clusters, use the complement: cache the all-present configuration and recompute along
  the paths of the absent tips. The work is $\min(n_k, N-n_k)$.

**Cost.**

- Let $\nu$ be the total number of non-default cells: $\sum_k\min(n_k,N-n_k)$.
- With $\sim3\times10^3$ clusters per genome and most of them either rare or near-core, $\nu$ is about
  $1$-$3\times10^8$ for $10^5\times10^5$.
- Each non-default tip adds at most $\bar\delta$ new nodes (mean depth above the point where its path
  merges into $\mathcal S_k$). For clade-clustered presence $\bar\delta$ is about 5-20; GTDB's mean
  depth of ~60 is the worst case.
- Each node costs $2\times2\times64$ multiply-adds per child, plus the $P(t)$ evaluation, which is a
  shared table per (branch, class) = $2N\times64$.
- Total: $\nu\bar\delta\times\sim600$ flops $\approx 10^{12}$-$10^{13}$ flops. That is a few to tens of
  CPU-hours, or minutes on a GPU, embarrassingly parallel over clusters. Memory per worker is $O(N)$
  plus the shared cache.
- The EM outer loop is 3-5 iterations for $r_c$. The NPMLE inner loop works on the $K\times64$ matrix of
  $Z_k$ only (seconds).

**Scoring a new MAG.**

- Place the MAG with GTDB-Tk.
- For every cluster, the outside message at the attachment edge combines two parts: the path from the
  attachment point up to where it meets $\mathcal S_k$, which uses cached all-absent messages, and the
  stored outside message at that meeting node.
- The cost is $O(K\cdot\text{depth}\cdot64)\approx6\times10^{8}$ flops, about 1 s per MAG.
- Storing outside messages only on $\mathcal S_k$ costs $\nu\bar\delta\times64$ floats, ≈ 0.1-1 TB at
  full scale. The cheaper route stores $\gamma_k$ and recomputes the Steiner pass for the clusters
  needed. Even cheaper: store the posterior state at each Steiner node (2 floats per class collapsed by
  $\gamma_k$, ≈ 10 GB) and accept the small error from collapsing classes.

**Approximations, in the order to try them.**

1. $R\times Q=64\to16$ classes.
2. Collapse each species to its crown with a species-level tip vector $\Lambda_{ks}(z_c)$: the
   detection-aware species likelihood, propagated within the species by a star approximation. This
   cuts $N$ to ≈ 1.4×10⁵ representatives and handles the 10⁶-genome case.
3. Score clusters with $n_k\le2$ by a closed form: two tips and their MRCA path.

### 4.3 The HBB as a fast special case

If the tree is replaced by the rank tree with heights $h_\ell$, and pruning by its moment match (§1.4),
the result is `model_spec.md` §3.4 A. Implementation A therefore stays a valid quick baseline. The chain
should beat it by the amount predicted below.

---

## 5. Falsifiable predictions and a light check

### 5.1 Check on existing Chesterfield outputs (done, 5 CPU-s)

- Script: `/data/scratch/metapanclust_data/theory_phylo/rank_decay.py`; output: `rank_decay.out`.
- Input: `catalogue_surprise/inputs.npz`. "Seen" means I or P.
- For bin pairs grouped by deepest shared GTDB rank, the quantity is a completeness-corrected
  conditional,
  $\widehat{\Pr}(z_j=1\mid \text{seen}_i)=2\sum\text{both}/\sum(\text{seen}_i c_j+\text{seen}_j c_i)$.
  Families are stratified by the number of bins they are seen in (out of 93). All 93 bins; pair counts
  in brackets.

| Family seen in | domain only (3093) | phylum (266) | class (422) | order (91) | family (157) | genus (66) |
|---|---|---|---|---|---|---|
| 2-3 bins | 0.01 | 0.03 | 0.05 | 0.14 | 0.36 | 0.28 |
| 4-9 | 0.03 | 0.10 | 0.24 | 0.52 | 0.73 | 0.72 |
| 10-29 | 0.14 | 0.45 | 0.66 | 0.81 | 0.87 | 1.00 |
| ≥ 30 | 0.75 | 0.94 | 0.95 | 1.06 | 1.01 | 1.08 |

**Reading** (restricting to bins ≥ 50% complete gives the same pattern):

- **The informative rank moves with frequency, as rate heterogeneity predicts.**
  - Common families are already at ≈ 0.95 when the pair shares only a phylum or class (slow).
  - Rare families inform only within order or family (fast).
  - For the rarest stratum, $-\ln\Pr$ falls from 3.0 (class) to 2.0 (order) to 1.0 (family). With
    $d=2h$ from the RED heights of §2, this gives $\mu\approx 2.9, 2.7, 2.2$ per RED unit: close to a
    single exponential.
  - The middle strata give $\mu d$ ratios that fall with $d$, the log-convex signature of a rate
    mixture (§1.2).
- **Genus pairs are not more concordant than family pairs** (0.28 vs 0.36, 0.72 vs 0.73). On a tree,
  closer pairs must be at least as concordant. Two explanations remain, and both are §2.2 reasons:
  - within Chesterfield, genus-pair distances are not shorter than family-pair distances (deep marine
    genera);
  - the cross-genus construction filter removed the genus-specific clusters.

  Rank is not a good proxy for distance at the bottom of the hierarchy.
- **Caveats.**
  - The ≥ 30 stratum exceeds 1 because the completeness correction is crude.
  - Pairs are not independent: one bin contributes many pairs.
  - "Number of bins seen" is a noisy, sample-dependent proxy for rate.

### 5.2 Predictions to test next (existing data, no new tools beyond the container)

- **P1. Distance, not rank, governs concordance** (Chesterfield, then PanGBank).
  - Recompute the table above against a continuous distance: GTDB-Tk tree patristic distance from the
    classify-step tree, or AAI.
  - Prediction: within each frequency stratum, genus and family pairs fall on one decreasing curve in
    $d$.
  - Within one rank, the closer half of pairs is more concordant than the farther half, by a factor
    $\approx e^{-\mu\Delta d}$.
  - **Falsified if** within-rank concordance does not depend on $d$.
- **P2. Screening off in PanGBank B** (`worker.py genus` outputs already hold the needed splits).
  - Re-run the comparison at $m\in\{1,2\}$ as well as at 5, 10, 20.
  - Prediction: genus information beats `sp_eb` at $m=1$ (and likely 2) in the shared stratum and loses
    at $m\ge5$. The crossover $m^*$ is where the genus pseudo-count $\kappa(1-\bar p)/(m+\kappa)$ falls
    below the discrete-state error, about $m^*\approx 3$ at $\kappa=0.5$.
  - A discrete-state predictor should match or beat `sp_eb` in focal_only and congener_only at every
    $m$:
    - the species-crown state is inferred from the $m$ genomes;
    - then $P_{s1}(h_s)$ is applied;
    - with genus congeners entering only $\Pr(z_c)$.
  - **Falsified if** the chain-like predictor loses to `sp_eb` at $m\ge5$, or if genus helps
    increasingly with $m$.
- **P3. Rank $\kappa$ tracks depth ratios** (40-genus full run, plus families).
  - Per-genus fitted $\kappa_{\text{species}|\text{genus}}$ should correlate positively with
    $h_s/(h_g-h_s)$ from the GTDB tree (species-crown height over the genus crown-to-species gap).
  - The same holds for per-cluster $\kappa$ against within-species polymorphism, which proxies $\mu h_s$:
    clusters polymorphic in many species have $\kappa\gg1$, fixed clusters have $\kappa\approx$ the
    geometric value.
  - **Falsified if** $\kappa$ is unrelated to depth ratio for fixed clusters.
- **P4. Level-1 vs level-2 overlap.** The chain on a within-species ANI tree should capture most of the
  "Jaccard-neighbour" gain in `context_model.md` (0.04-0.06 bits) and part of the logistic gain. The
  rest of the logistic gain is co-inheritance of islands, which A2 excludes.

---

## 6. Implications for metapanclust

1. **Keep the HBB (implementation A) as the fast baseline, but fix its $\kappa$ by geometry.** Start
   from $\kappa_\ell=h_\ell/(h_{\mathrm{pa}}-h_\ell)$ computed from GTDB RED. Let per-cluster $\kappa$
   rise only through a turnover covariate. Do not fit a free global $\kappa_{\text{genus}}$ on
   catalogues filtered for cross-genus support: that filter biases it upward.
2. **Replace "shrink the species towards the genus" with "infer the species-crown state, then
   propagate".**
   - For PanGBank-scale species ($m\ge5$), species data dominate and genus enters only through the
     crown-state posterior.
   - For novel species or genera (most Chesterfield bins), the query attaches to a deep edge. Prediction
     then comes from the nearest sampled clades weighted by $e^{-\mu d}$, which is where taxonomy helps
     and where the 16% came from.
3. **Build level 1 as the gain/loss chain on the GTDB tree** with detection tip vectors (§3.1), 64
   rate×prevalence classes pooled by NPMLE (§3.3), family-level lineage multipliers, and Steiner-tree
   sparsity with absent-subtree caching (§4.2). The up-down pass also gives exact leave-one-out (and
   leave-ANI-group-out), which retires the §3.4 A leak and its per-bin refit.
4. **Report surprise with its rate class.** A present cluster with posterior mass on fast classes
   carries little phylogenetic surprise, whatever its prevalence. Absence of a slow cluster in a
   high-detection MAG is the informative event. The class posterior $\gamma_k$ is a principled
   core/shell/cloud label tied to turnover rather than to frequency thresholds.
5. **Use distances, not ranks, wherever possible.** Even inside implementation A, use skani ANI
   (within genus) and GTDB-Tk patristic distances (above it) to weight relatives by $e^{-\mu_r d}$. This
   is the cheapest way to get most of the chain's benefit.
6. **What level 1 cannot do** (unchanged): never-seen content, and correlated islands. Those stay with
   level 2 and with block scoring (`context_model.md`).

## Open points

- RED-to-height conversion and species-crown heights are rough. Compute them from the r226 tree before
  using §2's numbers quantitatively.
- HGT gains are not tree-like. Gains by close-relative transfer inflate concordance at short distance in
  a way that mimics low $\mu$. A "gain from a neighbour" (spatial) term would be a level-2 feature.
- The detection model's $\beta$ and $\lambda_k$ interact with rate estimation. Low-completeness MAG tips
  look like losses on terminal branches and inflate $l$. Fit detection parameters jointly (§3.2), or fit
  rates on isolate and high-quality tips first.
