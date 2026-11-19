# Surprise model: statistical specification (draft, 2026-11-18)

**What this is.** This is an implementable specification of the surprise model. It defines:

- **Level 1:** prevalence of each conserved gene cluster along the GTDB taxonomy, using a hierarchical
  beta-binomial with detection-corrected counts.
- **A detection model:** how a cluster that is present or absent shows up in a MAG, given completeness,
  contig ends, composition-driven binning loss, depth and contamination.
- **Level 2:** a context model, one masked logistic model per cluster with a level-1 offset.
- **A population version:** cluster frequency estimated from coverage ratios to the core.
- **Surprise scores:** per cluster, per bin, and per Mb.
- **A binning rule:** a Bayes factor between "one population" and "two populations" for merges and
  contig moves.

Two things are new relative to `first_tests.md` and `fragmentation.md`:

- The detection model has an explicit three-state observation (intact / partial / absent) and a
  false-presence term.
- Level 1 is fitted *through* the detection model, because most GTDB genomes in our clades are MAGs.

Related methods and what each contributes are in [related_work.md](related_work.md).

The prototype in `scripts/module_surprise.py` is a special case: one species level, Beta(1,1) prior, no
fragmentation, no false presence. Present scores -log2 p; absent scores -log2(1 - p c).

---

## 1. Notation

| Symbol | Meaning |
|---|---|
| k = 1..K | conserved gene clusters (the catalogue; each k is a set of gene families F_k, m_k = abs(F_k) ≥ 2, typical span s_k bp) |
| g | a reference genome (GTDB representative or member); j a query genome or bin; b a population (co-assembly) bin |
| t, pa(t), ch(t) | node of the GTDB taxonomy (domain → phylum → class → order → family → genus → species), its parent, its children |
| ℓ(t) | rank of node t |
| s(g) | species node of genome g |
| z_gk ∈ {0,1} | true presence of cluster k in genome g (at least τ_f of F_k co-located; default τ_f = 0.5, as in the PPanGGOLiN `completion` cut-off used so far) |
| O_gk ∈ {I, P, A} | observed state: **I**ntact (≥ τ_f of F_k on one contig, not touching a contig end); **P**artial (some families found, but truncated at a contig end or split across contigs); **A**bsent (no family of F_k in the genome) |
| n_gk | number of observed copies (for duplication terms) |
| θ_kt | prevalence of k in clade t: P(z = 1) for a random genome of clade t |
| κ_kℓ | concentration at rank ℓ: how tightly child-clade prevalences follow their parent |
| c_j, r_j | completeness and contamination of j (CheckM2; fractions) |
| L_j, G_j | binned length of j; estimated true genome size G_j = L_j / c_j |
| f_j(s) | P(an interval of span s lies inside one binned contig, given it lies on binned DNA), from j's contig lengths |
| λ_k | composition-driven binning loss of k (the contig carrying k is not binned with its genome) |
| ε_jk | false presence: P(k observed in j, given j's genome lacks k) |
| p̂_jk | predictive prevalence for j (level 1: p̂^(1); level 2: p̂^(2)) |
| s_jk | surprise in bits, -log2 P(O_jk given context) |
| i = 1..n | samples, for population bins |
| R_ik, R_iC | reads mapped to cluster k, and to the bin's core reference set C, in sample i |
| φ_k | fraction of cells of population b carrying k (cluster frequency; may exceed 1 for multi-copy) |

Logs are base 2 for surprise and natural logs for likelihoods. ln BF = ln2 × (difference in bits).

---

## 2. Observation (detection) model for one genome or MAG

The detection model is defined first because level 1 is fitted through it.

### 2.1 Probabilities given true presence (z = 1)

From `fragmentation.md`, which has the closed form and simulations:

```
P_intact(j,k) = c_j * (1 - β (1 - f_j(s_k)))                     β = 1 default (random breaks)
P_lost(j,k)   = (1 - c_j) f_j(s_k) + (1 - c_j)^2 (1 - f_j(s_k))   (approximation; see fragmentation.md)

f_j(s) = Σ_i max(ℓ_i - s + 1, 0) / Σ_i ℓ_i                       ℓ_i = binned contig lengths
```

**Composition loss.** Add the cluster-specific loss λ_k: the cluster's contig was not binned even though
the genome is "complete" by markers. Mise & Iwasaki 2022 show this for ribosomal-protein operons, and
Maguire et al. 2020 for genomic islands. Add a depth-dependent assembly term a_j: the probability that a
region at the bin's coverage C_j assembles at all.

```
P(I | z=1) = (1 - λ_k) a_j P_intact(j,k)
P(A | z=1) = λ_k + (1 - λ_k) [ a_j P_lost(j,k) + (1 - a_j) ]
P(P | z=1) = 1 - P(I | z=1) - P(A | z=1)
```

**Estimating λ_k.** Within GTDB species that have both isolates and MAGs, compare the
detection-corrected prevalence among MAGs with that among complete isolates:

```
λ̂_k = 1 - θ̂_k(MAGs) / [ θ̂_k(isolates) × mean_g P_intact-or-partial(g,k) ]      clipped to [0, 1)
```

Pool across species with a logistic regression on cluster covariates, so that unseen clusters get a
value:

```
logit λ_k = γ0 + γ1 |GC_k - GC_genome| + γ2 |TNF distance_k| + γ3 [MGE_k] + γ4 [rRNA/rp-adjacent_k]
```

**Assembly term.** a_j = 1 - exp(-(C_j / C_50)^h). Fit C_50 and h from danaSeq bins: the fraction of
level-1 core clusters (p̂ ≥ 0.95) observed as I or P, as a function of bin coverage. Default a_j = 1
when C_j ≥ 10×.

### 2.2 Probabilities given true absence (z = 0)

False presence comes from contaminating contigs and from annotation or clustering errors:

```
ε_jk = ε0 + (1 - exp(-r_j G_j ρ_k))       ρ_k = copies of k per bp of "donor" DNA
```

- **ρ_k.** Estimate it from the other bins of the same assembly: Σ_{j'≠j} n_{j'k} / Σ_{j'≠j} L_{j'}.
  This is the density of k in the contaminant pool. It handles the Majernik et al. 2025 case, where a
  contaminant gene tracks a different abundant taxon.
- **ε0.** A floor for annotation and family-assignment false positives. Default 1e-4; estimate from
  clusters "observed" in complete isolates of clades where the clade prevalence is 0.
- **Splitting ε_jk between I and P.** Use the same fragmentation logic on the donor contig lengths. A
  contaminant is usually a short contig, so most false presences are P.

```
P(I | z=0) = ε_jk × f_j(s_k) (1 - λ_k)     (contaminant contig long enough to hold the cluster)
P(P | z=0) = ε_jk - P(I | z=0)
P(A | z=0) = 1 - ε_jk
```

### 2.3 Marginal over z

For any predictive prevalence p (level 1 or level 2):

```
P(O = o | p) = p · P(o | z=1) + (1 - p) · P(o | z=0)                                    (D1)
```

This is linear in p, which is what makes level-1 fitting with detection tractable (§3.3).

Special cases that serve as checks:

- **Complete isolate** (c = 1, f = 1, λ = 0, ε → 0): P(I) = p and P(A) = 1 - p. Plain Bernoulli.
- **Prototype:** P(A) = 1 - p c.
- **`fragmentation.md` absence score, -log2(1 - p P_intact):** this treats every non-intact observation
  as absent. With the three states above, an observed A uses P(A) = p P(A|z=1) + (1-p)(1-ε), and
  P(A|z=1) = P_lost ≤ 1 - P_intact. So fragmentation.md's form is the conservative upper bound on P(A)
  (fewer bits). Use (D1) once λ and ε are estimated; keep the conservative form until then.
- **mOTUpan core likelihood:** p = 1, f = 1, λ = 0, ε = 0, so P(A) = 1 - c.

### 2.4 Copies

The duplication term is used in contamination and merge scoring. Given z = 1 in one genome:

```
P(n ≥ 2 | z=1) = δ_k
```

δ_k is the within-genome duplication rate of k in complete isolates of the clade, with Beta(0.5, 20)
smoothing. For level-1 core clusters (p̂ ≥ 0.95, δ_k < 0.02), extra copies are CheckM-style
contamination evidence.

---

## 3. Level 1: prevalence over the taxonomy (hierarchical beta-binomial)

### 3.1 Model

For each cluster k independently:

```
θ_k,root ~ Beta(a0, b0)                                       per domain; a0, b0 fitted across clusters
θ_k,t | θ_k,pa(t) ~ Beta(κ_kℓ θ_k,pa(t), κ_kℓ (1 - θ_k,pa(t)))   ℓ = ℓ(t)
z_gk | θ_k,s(g) ~ Bernoulli(θ_k,s(g))
O_gk | z_gk     ~ detection model of genome g (§2)
```

Properties of the parent-to-child step:

- **The mean is preserved:** E[θ_t | θ_pa] = θ_pa.
- **Variance:** Var[θ_t | θ_pa] = θ_pa (1 - θ_pa) / (κ + 1).

How to read κ:

- **κ_kℓ small (≲ 1):** child clades are near 0 or 1. The cluster is a clade marker, gained or lost once
  and then inherited (for example a photosystem).
- **κ_kℓ large (≳ 100):** every child clade has about the parent's prevalence. Typical of shell or
  cloud content spread by HGT.
- **Sharing κ across genera.** Henoch et al. 2026 report that a region's "propensity to vary is
  conserved across genera". That supports one κ_k per rank shared by all clades, rather than per clade.

**Why taxonomy ranks and not the tree.** GTDB ranks are RED-normalised, so "genus" means roughly the
same divergence everywhere. A rank-indexed κ is then interpretable and needs only 6 numbers per cluster
(and only 6 shared numbers for most clusters). A continuous-time gain/loss model on the GTDB tree is the
alternative (castor-style, or the phylogenetic occupancy model of Mattick et al. 2026); see open choice
O1.

### 3.2 Genome weights (clonal oversampling)

GTDB species are sampled very unevenly, and near-identical genomes are not independent draws. Weight
genome g by:

```
w_g = 1 / m_g^γ
```

- m_g is the size of g's 99.5% ANI (skani) cluster within its species.
- γ = 1 by default: each near-clonal group counts once.
- Weights enter as likelihood powers: Π_g P(O_g | θ)^{w_g}.

### 3.3 Species-level likelihood with detection

Collapsing z with (D1), the data for species s and cluster k give a likelihood in θ = θ_ks:

```
Λ_ks(θ) = Π_{g ∈ s} [ θ a_gk + (1 - θ) b_gk ]^{w_g}        a_gk = P(O_gk | z=1), b_gk = P(O_gk | z=0)
```

- For complete genomes with ε → 0 this reduces to θ^x (1 - θ)^{n - x}, the conjugate beta-binomial.
- For MAGs it is a polynomial in θ. It can be evaluated exactly on a grid, or approximated by a Beta
  with matched mean and variance.

### 3.4 Inference

Two implementations: **A** for production at catalogue scale, **B** for validation and for fitting κ.
They must agree on held-out log loss within Monte Carlo noise before A is trusted.

**B. Exact belief propagation on a grid.**

- Discretise θ on a logit-spaced grid u_1..u_G (G = 128, from 1e-4 to 1 - 1e-4). The tails matter
  because surprise lives there.
- Precompute one Beta transition kernel per (κ, rank): K_κ[i, p] = Beta(u_i; κ u_p, κ(1 - u_p)) Δu_i.
- Upward pass:

```
λ_t(u)    = Λ_kt(u) × Π_{c ∈ ch(t)} m_c(u)        (Λ_kt = 1 if no genome attaches directly to t)
m_t(u_p)  = Σ_i K_κ(ℓ(t))[i, p] λ_t(u_i)            message from t to its parent
```

- Downward pass:

```
post_root(u) ∝ Beta(u; a0, b0) λ_root(u)
post_t(u)    ∝ λ_t(u) Σ_p K[i(u), p] post_pa(u_p) / m_t(u_p)
```

- **Predictive prevalence** for a new genome placed at t: p̂_kt = Σ_u u post_t(u).
- **Marginal likelihood** of cluster k: Z_k = Σ_u Beta(u; a0, b0) λ_root(u), used to fit κ (§3.5).
- **Cost.** O(N_nodes G²) per cluster. Run it only on the subtree spanning genomes where k is observed
  (I or P), plus their siblings. All-absent subtrees produce identical messages that depend only on
  (n, rank) and can be cached.

**A. Sequential empirical-Bayes shrinkage** (default at scale; closed form, O(nodes) per cluster).

1. **Detection-corrected counts at each species** (one EM step from a starting p̂):

```
x_s = Σ_g w_g P(z_gk = 1 | O_gk, p̂) = Σ_g w_g p̂ a_gk / (p̂ a_gk + (1 - p̂) b_gk)
n_s = Σ_g w_g
```

2. **Top-down shrinkage.** Pooled subtree counts (X_t, N_t) are sums over species of x_s and n_s, with
   each species capped at n_max = 20 effective genomes, so that one huge species cannot set its genus's
   value.

```
p̂_root = (X_root + a0) / (N_root + a0 + b0)
p̂_t    = (X_t + κ_ℓ(t) p̂_pa(t)) / (N_t + κ_ℓ(t))                    for t below root, in BFS order
```

   For a species with no genomes, p̂_s = p̂_genus. For a novel genus, p̂ = p̂_family, and so on.

3. **Iterate** steps 1-2 two or three times. The posterior z probabilities depend on p̂ through the
   detection model.

**Leave-one-out.** When scoring a genome that is itself in the reference set, subtract its (and its
99.5% ANI group's) contribution from x_s and n_s before computing p̂. Otherwise our own bins inform their
own expectation. The PanGBank GTDB_all pangenomes contain GTDB MAGs and so need this too.

### 3.5 Concentration parameters

**Global, per rank** (always). Choose κ_ℓ to maximise the summed marginal likelihood over clusters:
Σ_k ln Z_k(κ) with B, or with A plus a beta-binomial approximation. Search a log grid
κ ∈ {0.1, 0.2, 0.5, …, 2000}, coordinate-wise from species upward.

- **Species level, starting point.** With complete genomes only, κ_species maximises
  Σ_{s,k} ln BetaBin(x_sk | n_sk, κ p̂_genus,k, κ(1 - p̂_genus,k)).
- **Bias in genus p̂.** p̂_genus is itself estimated from those species. Use leave-species-out p̂_genus
  so that κ is not biased upward.

**Per cluster** (when there is data). Put a lognormal hyperprior on κ:

```
ln κ_kℓ = μ_ℓ + βᵀ x_k + e_kℓ,     e_kℓ ~ N(0, τ_ℓ²)
x_k = (MGE flag, ln m_k, logit of root prevalence, fraction of the cluster's families in persistent partitions)
```

- Fit μ, β and τ by maximising the marginal likelihood. Then take each cluster's posterior mode of
  ln κ_kℓ, given its own Z_k(κ) on the grid.
- Clusters seen in fewer than 5 clades at rank ℓ keep the regression prediction.

**Hyperprior at the root.** a0 and b0 come from method of moments on the domain-level p̂ across
clusters. Expect a U shape (a0, b0 < 1).

### 3.6 Placing a query bin

- GTDB-Tk gives the lineage of j. If the species is a GTDB species, use p̂_k,s. If the bin is a novel
  species, use the deepest assigned rank.
- Our own species clusters from 95% ANI dereplication across danaSeq assemblies can be added as new
  species nodes under their genus, with their own bins as data. Then the second Chesterfield bin informs
  the first, but not itself (leave-one-out, §3.4).
- **Placement uncertainty.** If GTDB-Tk reports a placement with low RED or ANI support, mix:
  p̂ = Σ_t P(t) p̂_kt over candidate nodes (open choice O4).

### Pseudocode (implementation A)

```python
def fit_level1(clusters, taxonomy, genomes, kappa, a0, b0, n_iter=3):
    # genomes: dict g -> {species, weight w_g, detection params (c, f(.), C, r, ...)}
    # obs[g][k] in {'I','P','A'}; stored sparsely: only I/P listed, A implied
    p_hat = {}                                   # (k, node) -> predictive prevalence
    for k in clusters:
        p_g = {g: 0.5 for g in genomes}          # starting prevalence per genome
        for it in range(n_iter):
            X, N = defaultdict(float), defaultdict(float)
            for g in genomes:
                a = P_obs_given_present(g, k)    # §2.1, for the observed state of k in g
                b = P_obs_given_absent(g, k)     # §2.2
                post = p_g[g] * a / (p_g[g] * a + (1 - p_g[g]) * b)
                s = genomes[g].species
                X[s] += genomes[g].w * post
                N[s] += genomes[g].w
            cap_species(X, N, n_max=20)
            pool_up(X, N, taxonomy)              # subtree sums
            for t in taxonomy.bfs():
                prior = (a0, b0) if t.is_root else kappa[k, t.rank] * p_hat[k, t.parent]
                p_hat[k, t] = (X[t] + prior_mean_mass(prior)) / (N[t] + prior_total(prior))
            p_g = {g: p_hat[k, genomes[g].species] for g in genomes}
    return p_hat
```

Two notes on running this for real:

- **Sparsity.** Loop only over genomes where k is I or P. Handle all-A genomes per species in closed
  form: A depends only on each genome's P(A|z), which is summed per species once per cluster, not per
  genome.
- **Leave-one-out at scoring time.** Subtract the scored genome's terms from X and N along its path and
  recompute p̂ along that path only, which costs O(depth).

---

## 4. Surprise

### 4.1 Per cluster

```
P_jk(o) = p̂_jk P(o | z=1, j) + (1 - p̂_jk) P(o | z=0, j)           o ∈ {I, P, A}
s_jk    = -log2 P_jk(O_jk)                                          bits
h_jk    = -Σ_o P_jk(o) log2 P_jk(o)                                 expected surprise (entropy)
v_jk    = Σ_o P_jk(o) (log2 P_jk(o))² - h_jk²                       variance of surprise
e_jk    = s_jk - h_jk                                               excess surprise
```

**Direction:**

- **Surprisingly present:** O ∈ {I, P} and e > 0.
- **Surprisingly absent:** O = A and e > 0.

**Scoring P observations.** P is scored with its own probability. A cluster truncated at a contig edge
is not an absence, and its surprise is small when p̂ is moderate. When p̂ ≈ 0 it is a strong
"present" surprise, discounted only by ε. That is correct: half of a foreign operon on a short contig is
contamination evidence.

**Which clusters to sum over.** Clusters with p̂ < 1e-4 that are not observed contribute about 0 bits,
so only clusters observed in j, or with p̂ ≥ 1e-4 at j's clade, need computing.

### 4.2 Per-bin totals

```
S_j+ = Σ_{k: O ∈ {I,P}} e_jk          presence excess (foreign or unexpected content)
S_j- = Σ_{k: O = A}     e_jk          absence excess (missing expected content)
Z_j  = Σ_k e_jk / sqrt(Σ_k v_jk)      standardised total (independence approximation)
```

**Normalisation:**

- **Presence per Mb of bin:** D_j+ = S_j+ / (L_j / 1e6). Contamination adds sequence and presence bits
  in proportion, so bits per Mb is comparable across bin sizes.
- **Absence per expected Mb:** D_j- = S_j- / (Σ_k p̂_jk s_k / 1e6). The denominator is the amount of
  conserved-cluster sequence the clade leads us to expect. Missing content scales with what was
  expected, not with what was binned.
- **Per contig, for localisation:** d_c = Σ_{k on contig c} e_jk / (len(c) / 1e6), plus a GUNC-style
  check of whether presence surprise is concentrated on few contigs (open choice O7).

**Dependence.** Clusters are not independent. Overlapping clusters (shared families) and clusters that
co-occur through co-selection make Z_j overconfident. Mitigations, in order:

1. Score a non-overlapping cluster set. Resolve nested and overlapping clusters greedily by m_k, keeping
   the larger.
2. Use level-2 conditionals (§5) in the sum. That is a pseudo-likelihood, not a joint: label it
   "pseudo-surprise".
3. Calibrate the null distribution of Z_j empirically on clean complete genomes fragmented to MAG-like
   contig sets (§9, E4), and report it as a quantile or empirical P value, not as N(0,1).

### 4.3 MGE clusters

Clusters on geNomad plasmid or virus contigs, or with the MGE flag in x_k, are reported in a separate
column and excluded from S_j+ and S_j-. Horizontal transfer is supposed to be surprising.

---

## 5. Level 2: context model (masked logistic with clade offset)

**Goal.** P(z_gk = 1 | clade, all other clusters of g). The model has to learn only co-variation *beyond*
the clade expectation; otherwise inheritance is counted twice.

### 5.1 Model

One model per cluster k:

```
logit P(z_gk = 1 | ·) = logit p̂^(1)_{k,s(g)}^{(-g)}  +  b_k  +  Σ_{k' ∈ N_k} w_kk' (ẑ_gk' - p̂^(1)_{k',s(g)})
```

- **Offset.** The level-1 leave-one-out prediction. It carries the clade, and it carries all inheritance
  at taxonomic resolution.
- **Features.** Residuals ẑ - p̂^(1). Here ẑ_gk' = P(z_gk' = 1 | O_gk') from level 1 plus detection: soft
  labels, so MAG incompleteness does not masquerade as context. A cluster that is simply typical of the
  clade contributes zero.
- **Candidate set N_k.** Up to M = 200 clusters, chosen by phylogenetically corrected association with k:
  1. Restrict to clades where k is variable (0.05 < p̂^(1) < 0.95).
  2. Rank k' by |correlation of residuals|.
  3. Optionally keep only pairs that pass Goldfinder's simulation null (simultaneous score, FDR 0.05) in
     at least one species where both are variable.
  4. Exclude clusters that share families with k, and clusters physically adjacent to k (within 5 kb) in
     more than 50% of genomes carrying both. These are the same locus, and should be merged in the
     catalogue rather than used as context.
- **Within-species population structure** (open choice O5). Add the top 5 principal coordinates of the
  species' core-genome distance matrix as covariates in species with ≥ 30 genomes. Or rely on the
  Goldfinder filter alone.

### 5.2 Fitting

- **Training genomes.** Complete isolates plus MAGs with c ≥ 0.95 and r ≤ 0.02. The target is ẑ_gk, a
  soft label, using cross-entropy with soft targets.
- **Regularisation.** Elastic net (α = 0.9), with the strength chosen by grouped cross-validation where
  folds are genera. Genus folds prevent the context from memorising clade identity that level 1 already
  has.
- **Masking** is implicit: k is never in its own feature set.
- **Scope.** One global model per cluster with clade offsets (default), or per-family models for clusters
  whose role changes across families (open choice O6).
- **Keep only models that help.** Keep the level-2 model for k only if its grouped-CV log loss beats
  level 1 by at least 0.01 bits per genome. Otherwise p̂^(2) = p̂^(1).

### 5.3 Use on a bin

- p̂^(2)_jk replaces p̂^(1)_jk in §4. Report both, and the context gain s^(1) - s^(2).
- **Contaminated context.** If j is contaminated, its features are too. Compute ẑ for the features only
  from *anchor contigs*: contigs carrying at least one cluster with p̂^(1) ≥ 0.9 observed I, or the
  largest contigs covering 50% of L_j.
- **Iterate after a move.** After a contig move (§7), recompute the context features for that bin.

### Pseudocode

```python
def fit_level2(k, genomes, p1, zsoft, folds_by_genus, M=200):
    var_clades = [t for t in clades if 0.05 < p1[k, t] < 0.95]
    G = [g for g in genomes if species(g) in var_clades and is_training(g)]
    resid = {kk: [zsoft[g][kk] - p1[kk, species(g)] for g in G] for kk in candidate_pool(k)}
    Nk = top_by_abs_corr(resid, target=[zsoft[g][k] - p1[k, species(g)] for g in G], M=M,
                         exclude=same_locus(k) | shares_families(k))
    Nk = goldfinder_filter(k, Nk) if USE_GOLDFINDER else Nk
    X = matrix(resid[kk] for kk in Nk).T
    offset = [logit(p1_loo(k, g)) for g in G]
    y = [zsoft[g][k] for g in G]
    model = elastic_net_logistic(X, y, offset, alpha=0.9, cv=folds_by_genus)
    return model if model.cv_logloss < level1_cv_logloss(k) - 0.01 else None
```

---

## 6. Population level: cluster frequency from coverage ratios

A co-assembly bin b is a population. The question is no longer "is k present?" but "what fraction φ_k
of cells carry k?". This is answered across samples from read counts relative to the bin's core.

### 6.1 Core reference

C_b is the set of clusters in b with p̂^(1) ≥ 0.95, observed I, single copy (n = 1), and not on MGE
contigs. Reads mapping to C_b in sample i: R_iC, over effective length L_C. The effective mappable
length excludes the 150 bp at each contig end and low-MAPQ positions.

### 6.2 Likelihood

Conditioning on the total reads of cluster plus core in each sample removes depth and library size:

```
N_ik = R_ik + R_iC
R_ik | N_ik, φ_k ~ BetaBinomial(N_ik, mean π_k, concentration ω_b)
π_k = φ_k L_k / (φ_k L_k + L_C)
```

- **ω_b**, the overdispersion from mapping, GC and coverage noise, is calibrated inside the bin. Split C_b
  at random into two halves. Their read ratio has known φ = 1, so ω_b is fitted by maximum likelihood
  over samples and 20 random splits.
- **The beta-binomial on read counts** is the answer to "beta-binomial on read counts?": yes. The
  conditional form is the same one MIDAS2 uses implicitly (copy number = gene coverage / marker
  coverage). Here it is made into a likelihood with overdispersion, pooled over samples.
- **Cross-mapping.** Reads from a related population inflate R_ik. Use MAPQ ≥ 20 and an identity filter
  of at least 97% (or the bin's popANI minus 2 points). Add a background term ν_ik if needed (open
  choice O8).

**Constancy across samples.** Two hypotheses:

- **H_const:** φ_ik = φ_k in all samples.
- **H_var:** φ_ik ~ Beta(φ_k ψ, (1 - φ_k) ψ) per sample. Strain composition shifts, or the contig is not
  part of this population.

The log Bayes factor between them comes from 1-D and 2-D grid integration. It feeds the binning rule
(§7).

### 6.3 Prior on φ_k from levels 1 and 2

A population is a mixture of a few strains drawn from the species. With an effective strain number N_eff
and species-level prevalence p = p̂^(1 or 2)_{k,s(b)}:

```
φ_k ~ π0 δ(0) + π1 δ(1) + (1 - π0 - π1) Beta(p η, (1 - p) η)
π0 = (1 - p)^N_eff         (no strain carries it)
π1 = p^N_eff               (all strains carry it)
```

- **N_eff:** from inStrain nucleotide diversity of the bin's core, mapped to a strain count. Default 3
  when unknown (open choice O9). Default η = 2.
- **Multi-copy clusters** (δ_k > 0.1 in references) get support on (0, 3]. Everything else is capped at
  1.2, to allow for noise.

### 6.4 Posterior, states and detection limit

- Compute the posterior of φ_k on a grid: {0} ∪ 200 points on (0, 1.2] ∪ {1}.
- Define three states:
  - **absent:** φ < 0.05
  - **partial:** 0.05 ≤ φ < 0.9
  - **fixed:** φ ≥ 0.9

```
P_prior(state), P_post(state)
s_bk = Σ_state P_post(state) (-log2 P_prior(state))          expected surprise under the posterior
```

- **Reporting.** Report the MAP state when P_post(MAP) ≥ 0.9. Otherwise label it "unresolved". A
  low-depth bin gives a posterior close to the prior and is labelled unresolved, so depth never
  manufactures surprise.
- **Detection limit.** The smallest φ distinguishable from 0 at 95%. Roughly, with total core reads
  R_C = Σ_i R_iC, a cluster of length L_k gets an expected φ L_k R_C / L_C reads, and
  P(0 reads) = exp(-that). So:

```
φ_min ≈ 3 L_C / (L_k R_C)
```

  Report φ_min per bin. "Absent" calls require φ_min < 0.05.

**Clusters missing from the co-assembly.** A cluster absent from the assembly may still be in the reads
at low φ, below the assembly threshold. To test it, recruit reads to reference sequences of k from GTDB
genomes of the species (MIDAS-style), at an identity matched to the bin's ANI to those genomes, and run
the same model. Without that recruitment, an assembly absence is only P(A) from §2 with a_j derived from
φ C̄_b.

### 6.5 Population-level reading of the first test

For the UBA9145 flagellar modules (36 modules, none in 52 public genomes):

- The check above estimates φ for each module in each sample, and asks whether it is constant (H_const)
  and close to 1.
- **φ ≈ 1 and constant:** a motile local population. The surprise is real, biology.
- **φ varies, tracking another bin's coverage:** contamination. That is the Majernik et al. gene-to-species
  correlation in likelihood form.

---

## 7. Binning decision rule

**Question.** For a proposed edit, a merge of bins A and B or a move of contig set X from bin A to B
(or to unbinned), compare:

- **H1:** the result is one population.
- **H2:** they are separate (two populations, or X is not part of A).

```
ln BF = ln P(data | H1) - ln P(data | H2) = ℓ_cov + ℓ_content [+ ℓ_comp]
accept edit if ln BF + ln prior_odds > τ
```

### 7.1 Coverage term ℓ_cov

**Merge.** Let R_iA and R_iB be the reads on the cores of A and B, with core lengths L_A and L_B.

- **Under H1:** both cores are the same population's core, so every sample has the same ratio:

```
R_iA | R_iA + R_iB ~ BetaBin(N_i, π = L_A / (L_A + L_B), ω)
```

  If B has no core of its own (an accessory-only fragment), use π = φ L_B / (φ L_B + L_A) with φ
  integrated over its prior (§6.3).

- **Under H2:** π_i ~ Beta(π̄ ψ, (1 - π̄) ψ), free per sample, with π̄ and ψ integrated on a grid. This
  is the null that the two populations' abundances vary independently.

- ℓ_cov = ln marginal(H1) - ln marginal(H2).

**Move of X into A.** H1 uses the §6 likelihood for X with the cluster-informed prior on φ_X. H2 is the
same likelihood with X's own free per-sample profile (H_var with ψ → small), or with the best other
bin's core as reference.

**Accessory contigs.** This term implements "accept accessory contigs whose coverage follows the core at
a lower ratio". A contig at a constant φ ≈ 0.4 of A's core across 10 samples has a high ℓ_cov under H1
even though its raw coverage profile differs from A's.

### 7.2 Content term ℓ_content

```
ln P(content | H1) = Σ_k ln P_{A∪B}(O_{A∪B,k})   with p̂ from the union's placement, detection for the union,
                     + Σ_{k: n_{A∪B,k} ≥ 2} ln δ'_k   (duplication: δ'_k = δ_k for one genome)
ln P(content | H2) = Σ_k [ ln P_A(O_Ak) + ln P_B(O_Bk) ]   each with its own placement and detection
ℓ_content = ln2 × ( [S_A + S_B] - S_{A∪B} )         S = Σ_k s (bits), duplication terms included
```

- The union needs its own c and r. Run CheckM2 on the union, which takes seconds per bin.
- **Complementary bins** (A has half the core, B the other half, no duplicates) give S_{A∪B} < S_A + S_B,
  because absences are resolved. Merge is favoured.
- **Overlapping bins** (both carry the same core clusters) pay ln δ_k per duplicated core cluster.
  Merge is penalised. This is CheckM contamination as a likelihood.
- **Different placements** (A and B in different genera) carry clade-specific content that is
  surprising in the union's clade. Merge is penalised.

**Moves.** ℓ_content = ln2 × [S_A + S_rest(X) - S_{A+X}]. S_rest(X) is X scored as its own fragment
under the clade of the bin it came from (or "unbinned", scored at domain level).

**Composition term ℓ_comp** (optional; open choice O10). A tetranucleotide or GC likelihood ratio.
Binners have already used composition, so adding it double-counts evidence. Leave it out by default and
use it only as a tie-breaker.

### 7.3 Prior odds and threshold

- **prior_odds:** from binner agreement. With DAS Tool inputs, (number of binners co-assigning A and B
  + 0.5) / (number separating them + 0.5). Default 1.
- **τ:** calibrated to a target false-merge rate (FDR 5%) on CAMI II strain-madness and on simulated
  mixtures (§9, E6). Do not use τ = 0 (the raw Bayes rule) until the detection model is calibrated:
  overconfident content terms would dominate.

### Pseudocode: greedy refinement

```python
def refine(bins, contigs, samples, max_rounds=5):
    for rnd in range(max_rounds):
        proposals = []
        # merges: only pairs that share a species placement or are adjacent in the assembly graph
        for A, B in candidate_pairs(bins):
            lbf = ell_cov_merge(A, B) + ell_content(A, B) + log(prior_odds(A, B))
            proposals.append((lbf, 'merge', A, B))
        # moves: contigs whose presence-surprise density is high in their bin,
        # plus unbinned contigs carrying clusters with p̂ ≥ 0.2 in some bin's clade
        for X in suspicious_or_unbinned(contigs, bins):
            for B in target_bins(X):
                lbf = ell_cov_move(X, B) + ell_content_move(X, src(X), B)
                proposals.append((lbf, 'move', X, B))
        accepted = [p for p in sorted(proposals, reverse=True) if p[0] > TAU]
        apply_non_conflicting(accepted, bins)    # each bin edited at most once per round
        recompute(bins)                          # CheckM2, placement, level-1/2 p̂, surprise
        if not accepted: break
    return bins
```

---

## 8. Practical defaults

| Parameter | Default | Source or estimation |
|---|---|---|
| τ_f (fraction of families for "present") | 0.5 | PPanGGOLiN completion convention; sensitivity 0.33-0.8 |
| β (fragmentation coupling) | 1 (conservative) | fragmentation.md; 0.25 suggested for conserved clusters |
| ε0 | 1e-4 | complete isolates in clades with zero prevalence |
| γ (genome weighting) | 1 at 99.5% ANI | sensitivity 0, 0.5 |
| n_max per species in pooling | 20 | sensitivity |
| κ grid | 0.1-2000, log-spaced | marginal likelihood |
| θ grid (implementation B) | 128, logit-spaced, 1e-4 .. 1 - 1e-4 | |
| p̂ floor and ceiling | [1e-4, 1 - 1e-4] | caps single-cluster surprise at about 13.3 bits |
| Level-2 M, α | 200, 0.9 | grouped CV by genus |
| Population state cuts | 0.05 / 0.9 | |
| N_eff, η | 3, 2 | inStrain microdiversity when available |
| MAPQ, identity for read counts | 20, 97% | |
| Merge or move FDR | 5% | CAMI II strain-madness and simulations |

---

## 9. Calibration and evaluation plan

Each experiment names its baseline and the number that decides it. **Leakage rule throughout:** a genome
or bin being scored never contributes to the level-1 counts or level-2 training used to score it. The
same applies to its 99.5% ANI group and, for held-out-clade tests, its whole clade.

**E1. Level-1 prevalence, held-out clades**

- **Design.** Within GTDB r226 families with at least 3 genera and at least 20 species, hold out:
  - (a) single genomes (leave-one-out);
  - (b) whole species: predict from the genus;
  - (c) whole genera: predict from the family.
- **Metric.** Held-out log loss in bits per cluster per genome, on complete isolates only, so that
  detection is not confounded. Also Brier score and reliability diagrams in 10 bins of p̂.
- **Baselines:**
  - species-only Beta(1,1), which is the prototype;
  - no-taxonomy global prevalence;
  - CheckM-style 0/1 at the nearest node with ≥ 0.97;
  - castor Mk hidden-state prediction on the GTDB tree;
  - the Mattick et al. phylogenetic occupancy model, once its code scales.
- **Decision.** Keep implementation A if it is within 2% of B's log loss. Adopt the occupancy model if it
  beats B by more than 5% on (b) and (c).

**E2. Detection model**

- **Fragmentation.** Already done (fragmentation.md). Re-check P_intact calibration for β = 0.25 versus
  β = 1.
- **λ_k.** Estimate in species that have both isolates and MAGs. Validate the regression on held-out
  species: predicted versus observed MAG detection of ribosomal-protein and rRNA-adjacent clusters (the
  Mise & Iwasaki case).
- **ε.** Contaminate complete genomes in silico with 1-20% of donor contigs (same genus, family, phylum,
  or same sample bins) and check that predicted false-presence rates match.

**E3. Level 2**

- **Design.** Mask clusters in held-out complete genomes, with genus-grouped folds.
- **Metric.** Log-loss gain over level 1, overall and for variable clusters (0.05 < p̂ < 0.95).
- **Baselines:**
  - a Beavan-style random forest on raw presence, without offset;
  - a logistic model without offset, to show the leakage of clade into context;
  - a linear probe on Bacformer genome embeddings, if feasible.
- **Decision.** Use level 2 only for clusters where it wins (§5.2).

**E4. Bin-level surprise as a quality measure**

- **Simulation.** Fragment complete genomes to danaSeq-like contig length distributions. Subsample to
  c ∈ {0.5, 0.7, 0.9, 1.0}. Add contamination r ∈ {0, 0.02, 0.05, 0.1, 0.2} from same-genus, same-family
  and other-phylum donors. This mirrors CheckM's simulation design.
- **Metrics:**
  - null calibration of Z_j and D_j+ on clean genomes (empirical quantiles stable across c);
  - AUROC for detecting contamination at each r and donor distance, versus CheckM2 contamination and
    GUNC CSS;
  - contig-level AUROC of d_c for identifying the donor contigs, versus MAGpurify.
- **Expected win:** same-genus contamination, where GUNC is weak by design.

**E5. Population frequencies**

- **Simulation.** CAMISIM mixtures of 2-5 strains per species, with known gene content and with 5-20
  samples of varying strain proportions. Co-assemble with danaSeq.
- **Metrics:**
  - RMSE of φ_k;
  - state accuracy (absent/partial/fixed);
  - calibration of φ_min (the claimed detection limit versus the empirical one).
- **Baselines.** MIDAS2 copy number ≥ 0.35; PanPhlAn 0.25/0.5.
- **Real-data check.** The UBA9145 flagellar modules (§6.5).

**E6. Binning**

- **Datasets:**
  - CAMI II marine and strain-madness, which have a gold standard;
  - CAMI III toy human gut, the 2026 release with 10 individuals × 2 samples;
  - later, the CAMI III gut challenge data once the gold standard is public;
  - Chesterfield co-assemblies for qualitative checks.
- **Metrics:**
  - AMBER completeness, purity and ARI;
  - number of near-complete bins;
  - **recall of accessory and GI contigs** in the correct bin (the Maguire et al. failure mode);
  - false-merge rate at the chosen τ.
- **Baselines.** DAS Tool output (our current input), DAS Tool + MAGpurify, and GUNC-filtered bins.

**E7. Catalogue sanity** (the units themselves)

- **Operons.** Recall of EcoCyc transcription-unit pairs (in use: `scripts/eval_pairs.py`).
- **Defence systems.** PADLOC/DefenseFinder system recovery, at least the 95% Spacedust reports.
- **Predicted operons.** Recall of DOOR/ProOpDB operons in clades with predictions.
- **Biological checks.** p̂ patterns for known clade markers (photosynthesis clusters in Cyanobacteria
  versus heterotrophs; flagellar operons) should show small κ at the right rank.

---

## 10. Open choices

| ID | Choice | Options | Lean |
|---|---|---|---|
| O1 | Phylogenetic structure for level 1 | rank-indexed beta shrinkage (here) / gain-loss Mk on GTDB tree / phylogenetic occupancy model | ranks first; compare in E1 |
| O2 | Cluster presence definition | τ_f of families co-located / any family / weighted by family conservation | τ_f = 0.5; test 0.33-0.8 |
| O3 | Cluster catalogue source | Spacedust loci grouped by GFRR / PPanGGOLiN modules + persistent runs / both | both; modules alone miss core operons |
| O4 | Placement uncertainty | point placement / mixture over candidate nodes | point first |
| O5 | Within-species structure in level 2 | core-genome PCs / Goldfinder filter / phylogenetic random effect | Goldfinder filter + PCs |
| O6 | Level-2 scope | one global model per cluster / per-family models | global with offsets |
| O7 | Contig-level surprise test | per-contig density threshold / GUNC-style entropy with contig-size null | entropy with null |
| O8 | Cross-mapping in population counts | identity filter only / explicit background ν_ik from related bins | filter only first |
| O9 | N_eff for the φ prior | fixed 3 / from inStrain microdiversity / learned per species | fixed, then inStrain |
| O10 | Composition term in the binning BF | none / TNF likelihood ratio | none by default |
| O11 | Bin totals under dependence | non-overlapping set + empirical null / pseudo-likelihood with level 2 / learned effective number of clusters | first option |
| O12 | MGE handling | exclude from totals and report separately / own prior with no taxonomic shrinkage | exclude and report |
| O13 | Reference prevalence source | PanGBank GTDB_all (MAGs included) / GTDB_refseq (isolates) / own PPanGGOLiN runs | PanGBank GTDB_all through the detection model, with GTDB_refseq to estimate λ |

---

## 11. What would falsify the approach

- **E1:** level 1 does not beat species-only Beta(1,1) on held-out species. The taxonomy carries no
  usable signal for cluster content.
- **E4:** at matched contamination, Z_j or D_j+ never beats CheckM2 contamination plus GUNC CSS.
  Content surprise would then add nothing to quality control.
- **E6:** the content term in the binning rule makes merges worse than the coverage term alone.
