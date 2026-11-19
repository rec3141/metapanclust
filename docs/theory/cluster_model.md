# A model for conserved gene clusters: phylogenetic null, generative units, estimator (2026-11-18)

Goal: replace the heuristic catalogue pipeline (`experiments/catalogue.md`: cross-genus support counts,
dip splitting at 0.34, Leiden CPM at 0.1/0.4) with one definition of a *conserved gene cluster* and one
estimator. Each heuristic should then be a special case or an approximation of it. Sections 1-2 are
derivations. Section 3 is the algorithm. Section 4 lists predictions; two were checked on the existing
`phylo_scaling` matrices.

## Summary

- **Definition.** A conserved gene cluster is a set of gene families whose adjacencies are broken by
  rearrangement *more slowly than neutral adjacencies* along the phylogeny. Each adjacency has a
  retention exponent $\phi<1$; the neutral case is $\phi=1$. Under this model, if a neutral adjacency
  survives with probability $P_0$ between two genomes, a selected one survives with probability
  $P_0^{\phi}$.
- **Null.** A neutral adjacency survives a path of length $t$ with probability $e^{-\rho t}$. The log
  likelihood ratio (LLR) carried by one conserved observation is $(1-\phi)\rho t$. It goes to 0 as
  $t \to 0$, so close genomes return the base rate, which is what `phylo_scaling` saw.
- **Checked on existing data.** For every panel genome $k$, set
  $x_k = -\ln P(\text{conserved}\mid\text{different TU})$. This is a self-calibrated divergence, no ANI
  needed. The model predicts
  $\ln\frac{P(c\mid\text{same TU})}{P(c\mid\text{diff TU})} = (1-\phi)\,x_k$, a line through the origin.
  - E. coli, 66 genomes: slope 0.469, intercept -0.018, r = 0.977.
  - B. subtilis, 90 genomes: slope 0.471, intercept -0.062, r = 0.967.
  - Both organisms give $\phi \approx 0.53$: operon adjacencies are lost at about half the neutral
    rate.
- **Multiple genomes.** The sufficient statistic is the **phylogenetic diversity** (total branch length)
  of the subtree over which the adjacency is retained, minus a penalty per inferred break. Raw genome
  counts are the wrong statistic. "Distinct partner genera" in `02_loci.py` is a crude version of it.
- **Recall grows about as $\log n$.** Two reasons: per-pair conservation probabilities are spread
  roughly log-uniformly, and PD grows logarithmically in the number of tips sampled from a tree. A
  log-uniform fit on $n \le 10$ reproduces the E. coli curve at $n = 20, 30$ (0.342 and 0.382 against
  0.36 and 0.40 observed).
- **Leiden is an approximation.** CPM-Leiden on content similarity is correlation clustering with a
  linearised binomial LLR. The resolution $\gamma$ is the similarity at which "same unit" and "different
  units" are equally likely. It depends on how often families match by chance (paralogue density):
  $\gamma \approx 0.25$ for typical families and $\approx 0.5$ for ABC-like hubs. That is why the
  catalogue needs both 0.1 and 0.4.

---

## 1. A phylogeny-aware null for neighbourhood conservation

### 1.1 One adjacency, two genomes

Take two genes $i,j$ adjacent in query genome $A$ (same strand, intergenic gap $g$). Genome $B$ is
separated from $A$ by a path of total branch length $t$ through their MRCA, $t = 2\tau$ for an
ultrametric tree. Rearrangements (inversions, transpositions, insertions of mobile elements, deletions)
place breakpoints as a Poisson process along the chromosome. An adjacency is broken at rate $\rho$ per
unit branch length. This is the standard breakpoint model (Nadeau-Taylor; the DCJ/breakpoint distance
literature). Selection for co-localisation lowers the rate at which breaks are *fixed* at that
adjacency to $\phi\rho$, with $0 \le \phi \le 1$. For a strongly co-regulated operon, $\phi = 1 - s$ is
the fraction of breakpoints that are tolerated.

Three more processes matter:

- **Homologue retention** $h_{ij}(t)$: both families must still be present and detectable in $B$. It
  decays roughly like $e^{-\eta t}$, with family-specific $\eta$.
- **Detection** $\varepsilon$: Spacedust or MMseqs2 misses or mis-segments some clusters even at
  $t = 0$. Eco_Sakai conserves only 81% of different-TU pairs, so $1-\varepsilon \approx 0.81$ for
  "co-membership in one Spacedust cluster".
- **Convergent re-formation** $c_0$: the adjacency is re-created independently, by HGT of the whole
  cluster or by chance placement. This is small per pair but sets a floor.

Then

$$
P(\text{cons}\mid t, H_0) = (1-\varepsilon)\,h_{ij}(t)\left[e^{-\rho t} + c_0\right], \qquad
P(\text{cons}\mid t, H_1) = (1-\varepsilon)\,h_{ij}(t)\left[e^{-\phi\rho t} + c_0\right].
$$

$H_0$: the adjacency is neutral (the two genes merely happen to be neighbours). $H_1$: selected
co-localisation (same transcription unit, or a conserved superoperon junction).

**Why close genomes carry no evidence.** If $B$ conserves the adjacency, the likelihood ratio for
$c_0 \ll e^{-\rho t}$ is

$$
\Lambda_+(t) = \frac{e^{-\phi\rho t}}{e^{-\rho t}} = e^{(1-\phi)\rho t}, \qquad
\ln\Lambda_+(t) = (1-\phi)\rho t \xrightarrow{t\to 0} 0 .
$$

With prior $\pi$ (the base rate of same-TU pairs), the precision of the call "conserved, so same TU"
from one genome is

$$
\mathrm{prec}(t) = \frac{\pi\, e^{-\phi\rho t}}{\pi\, e^{-\phi\rho t} + (1-\pi)\, e^{-\rho t}}
\xrightarrow{t\to 0} \pi .
$$

So precision equals the base rate at $t=0$. Measured: E. coli level 0, 0.62 against a base rate of
0.612; B. subtilis levels 0-1, 0.738-0.754 against 0.736. Conservation by descent means that every
adjacency, neutral or not, is still present in a close relative. Precision climbs towards 1 as $t$
grows. It is capped by $c_0$ (HGT'd or convergent clusters) and by the $\varepsilon$/annotation noise
in the labels:

$$
\mathrm{prec}(\infty) = \frac{\pi}{\pi + (1-\pi)\,\frac{c_0}{c_0 + e^{-\phi\rho t}}} .
$$

**Absence is evidence too.** A break seen in $B$ with both homologues present gives

$$
\Lambda_-(t) = \frac{1-e^{-\phi\rho t}}{1-e^{-\rho t}} \xrightarrow{t\to0} \phi ,
$$

that is, $\ln\phi \approx -0.63$ for $\phi = 0.53$. A break in a *close* relative is fairly strong
evidence *against* an operon. Breaks are rare at small $t$, so the expected information is small.
Union ("conserved in any genome") rules discard this negative evidence completely. A union rule
therefore gets *worse* as close genomes are added: every extra genome can only add positives, at
base-rate precision. That matches the drop from 0.850 to 0.69-0.74 in `phylo_scaling.md`.

**Expected information per genome.** One genome's expected evidence is the KL divergence between the
Bernoulli outcomes under $H_1$ and $H_0$:

$$
I(t) = p_1\ln\frac{p_1}{p_0} + (1-p_1)\ln\frac{1-p_1}{1-p_0}, \quad p_1 = e^{-\phi x},\; p_0 = e^{-x},\; x = \rho t .
$$

For $\phi = 0.53$:

| $x$ | 0.1 | 0.5 | 1 | 2 | 3 | 4 | 6 |
|---|---|---|---|---|---|---|---|
| $I(x)$, nats | 0.013 | 0.058 | 0.10 | 0.143 | 0.147 | 0.129 | 0.079 |

$I$ is linear in $x$ near 0, peaks at $x \approx 2$-3 (5-14% of neutral adjacencies retained), and
then decays as both outcomes become rare. Homologue loss $h(t)$ adds a further factor $e^{-\eta t}$
and pulls the optimum inwards. For E. coli, $x \approx 2$-3 means Vibrio, Haemophilus, Pseudomonas
(levels 4-5). A single genome there is the most informative *per genome*. The best *panels* are still
deeper and broader, because nearby genomes are not independent of one another (1.3).

### 1.2 Self-calibrated divergence (no ANI needed)

Most adjacencies in any genome are neutral, i.e. not selected (for E. coli, 39% of labelled same-strand
neighbours are different-TU, and every opposite-strand neighbour). Their retention therefore *measures*
$\rho t$:

$$
\hat x_{AB} = -\ln \hat P(\text{cons}\mid \text{neutral}) \approx \rho t_{AB} + \eta t_{AB} - \ln(1-\varepsilon).
$$

In practice, use all of $A$'s adjacencies, or those with a long gap or opposite strands as a proxy for
neutral ones. This is a genome-rearrangement distance on the same footing as the evidence. It is better
than ANI, which saturates beyond the genus level, and better than taxonomic rank. B. subtilis shows the
rank problem: "outside the GTDB order" is still too close there (`phylo_scaling_bsub.md`). The
model's form of the distance-floor rule is to **set the floor on $\hat x$, not on rank**.

### 1.3 Many genomes: dependence, phylogenetic diversity, effective number of lineages

Adjacency states across genomes are not independent: genomes that share a branch share its breaks.
Model the adjacency as a two-state character (present/broken) on the tree $T$. It is lost at rate
$\mu \in \{\rho, \phi\rho\}$ and regained at a negligible rate (Dollo-like; $c_0$ is handled as a small
gain rate). Condition on the families being present at every tip considered (tips where they are not
are pruned, since they are uninformative about rearrangement). The likelihood comes from Felsenstein
pruning. Under Dollo and parsimony it has a closed form:

$$
\ln L(\mu) = -\mu\,\mathrm{PD}_{+} + \sum_{b \in \mathcal{B}} \ln\!\left(1-e^{-\mu \ell_b}\right) ,
$$

- $\mathrm{PD}_+$ is the total branch length of the subtree, rooted at the MRCA, over which the
  adjacency is retained (Faith's phylogenetic diversity of the conserving tips, including the path to
  the query).
- $\mathcal{B}$ is the set of branches with an inferred break, of lengths $\ell_b$.

The log Bayes factor for selection is

$$
\boxed{\;\ln \mathrm{BF}_{10} = (1-\phi)\rho\,\mathrm{PD}_+ \;+\; \sum_{b\in\mathcal B}\ln\frac{1-e^{-\phi\rho\ell_b}}{1-e^{-\rho\ell_b}}\;}
\;\approx\; (1-\phi)\rho\,\mathrm{PD}_+ \;-\; |\mathcal B|\ln\frac1\phi \quad(\ell_b \text{ short}).
$$

Consequences:

1. **The sufficient statistic is PD, not a count.** Adding $n$ genomes from one genus adds the length
   of their crown, which is close to 0. Adding one genome from another phylum adds a path of the order
   of the tree depth. The catalogue's "distinct partner genera" approximates PD by counting tips of the
   tree cut at the genus level.
2. **Effective number of lineages.** The counting-based version is the phylogenetic effective sample
   size. If $C$ is the tip covariance matrix (shared path lengths divided by depth), then
   $n_\mathrm{eff} = \mathbf 1^\top C^{-1}\mathbf 1$, between 1 (a star of identical genomes) and $n$
   (a star with long branches). Weighting genome $k$ by $w = C^{-1}\mathbf 1$ is the
   independent-contrasts / GLS weighting. For a panel of distant lineages, $\mathrm{PD}_+ \approx \sum$
   of branch lengths to the deepest nodes, and $n_\mathrm{eff}$ is about the number of distinct
   classes or phyla. This is the "spread over classes and phyla" result.
3. **Majority votes over close genomes cannot help.** They all inherit the same state from one short
   stem. `phylo_scaling` saw exactly this (level 3, $k \ge 2$: precision 0.707 → 0.711).

### 1.4 Expected recall versus panel size

Recall here is the fraction of same-TU pairs conserved in at least one of $n$ distant panel genomes
(the union call, which has high precision in the distant regime). Pair $i$ has per-genome conservation
probability $p_i = (1-\varepsilon) h_i e^{-\phi_i x}$, averaged over the panel's $x$. For exchangeable
genomes:

$$
R(n) = \mathbb E_i\!\left[1-(1-p_i)^n\right].
$$

Two things spread the $p_i$ over orders of magnitude: the operon's age and strength ($\phi_i$, and
whether it formed after the panel floor), and the phyletic breadth of the families ($\eta_i$). If
$\ln p_i$ is roughly uniform on $[\ln p_\min, \ln p_\max]$ for a recoverable fraction $Q$ of pairs (the
rest have $p=0$: lineage-specific operons younger than the floor, or genes absent outside the clade),
then

$$
R(n) \approx \frac{Q}{\ln(p_\max/p_\min)}\left[\ln n + \gamma_E + \ln p_\max\right],
\qquad \frac{1}{p_\max} \ll n \ll \frac{1}{p_\min},
$$

with $\gamma_E$ the Euler-Mascheroni constant. Recall is therefore linear in $\log n$, with slope
$\beta = Q/\ln(p_\max/p_\min)$. It saturates at $R(\infty) = Q$ once $n \gtrsim 1/p_\min$.

PD gives a second route to the same curve. For coalescent- or Yule-like trees, the expected PD of $n$
random tips grows as $\sim 2\sum_{k<n}1/k \approx 2\ln n$. An adjacency is detected once
$(1-\phi)\rho\,\mathrm{PD}_+$ exceeds a threshold, so recall is again linear in $\log n$.

Fits (`theory_cluster/check.py`, per-pair $p_i$ from the distant pools):

| | E. coli (L7-9, 30 genomes) | B. subtilis (L6-9, 59 genomes) |
|---|---|---|
| empirical slope $\beta$ per $\ln n$ | 0.094 | 0.055 |
| log-uniform fit on $n\le10$ | $Q=0.82$, $p_\min=10^{-3.6}$, $p_\max=10^{-0.1}$ | $Q=0.32$, $p_\min=10^{-1.9}$, $p_\max=1$ |
| predicted $R(20), R(30)$ | 0.342, 0.382 | 0.246, 0.267 |
| observed (random panels, `*_curve.tsv`) | 0.364, 0.402 | 0.25 (20), 0.29 (40) |
| pairs conserved in none of the pool | 60% | 69% |

The slope is well determined. $Q$ and $p_\min$ are not separately identifiable from $n\le 30$: they
trade off along a ridge. The E. coli fit puts saturation beyond $10^3$ distant lineages, at a recall
ceiling near 0.8. The B. subtilis fit says long-spacer Bacillus operons mostly have $p=0$ outside the
Bacilli, so recall should plateau around 0.3 by $n \sim 100$. This fits the observation that
neighbourhood evidence barely helps B. subtilis accuracy (+0.003). The ~0.88 accuracy ceiling in
`phylo_scaling.md` follows: the pairs left at gap 51-300 bp mostly belong to the $1-Q$ class, which no
number of distant genomes recovers.

---

## 2. A generative model of gene clusters

### 2.1 Model

- **Units.** A latent set $\mathcal U$. Unit $u$ has:
  - a slot list $(s_1,\dots,s_{L_u})$; each slot has a distribution $\theta_{u,s}$ over gene families
    (one family for a specific operon, a vocabulary of paralogous families for a template such as
    "ABC importer = SBP + 2 TMD + NBD");
  - a consensus order and strands, with a Mallows dispersion $\beta_u$ for local permutations and
    inversions;
  - a per-slot inclusion probability $\kappa_{u,s}$ (optional or lost genes);
  - a junction retention exponent $\phi_u$ for its internal adjacencies.
- **Presence.** On the species tree, each genome has instance count $z_{g,u} \in \{0,1,2,\dots\}$.
  Presence follows a gain/loss CTMC on the tree (rates $\lambda_u, \mu_u$) plus HGT gains. Its
  stationary prevalence is $\pi_u$. This part is already the presence/pangenome side of the project
  (PPanGGOLiN partitions, the phylogenetic occupancy model).
- **Instantiation.** Each instance draws families per slot ($\theta$, $\kappa$) and an order (Mallows
  around the consensus). It is inserted as a contiguous block in the chromosome, or split across
  contigs by assembly (an observed truncation, not a biological event).
- **Background.** Genes outside units, and the order of units relative to each other, follow the
  neutral rearrangement process of section 1 ($\phi = 1$).

The likelihood of the observed genomes $G = \{g\}$ (gene-family sequences on contigs, with contig ends)
given tree $T$ is

$$
\mathcal L(\mathcal U \mid G, T) = \sum_{\text{assignments } a}\;\prod_{u}\Big[ P_\text{tree}(z_{\cdot,u}\mid \lambda_u,\mu_u,T)\;
\prod_{\text{instances } k\text{ of } u} P(\text{content}_k\mid\theta_u,\kappa_u)\,P(\text{order}_k\mid\beta_u)\Big]\;
P_\text{adj}(\text{junctions}\mid \phi, \rho, T)\;P_0(\text{background}),
$$

where the assignment $a$ maps each gene to a unit instance or to the background. $P_\text{adj}$ is the
product over adjacencies of the pruning likelihoods in 1.3. This is intractable as written: the sum is
over segmentations of every genome and the adjacency characters are coupled across genomes. The
estimator in section 3 is a **composite likelihood**:

1. Score adjacencies one at a time, from their PD Bayes factors (1.3).
2. Segment each genome given those scores. This gives *loci*: instances with uncertain unit labels.
3. Cluster loci into units with a pairwise (Potts / correlation clustering) approximation to the
   instance terms.

### 2.2 Loci from adjacency evidence (replaces support ≥ 1 and the dip split)

Along each contig, run a two-state HMM. Each adjacency $a$ emits evidence
$e_a = \ln \mathrm{BF}_{10}(a)$ from 1.3, plus the gap term (intergenic distance, strand). The states
are *inside a unit* (junction selected) and *boundary* (neutral junction). The posterior
$P(\text{inside}\mid e)$ gives locus boundaries, and Viterbi gives the segmentation.

- The current rule "support ≥ 1 cross-genus partner" is this HMM with $e_a$ replaced by an indicator
  of PD above a genus-depth threshold.
- The "split where support < 0.34 × the smaller side maximum" rule approximates the posterior dropping
  at a junction whose evidence is low *relative to* its flanks. Nested modules (abutting operons that
  form a superoperon in some lineages) are exactly the adjacencies with intermediate $\phi$. The HMM
  can carry three states (intra-operon / inter-operon-conserved / neutral), and its third state is the
  superoperon junction.

### 2.3 What Leiden-on-content-similarity approximates

Take loci $\ell, m$ with $n_\ell, n_m$ genes, of which $k$ are matched (same family or a Spacedust hit).
Under "same unit", each gene has a counterpart with probability $q_s$ (slot inclusion × homology
detection). Under "different units", a gene is matched only by chance, with probability $q_0$, which
depends on how many paralogous families of that kind exist. With $n = n_\ell + n_m$ and
$k = m_\ell + m_m$:

$$
\mathrm{LLR}_{\ell m} = k\ln\frac{q_s}{q_0} + (n-k)\ln\frac{1-q_s}{1-q_0}
= n\,\ln\!\frac{q_s(1-q_0)}{q_0(1-q_s)}\,\Big(\underbrace{\tfrac{k}{n}}_{\text{catalogue similarity}} - \gamma\Big),
\qquad
\gamma = \frac{\ln\frac{1-q_0}{1-q_s}}{\ln\frac{q_s(1-q_0)}{q_0(1-q_s)}} .
$$

The MAP partition of loci into units, under a pairwise model with independent same/different evidence
and a uniform prior, maximises $\sum_{\ell<m,\,c_\ell = c_m}\mathrm{LLR}_{\ell m}$. This is correlation
clustering. CPM-Leiden maximises $\sum_{\text{same}}(w_{\ell m}-\gamma)$ with $w = k/n$. That is the
same objective up to a per-pair factor $n\ln(\cdot)$, with **resolution $\gamma$ equal to the
similarity at which same/different are equiprobable**. The catalogue's CPM on (mA+mB)/(nA+nB) is
therefore correlation clustering with a linearised binomial LLR and one global $\gamma$. Three
simplifications follow from that:

1. **One $\gamma$ cannot fit all families.** Values of $\gamma(q_s,q_0)$:

   | $q_s$ | $q_0$ | $\gamma$ |
   |---|---|---|
   | 0.7 | 0.02 | 0.25 |
   | 0.7 | 0.10 | 0.36 |
   | 0.7 | 0.30 | 0.50 |

   Hub families (ABC ATPases, SDR, response regulators) have high $q_0$, so their loci need
   $\gamma \approx 0.4$-0.5 to split, while universal operons with lineage-split families need
   $\gamma \lesssim 0.25$ to stay together. Two Leiden levels (0.1 / 0.4) are an empirical bracket
   around this. The principled fix is a per-pair weight
   $\mathrm{LLR}_{\ell m}$ with $q_0$ estimated per family from its copy number in the genomes
   (expected chance-match rate). Then run CPM at $\gamma = 0$ on LLR weights, plus a prior term
   ($\ln$ of the prior odds of a new unit).
2. **The ≥ 2-matched-genes rule and the 90%-best-match filter** are crude ways to stop high-$q_0$
   single matches from adding positive weight. The LLR does this by construction.
3. **Pairs are not independent.** Loci from one genus are near-copies by descent, and their pairwise
   LLRs are redundant, exactly as in 1.3. Weight each locus by its phylogenetic weight $w_\ell$
   ($C^{-1}\mathbf 1$ over the bins carrying that unit), or collapse within-genus loci to one
   representative before clustering. The "same-bin edges dropped" and cross-genus rules are partial
   versions of this.

### 2.4 Split/merge criteria (Bayes factor ≈ MDL)

**Families/units.** With a Beta-Bernoulli prior on each slot's family inclusion, the marginal
likelihood of a unit's loci $S$ (phylogenetically weighted counts $\tilde n$, $\tilde k_f$) is

$$
M(S) = \prod_{f}\frac{B(a+\tilde k_f,\; b+\tilde n_S-\tilde k_f)}{B(a,b)}\;\times M_\text{order}(S),
$$

with a Mallows marginal for order (or omit order at first). Split $S$ into $S_1, S_2$ if

$$
\ln\mathrm{BF}_\text{split} = \ln M(S_1)+\ln M(S_2)-\ln M(S) + \ln\frac{P_\text{prior}(2\text{ units})}{P_\text{prior}(1\text{ unit})} > 0 .
$$

Use a Chinese restaurant or Pitman-Yor prior; the Pitman-Yor discount suits the heavy-tailed family
sizes (4,121 singletons). In MDL terms, $-\ln M$ is the code length of the loci given the unit, and
the prior is the cost of describing an extra unit, about $\tfrac12 |F_u| \ln \tilde n$ under BIC. The
BF prefers a split when the two halves use **mutually exclusive** families at the same slots:

- the Na⁺ and H⁺ ATP synthases (different β-subunit families at an otherwise shared slot structure)
  split;
- lineage-specific divergence of a family that is shared at a coarser level (30% families splitting
  uL2 by phylum) does not split, provided a coarser homology (Spacedust hit or profile) is used for
  matching.

**Nested and overlapping units.**

- *Superoperons* (S10-spc-alpha; alpha operon on its own as CF00030). Model a superoperon as a unit of
  units: blocks joined by junctions with their own $\phi_J$. Split into blocks if the junction's
  $\ln \mathrm{BF}_{10}$ (1.3) is below the prior log odds while the blocks' internal adjacencies are
  above it. Report both levels: the decision is about the junction, not about the content. The alpha
  operon is unit-level, and the S10-spc-alpha junction is conserved in some lineages and broken in
  others, which shows as many inferred breaks $|\mathcal B|$ on a long PD.
- *Combinatorial transporters* (ABC, TRAP). Model them as a template unit whose slots draw from
  vocabularies (SBP families × TMD families × NBD families), and instances as specific triples. Two
  levels:
  - **template**: share slot *superfamilies*, low $\gamma$;
  - **instance**: share actual families, higher $\gamma$ from high $q_0$.

  This is the 0.1/0.4 Leiden hierarchy, with a likelihood behind it. Use "split by slot" to decide
  whether the NBD is shared across systems (a shared-ATPase architecture; an overlapping unit) or
  system-specific.
- *Overlap.* A gene may belong to two units: shared NBD, or a gene sitting between two operons. A
  partition cannot express that. A feature-allocation prior (Indian buffet process) over units, with
  each locus a union of unit instances, can. A cheap approximation: after partitioning, allow a locus
  to join a second unit if the BF of adding it exceeds the threshold *for a sub-segment only*.

---

## 3. Practical algorithm (10⁵ genomes)

1. **Gene families.** MMseqs2 linclust cascade at 30-50% identity for counting. For matching, add a
   coarser homology level: Spacedust hits within the reference panel, or profile/Foldseek clusters for
   hub families. Per-family copy numbers give $q_0$.
2. **Lineage skeleton.** Use the GTDB tree (or one placed with GTDB-Tk), with branch lengths rescaled to
   rearrangement units. Calibrate $\rho$ (more exactly the composite rate $\rho+\eta$) once, on about
   $10^3$ sampled genome pairs, from $\hat x$ (1.2) against tree distance. Fit a monotone curve, per
   phylum if needed. Cost: Spacedust on $10^3$ pairs, already of the order of the existing runs.
3. **Representatives.** Close genomes add nothing, so cut the tree at $\hat x \approx 1$ (precision
   gain negligible below that). Keep one representative per lineage: about $10^3$-$10^4$ lineages for
   GTDB. Run Spacedust only among representatives. Within-lineage genomes feed the presence model.
   This cut is the principled form of "outside the order / outside the class".
4. **Adjacency table.** Each representative contributes about 4,000 family-pair keys
   $(f_1,f_2,\text{relative strand})$, about $4\times10^7$ rows for $10^4$ representatives. Bridging
   one gene (the catalogue's 1-gene bridge) adds keys for $(f_1, f_3)$. Sort by key (external sort). For
   each key, the set of representatives with the adjacency, and the set with both families present,
   define the character on the induced subtree. Spacedust clusters replace exact family keys for
   remote homology: each cluster contributes its member pairs' adjacencies.
5. **Score.** For each key (or each query adjacency), run Dollo pruning on the induced subtree:
   $O(k\log k)$ for $k$ tips with both families. This gives $\mathrm{PD}_+$, $|\mathcal B|$ and
   $\ln \mathrm{BF}_{10}$, and with a global or per-family-class $\phi$ prior the posterior
   $P(\text{selected})$. Fit $\phi$ by EM over all keys: a two-component mixture, neutral $\phi=1$
   plus selected $\phi\sim$ Beta.
6. **Loci.** HMM segmentation (2.2) on every genome, representatives and the others, using key scores.
   Genomes outside the representative set inherit scores by family-pair key: no new Spacedust.
7. **Units.** Candidate locus pairs come from shared keys (sparse; 210k for Chesterfield). Weight
   them by the binomial LLR (2.3) with per-family $q_0$ and phylogenetic locus weights. Run CPM-Leiden
   at $\gamma=0$ on LLR + prior, then refine with BF split/merge moves (2.4) on the largest and
   lowest-cohesion units.
8. **Outputs.** Units with $\hat\phi_u$, $\mathrm{PD}_u$, $n_\mathrm{eff}$, slot profiles, and
   prevalence $\pi_u$ from the presence model. For a MAG: the per-unit absence surprise, which now
   covers persistent operons too (PPanGGOLiN modules exclude them).

**Statistics to compute from existing Spacedust/MMseqs2 output now.** All of these are already
available for Chesterfield and the panels:

- per genome pair: $\hat x$ from the fraction of the query's adjacencies conserved (pairwise
  `clusters.tsv`);
- per adjacency: the conserving set and the families-present set, giving $\mathrm{PD}_+$ (with GTDB
  branch lengths) and $|\mathcal B|$;
- per family: copy-number distribution, giving $q_0$;
- per locus pair: matched counts $k$, $n$, giving the LLR.

---

## 4. Predictions

### P1: log-odds of conservation is linear in self-calibrated divergence, through the origin (checked)

Model: $\ln\frac{P(c\mid\text{op})}{P(c\mid\text{non})} = (1-\phi)\,x$, $x=-\ln P(c\mid\text{non})$,
with intercept 0 and one slope across the whole distance range. Alternative outcomes would falsify it:
a plateau (precision set by rank rather than divergence), a non-zero intercept (close genomes carrying
evidence), or different slopes per organism. Result
(`/data/scratch/metapanclust_data/theory_cluster/check_out.txt`):

| query | genomes | slope $1-\phi$ | intercept | r |
|---|---|---|---|---|
| E. coli K-12 (EcoCyc TUs) | 66 | 0.469 | -0.018 | 0.977 |
| B. subtilis 168 (SubtiWiki) | 90 | 0.471 | -0.062 | 0.967 |

Predicted vs observed single-genome precision for E. coli (slope 0.469, $\pi=0.612$):

| genome | Sakai | Salmonella | Yersinia | Vibrio | Pseudomonas | Neisseria |
|---|---|---|---|---|---|---|
| $x$ | 0.21 | 0.52 | 1.14 | 1.95 | 2.68 | 4.11 |
| predicted | 0.631 | 0.664 | 0.726 | 0.795 | 0.845 | 0.914 |
| observed | 0.624 | 0.658 | 0.715 | 0.790 | 0.871 | 0.913 |

$\phi \approx 0.53$ in both organisms: TU-internal adjacencies are retained as $P_0^{0.53}$. Because $x$
includes homologue loss, $1-\phi$ here is the combined retention elasticity, an upper bound on the
rearrangement-only selection effect if TU genes are also more often retained as homologues. The most
distant genomes (x > 5) sit slightly below the line, as expected from a convergence/HGT floor $c_0$.

**Remaining test:** compute $x$ from AAI or GTDB branch length instead of the same data. If
$-\ln P(c\mid\text{non})$ is linear in tree distance, $\rho$ is a clock; if it saturates, use $\hat x$
rather than the tree.

### P2: evidence scales with PD and lineages, not genomes (testable with existing outputs)

At a fixed panel size $n$, compare panels drawn as $n$ genomes from one class with $n$ genomes from
$n$ classes, both at $\hat x > 2$. The model predicts that the union call's recall follows $R(n_\mathrm{eff})$
in both cases. For the one-class panel, recall stays near the single-genome value plus crown PD. The
sum-of-LLR estimator ($\ln\mathrm{BF}$ via $\mathrm{PD}_+$) should beat the union rule at matched
precision, most of all when close genomes are included: the negative evidence from breaks in close
relatives ($\ln\phi$ per break) should *raise* accuracy above the gap-only baseline, where the union
rule lowered it to 0.69-0.74. Concrete test on `ext_matrix.tsv`: score each pair by
$\sum_k [c_k (1-\phi)x_k + (1-c_k)\ln\frac{1-e^{-\phi x_k}}{1-e^{-x_k}}]$ with $\phi=0.53$ and no fitted
parameters, then combine with the gap prior. Predicted: all 68 genomes give ≥ 0.878 accuracy, against
0.689 for the union rule. Falsified if adding levels 0-3 lowers accuracy under this score. Caveat: this
sum treats genomes as independent, so the close-genome terms are overcounted; collapse each level to
one effective genome as a first correction.

### P3: recall saturates at $Q$, set by operon age relative to the panel floor

Log-uniform fits on $n\le10$ predict a B. subtilis plateau near 0.31-0.32 by $n\approx100$ distant
lineages, while E. coli keeps rising (≈0.50 at 100). This can be tested with a larger panel
(GTDB-representative Spacedust runs). It also predicts the size of the gain over the gap-only baseline
for a new query organism: gain ∝ $Q \times$ (fraction of TU pairs above the gap threshold).

### P4: per-family resolution

$\gamma$ estimated from $q_0$ (family copy number) should match the CPM resolution at which each
catalogue family is most stable across seeds. Hub families should be stable at high $\gamma$, universal
operons at low $\gamma$. Testable on the catalogue by re-running Leiden on a $\gamma$ grid and taking
per-family ARI stability.

---

## Implications for the catalogue

- Replace "distinct partner genera ≥ 1" by $\ln\mathrm{BF}_{10}$ from PD with GTDB branch lengths and
  $\phi\approx0.5$. Same-genus partners then contribute their (small) crown PD instead of being
  dropped, which addresses the "genus-level cassettes are missed" caveat. Breaks in close relatives
  become negative evidence.
- Replace the 0.34 dip rule with HMM posterior boundaries. Superoperon junctions become a reportable
  third state.
- Replace CPM on similarity with CPM on binomial LLR (per-family $q_0$, $\gamma=0$), then BF
  split/merge. The 0.1/0.4 two-level output becomes template versus instance.
- Panel and representative choice: floor on $\hat x \approx 1$-2, maximise PD beyond it. Expect recall
  ∝ log(lineages) and an organism-specific ceiling $Q$.

## Files

- Check script: `/data/scratch/metapanclust_data/theory_cluster/check.py`. Output: `check_out.txt`
  (per-genome $P(c\mid\text{op})$, $P(c\mid\text{non})$, $x$, LLR; fits; recall curves). Input
  matrices: `phylo_scaling/results/ext_matrix.tsv`, `phylo_scaling_bsub/results/bsub_matrix.tsv`
  (values ≥ 1 = co-membership in one Spacedust cluster). Runtime < 1 CPU-min, nice 15, 1 thread.
