"""Surprise of a genome or bin against a species reference.

Input: a gene table for the query (contig, start, stop, fam = reference family index or -1 if unmapped,
pident = identity of the best reference hit) and the query's contig lengths.

For every reference unit u (module block or lone family):
  p_u      P(unit present in this genome) = mean of level 1 (gene-content tree) and level 2 (context) where
           level 2 is defined, else level 1 alone (the 50/50 average was best in tree_vs_context and xsp)
  d_u      P(detected intact | present) = p_intact(span_u, contig lengths, completeness) (fragmentation model)
  q_u      P(observed) = p_u d_u + (1 - p_u) EPS
  bits_u   -log2 q_u if observed, -log2(1 - q_u) if not
  excess   e_u = bits_u - H(q_u) (calibrated surprise); cells with q_u < 1/2 are the presence direction, the rest
           the absence direction. Z+ / Z- = sum e / sqrt(sum Var) over each direction.
Unmapped genes are novel to the species. Their expected number is the reference's own leave-one-out rate
(families carried by exactly one reference genome, per genome, scaled by completeness); the novel count is
scored as a negative-binomial upper tail.
Completeness, if not given, is the detected fraction of the reference's persistent families (weighted by
how often each is carried).
"""
import numpy as np, pandas as pd
from scipy.stats import poisson
from .detect import p_intact
from .tree import level1, platt
from .context import level2

EPS = 0.002
GENE_SPAN = 1000     # completeness is measured on single genes, so detection is relative to this span


def observe(ref, genes):
    """families and units observed in the query"""
    q_fam = np.zeros(len(ref.fams), bool)
    f = genes.fam.values
    q_fam[f[f >= 0]] = True
    q_unit = np.array([q_fam[m].mean() >= 0.5 if len(m) > 1 else q_fam[m[0]] for m in ref.unit_members])
    return q_fam, q_unit


def estimate_completeness(ref, q_fam):
    pers = np.flatnonzero(ref.part == 'P')
    w = ref.R[pers].mean(1)
    return float(np.clip((w * q_fam[pers]).sum() / max(w.sum(), 1e-9), 0.05, 1.0))


def novelty_counts(ref):
    """per-genome count of families carried only by that genome (leave-one-out novelty), genomes with
    completeness >= 0.9 (all if fewer than 5)"""
    single = ref.R.sum(1) == 1
    n = ref.R[single].sum(0)
    ok = ref.comp >= 0.9
    return n[ok] if ok.sum() >= 5 else n


def novelty_tail(n_obs, counts, c):
    """-log2 P(novel >= n_obs) under a negative binomial fitted by moments to the reference's leave-one-out
    novelty counts scaled to completeness c (singletons concentrate in a few genomes, so Poisson is wrong)"""
    from scipy.stats import nbinom
    m = max(counts.mean() * c, 0.5); v = max(counts.var() * c * c, m * 1.01)
    r = m * m / (v - m); pr = r / (r + m)
    return m, float(-np.log2(max(nbinom.sf(n_obs - 1, r, pr), 1e-300))) if n_obs > 0 else 0.0


def score(ref, genes, contig_lengths, comp=None, folds=5, seed=0, use_context=True, recalibrate=True):
    q_fam, q_unit = observe(ref, genes)
    c_est = estimate_completeness(ref, q_fam)
    c = float(comp) if comp is not None else c_est
    N = len(ref.genomes)
    nU = len(ref.unit_members)
    ku = ref.U.sum(1)
    units = np.flatnonzero(ku > 0)                     # units with no reference carrier left are novel
    calib = []
    p1_raw = level1(ref, q_fam, c, units, folds=folds, seed=seed, calib=calib)
    # recalibrate the tree posterior on the reference genomes' own leave-one-out predictions
    cal = platt(calib, ref.comp) if recalibrate else (lambda x: x)
    p1 = cal(p1_raw)
    p2 = level2(ref, q_fam, c, units) if use_context else np.full(len(units), np.nan)
    p = np.where(np.isnan(p2), p1, 0.5 * p1 + 0.5 * p2)
    L = np.asarray(list(contig_lengths), float)
    # completeness (estimated from detected genes, or CheckM2 markers) already includes genes lost at contig
    # breaks, so only the extra loss of spans longer than a gene is applied
    f0 = p_intact(GENE_SPAN, L, 1.0)
    d = np.clip(c * p_intact(np.maximum(ref.unit_span[units], GENE_SPAN), L, 1.0) / max(f0, 1e-9), 0, c)
    q = np.clip(p * d + (1 - p) * EPS, 1e-9, 1 - 1e-9)
    y = q_unit[units]
    bits = np.where(y, -np.log2(q), -np.log2(1 - q))
    H = -(q * np.log2(q) + (1 - q) * np.log2(1 - q))
    V = q * np.log2(q) ** 2 + (1 - q) * np.log2(1 - q) ** 2 - H ** 2
    e = bits - H
    pres_dir = q < 0.5
    U = pd.DataFrame({'unit': ref.unit_names[units], 'kind': ref.unit_kind[units], 'n_fam': [len(ref.unit_members[u]) for u in units],
                      'ref_carriers': ku[units], 'span': ref.unit_span[units].round(0), 'p_tree_raw': p1_raw.round(4), 'p_tree': p1.round(4),
                      'p_context': np.round(p2, 4), 'p': p.round(4), 'detect': np.round(d, 4), 'q': q.round(5),
                      'observed': y, 'bits': bits.round(3), 'excess': e.round(3), 'direction':
                      np.where(pres_dir, 'presence', 'absence')})
    # novel genes
    novel = genes.fam.values < 0
    n_nov = int(novel.sum())
    lam, nov_bits = novelty_tail(n_nov, novelty_counts(ref), c)
    # per contig: which contigs carry the presence-direction surprise
    unit_of_gene = np.where(genes.fam.values >= 0, ref.fam2unit[np.maximum(genes.fam.values, 0)], -1)
    eu = pd.Series(np.where(pres_dir & y, e, 0.0), index=units)
    g = genes.assign(unit=unit_of_gene, novel=novel,
                     persistent=(genes.fam.values >= 0) & (ref.part[np.maximum(genes.fam.values, 0)] == 'P'))
    g['unit_excess'] = g.unit.map(eu).fillna(0.0)
    # a unit's excess is shared among its genes in the query
    g['unit_excess'] /= g.groupby('unit').unit.transform('size').where(g.unit >= 0, 1)
    ctg = g.groupby('contig').agg(genes=('fam', 'size'), novel=('novel', 'sum'), persistent=('persistent', 'sum'),
                                  presence_excess=('unit_excess', 'sum'), median_pident=('pident', 'median'))
    ctg['novel_frac'] = (ctg.novel / ctg.genes).round(3)
    ctg['presence_excess'] = ctg.presence_excess.round(3)
    sp = pres_dir; sa = ~pres_dir
    summary = dict(
        genomes_in_ref=N, platt_a=round(getattr(cal, 'params', (0, 1))[0], 4), platt_b=round(getattr(cal, 'params', (0, 1))[1], 4), units_scored=len(units), completeness=round(c, 4), completeness_est=round(c_est, 4),
        contigs=len(L), Mb=round(L.sum() / 1e6, 3), genes=len(genes),
        bits_presence=round(bits[sp & y].sum(), 2), bits_absence=round(bits[sa & ~y].sum(), 2),
        excess_presence=round(e[sp].sum(), 2), Zplus=round(e[sp].sum() / np.sqrt(max(V[sp].sum(), 1e-9)), 3),
        excess_absence=round(e[sa].sum(), 2), Zminus=round(e[sa].sum() / np.sqrt(max(V[sa].sum(), 1e-9)), 3),
        rare_present=int((y & (q < 0.05)).sum()), expected_absent=int((~y & (q > 0.95)).sum()),
        novel_genes=n_nov, novel_expected=round(lam, 1), novel_bits=round(nov_bits, 2),
        median_pident=round(float(np.nanmedian(genes.pident.values)) if (~novel).any() else np.nan, 2))
    return summary, U, ctg.reset_index()
