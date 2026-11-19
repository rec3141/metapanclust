"""Detection probability of a gene cluster in a fragmented, incomplete MAG.

If the bin's contigs are pieces of the true genome (some pieces missing), and a cluster of span s bp
sits at a position independent of where the assembly broke, then

    P_intact(s) = sum_i max(L_i - s + 1, 0) / G,        G = sum_i L_i / completeness

exactly: of the G possible start positions, those whose whole span falls inside one binned contig.
It does not matter which pieces were lost (short-first or at random), only the observed contigs and
the true genome size. P_lost and P_partial below are approximations (checked by simulate.py).

Use in the surprise model: an absent cluster of span s is evidence of absence only with weight
P_intact; a cluster seen truncated at a contig end counts as present-but-unresolved.
    absent:  P(not seen intact | present) = 1 - p * P_intact   (p = prior presence probability)
"""
import numpy as np


def _prep(contig_lengths, completeness):
    L = np.asarray(contig_lengths, dtype=float)
    c = float(np.clip(completeness, 1e-3, 1.0))
    return L, L.sum() / c


def p_intact(span_bp, contig_lengths, completeness=1.0, beta=1.0):
    """P(every gene of a cluster spanning span_bp lies on one binned contig | cluster in genome).
    beta scales the split probability: 1 = assembly breaks independent of the cluster (exact under
    that assumption, and the conservative choice for scoring absence); Chesterfield bins suggest
    conserved clusters are split ~0.25x as often as random breaks would (breaks sit in variable DNA)."""
    L, G = _prep(contig_lengths, completeness)
    s = np.atleast_1d(np.asarray(span_bp, dtype=float))
    Ls = np.sort(L)
    cs = np.concatenate([[0.0], np.cumsum(Ls[::-1])])[::-1]  # cs[k] = sum of Ls[k:]
    k = np.searchsorted(Ls, s - 1, side='right')             # contigs with L > s - 1
    n = len(Ls) - k
    f = (cs[k] - n * (s - 1)) / L.sum()          # P(span inside one contig | on assembled DNA)
    out = (L.sum() / G) * (1 - beta * (1 - f))
    return out if np.ndim(span_bp) else float(out[0])


def p_lost(span_bp, contig_lengths, completeness=1.0):
    """Approx. P(no gene of the cluster in the bin): inside one missing piece, else two missing pieces."""
    L, G = _prep(contig_lengths, completeness)
    c = L.sum() / G
    a = p_intact(span_bp, L, 1.0)  # P(span within one piece), pieces drawn like the binned contigs
    return (1 - c) * a + (1 - c) ** 2 * (1 - a)


def detection(span_bp, contig_lengths, completeness=1.0, beta=1.0):
    """dict of intact / partial / lost probabilities for cluster span(s) in bp."""
    pi = p_intact(span_bp, contig_lengths, completeness, beta)
    pl = p_lost(span_bp, contig_lengths, completeness)
    return {'intact': pi, 'partial': np.clip(1 - pi - pl, 0, 1), 'lost': pl}


def absent_logprob(span_bp, contig_lengths, completeness, prior_present, beta=1.0):
    """-log2 P(cluster not observed intact | prior presence p): the 'surprisingly absent' score with
    fragmentation and completeness folded into the detection probability."""
    q = prior_present * p_intact(span_bp, contig_lengths, completeness, beta)
    return -np.log2(1 - q)


if __name__ == '__main__':
    import sys
    from common import fasta_lengths, checkm2, BINS
    b = sys.argv[1] if len(sys.argv) > 1 else 'dastool-lorbin_000'
    L = list(fasta_lengths(f'{BINS}/{b}.fna').values())
    c = checkm2()[b][0]
    for s in (1000, 3000, 5000, 10000, 20000, 50000):
        d = detection(s, L, c)
        print(b, s, {k: round(float(v), 3) for k, v in d.items()})
