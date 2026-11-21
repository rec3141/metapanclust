"""Level 2: presence of a unit given the rest of the genome (context), one masked L2 logistic per unit.

Features are the reference's variable accessory families (carried by 2 .. N-2 genomes) minus the unit's own
families, plus completeness; C = 0.1 (context_model, pangbank_scale, tree_vs_context). Units carried by fewer
than 2 or more than N-2 reference genomes are not modelled (returns NaN).
"""
import numpy as np


def level2(ref, q_fam, q_comp, units, min_genomes=20, C=0.1):
    from sklearn.linear_model import LogisticRegression
    N = len(ref.genomes)
    out = np.full(len(units), np.nan)
    if N < min_genomes: return out
    k = ref.R.sum(1)
    feat = np.flatnonzero((ref.part != 'P') & (k >= 2) & (k <= N - 2))
    Xr = ref.R[feat].T.astype(np.float32)
    xq = q_fam[feat].astype(np.float32)
    cm, cs = ref.comp.mean(), ref.comp.std() + 1e-9
    cr = ((ref.comp - cm) / cs)[:, None]; cq = np.array([(q_comp - cm) / cs])
    fpos = {f: j for j, f in enumerate(feat)}
    ku = ref.U.sum(1)
    for i, u in enumerate(units):
        if not (2 <= ku[u] <= N - 2): continue
        keep = np.ones(len(feat), bool)
        keep[[fpos[f] for f in ref.unit_members[u] if f in fpos]] = False
        y = ref.U[u].astype(int)
        m = LogisticRegression(C=C, max_iter=300).fit(np.hstack([Xr[:, keep], cr]), y)
        out[i] = m.predict_proba(np.r_[xq[keep], cq][None])[0, 1]
    return out
