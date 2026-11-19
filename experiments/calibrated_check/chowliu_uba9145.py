"""Chain-rule / Chow-Liu joint code vs independent sum, UBA9145 pangenome (56 genomes, PPanGGOLiN modules).
Leave-one-genome-out (LEAVE_BOTH_OUT=1 also drops the other Chesterfield genome, a near-duplicate, when scoring either): marginals, pairwise tables and the tree are fitted on the other 55 genomes, then the
held-out genome is scored. Beta(1/2,1/2) pseudocounts on marginals, +1/2 per cell on 2x2 tables (PAIR_PRIOR=indep: pseudocount 2a spread as the product of the marginals).
Code: s_ind = sum_k -log2 Pr(o_k); s_T = -log2 Pr(o_root) - sum_{k != root} log2 Pr(o_k | o_parent(k))."""
import numpy as np, pandas as pd
from scipy.sparse.csgraph import minimum_spanning_tree, breadth_first_order
D = '/data/scratch/metapanclust_data/pangenomes/uba9145/ppanggolin'; OUT = '/data/scratch/metapanclust_data/calibrated_check'
rt = pd.read_csv(f'{D}/gene_presence_absence.Rtab', sep='\t', index_col=0)
fm = pd.read_csv(f'{D}/functional_modules.tsv', sep='\t'); mg = pd.read_csv(f'{D}/modules_in_genomes.tsv', sep='\t')
G = list(rt.columns)
Xfam = (rt.loc[fm.family_id.unique()] > 0).T.values.astype(int)                     # genomes x families
mods = sorted(fm.module_id.unique()); Xmod = np.zeros((len(G), len(mods)), int)
for m, g, cpl in mg.itertuples(index=False): Xmod[G.index(g), mods.index(m)] = int(cpl >= 0.5)
fmod = fm.drop_duplicates('family_id').set_index('family_id').module_id.loc[fm.family_id.unique()].values

import os
LBO = os.environ.get('LEAVE_BOTH_OUT') == '1'
PAIR_INDEP = os.environ.get('PAIR_PRIOR') == 'indep'
TWINS = ['chesterfield_marine_maxbin_019', 'chesterfield_maxbin_040']

def score(X, a=0.5):
    N, K = X.shape; res = []
    for j in range(N):
        drop = [j] + ([G.index(t) for t in TWINS if t != G[j]] if G[j] in TWINS and LBO else []); Y = np.delete(X, drop, 0); n = len(Y); x = X[j]
        p1 = (Y.sum(0)+a)/(n+2*a); lp = np.where(x == 1, np.log2(p1), np.log2(1-p1))
        s_ind = -lp.sum()
        # 2x2 joint with +a per cell
        n11 = Y.T @ Y; n10 = Y.T @ (1-Y); n01 = (1-Y).T @ Y; n00 = (1-Y).T @ (1-Y)
        if PAIR_INDEP:   # pair prior = product of the smoothed marginals (total weight 2a): no data -> no association
            q0 = 1-p1; prior = np.stack([np.outer(q0, q0), np.outer(q0, p1), np.outer(p1, q0), np.outer(p1, p1)])
            T = np.stack([n00, n01, n10, n11]).astype(float)+2*a*prior
        else:
            T = np.stack([n00, n01, n10, n11]).astype(float)+a
        T /= T.sum(0)
        pa = np.stack([T[0]+T[1], T[2]+T[3]])   # Pr(o_k) rows k
        pb = np.stack([T[0]+T[2], T[1]+T[3]])   # Pr(o_l)
        MI = sum(T[i]*np.log2(T[i]/(pa[i//2]*pb[i % 2])) for i in range(4)); np.fill_diagonal(MI, 0)
        W = MI.max()+1e-6-MI; np.fill_diagonal(W, 0)
        mst = minimum_spanning_tree(W); mst = mst+mst.T
        order, pred = breadth_first_order(mst, 0, directed=False)
        s_T = -lp[0]; pmi = 0.0
        for k in order[1:]:
            l = pred[k]; i = 2*x[k]+x[l]                       # T index: (o_k, o_l)
            cond = T[i][k, l]/pb[x[l]][k, l]                    # Pr(o_k | o_l)
            s_T -= np.log2(cond)
            pmi += np.log2(T[i][k, l]/(pa[x[k]][k, l]*pb[x[l]][k, l]))
        s_tree_marg = s_ind - pmi                             # pmi form with smoothed-pair marginals
        res.append((G[j], int(x.sum()), s_ind, s_T, s_tree_marg, mst.nnz//2, MI[np.triu_indices(K, 1)][mst.toarray()[np.triu_indices(K, 1)] > 0].mean()))
    return pd.DataFrame(res, columns=['genome', 'n_present', 'bits_indep', 'bits_chowliu', 'bits_indep_minus_pmi', 'edges', 'mean_edge_MI'])

for tag, X in [('modules', Xmod), ('module_families', Xfam)]:
    keep = (X.sum(0) > 0); X = X[:, keep]
    r = score(X); r['ratio'] = r.bits_chowliu/r.bits_indep
    r.round(3).to_csv(f'{OUT}/chowliu_{tag}' + ('_lbo' if LBO else '') + ('_pindep' if PAIR_INDEP else '') + '.tsv', sep='\t', index=False)
    print('==', tag, 'K =', X.shape[1], 'genomes =', X.shape[0])
    print(r[r.genome.str.startswith('chesterfield')].round(2).to_string(index=False))
    print('all genomes: median bits indep %.1f, Chow-Liu %.1f, ratio median %.3f (IQR %.3f-%.3f); sum indep %.0f vs CL %.0f' % (
        r.bits_indep.median(), r.bits_chowliu.median(), r.ratio.median(), r.ratio.quantile(.25), r.ratio.quantile(.75),
        r.bits_indep.sum(), r.bits_chowliu.sum()))
    print('rank of chesterfield_maxbin_040 by bits: indep %d, CL %d of %d' % (
        r.bits_indep.rank(ascending=False)[r.genome == 'chesterfield_maxbin_040'].iloc[0],
        r.bits_chowliu.rank(ascending=False)[r.genome == 'chesterfield_maxbin_040'].iloc[0], len(r)))
    from scipy.stats import spearmanr
    print('rho(indep, CL) across genomes %.2f' % spearmanr(r.bits_indep, r.bits_chowliu)[0])
