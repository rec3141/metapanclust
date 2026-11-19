"""Aggregate the tree-vs-context run.

usage: aggregate.py OUTDIR SPECIES_TVC.tsv
Reads OUTDIR/res/*.{summary,matcor,diag}.tsv and *.genomes.tsv.gz; writes OUTDIR/agg/*.tsv and prints markdown.

  1. % of the prevalence log loss removed by each predictor, median over species, per N and target kind.
  2. Decomposition of the context model's gain (per species, then median):
       ctx  = gain of lr alone; rel = gain of the relatedness stack alone; both = gain of lr with relatedness offset
       unique context = both - rel   (the answer: what gene-gene context adds on top of relatedness)
       shared         = ctx + rel - both
     for sequence relatedness only (relseq: ANI + AAI) and for all (relall: + gene-content Jaccard).
  3. Which relatedness measure predicts best: per species x N, the winner among ker/tree x ani/aai/jac.
  4. Disagreement: Spearman between distance matrices; |GC deviation| vs outlyingness per measure; genomes
     outlying under ANI/AAI but typical under gene content (and vice versa), with GC deviation and whether
     sequence- or content-based relatedness predicts them better.
"""
import glob, os, sys
import numpy as np, pandas as pd

out, spf = sys.argv[1], sys.argv[2]
os.makedirs(f'{out}/agg', exist_ok=True)
sp = pd.read_csv(spf, sep='\t'); sp = sp[[c for c in ('id', 'phylum', 'species', 'genomes', 'stratum') if c in sp]]
pd.set_option('display.width', 250)


def md(df, floatfmt='.1f'):
    try: return df.to_markdown(floatfmt=floatfmt)
    except ImportError: return df.to_string()


def load(pat, **kw):
    fs = sorted(glob.glob(f'{out}/res/{pat}'))
    return pd.concat([pd.read_csv(f, sep='\t', **kw) for f in fs], ignore_index=True) if fs else pd.DataFrame()


A = load('*.summary.tsv')
A = A[A.n_targets > 0]
print(f'# tree vs context: {A.id.nunique()} species, {len(A)} summary rows, '
      f'{len(glob.glob(f"{out}/res/*.done"))} done\n')
preds = [c[5:] for c in A.columns if c.startswith('bits_')]
A['Nlab'] = np.where((A.N == A.pool) & (A.pool < 400), 'all', A.N.astype(str))
P = A.groupby(['id', 'N', 'Nlab', 'kind'])[[f'bits_{m}' for m in preds]].mean().reset_index()
for m in preds[1:]:
    P[f'pct_{m}'] = 100 * (P.bits_prev - P[f'bits_{m}']) / P.bits_prev
P = P.merge(sp, on='id', how='left')
# N bins: subsample sizes; 'all' rows (pool < 400) go to the largest listed size they reach
P['Nbin'] = np.where(P.Nlab == 'all', P.N.map(lambda n: max([s for s in (25, 50, 100, 200, 400) if s <= n] or [25])),
                     P.N).astype(int)
P.to_csv(f'{out}/agg/per_species.tsv', sep='\t', index=False)

print('## 1. % of prevalence loss removed (median over species; mean over subsamples within species)\n')
for kind in ('module', 'family'):
    t = P[P.kind == kind].groupby('Nbin')[[f'pct_{m}' for m in preds[1:]]].median()
    t.columns = [c[4:] for c in t.columns]
    t.insert(0, 'species', P[P.kind == kind].groupby('Nbin').id.nunique())
    print(f'### {kind} targets\n'); print(md(t.T)); print()
    t.to_csv(f'{out}/agg/pct_{kind}.tsv', sep='\t')

print('## 2. Decomposition of the context gain (% of prevalence loss, median over species)\n')
rows = []
for rel in ('relseq', 'relall'):
    D = P.copy()
    D['ctx'] = D.pct_lr; D['rel'] = D[f'pct_{rel}']; D['both'] = D[f'pct_lr_{rel}']
    D['unique_ctx'] = D.both - D.rel; D['unique_rel'] = D.both - D.ctx; D['shared'] = D.ctx + D.rel - D.both
    D['ctx_beats_rel'] = (D.both > D.rel + 1).astype(float)
    for (kind, N), g in D.groupby(['kind', 'Nbin']):
        rows.append(dict(relatedness=rel, kind=kind, N=N, species=g.id.nunique(),
                         **{c: g[c].median() for c in ('ctx', 'rel', 'both', 'unique_ctx', 'unique_rel', 'shared')},
                         frac_species_ctx_adds_gt1pt=g.ctx_beats_rel.mean(),
                         avg_pct=g[f'pct_avg_{rel}'].median()))
Dec = pd.DataFrame(rows)
Dec.to_csv(f'{out}/agg/decomposition.tsv', sep='\t', index=False)
print(md(Dec.set_index(['relatedness', 'kind', 'N']), '.2f')); print()

print('## 3. Best single relatedness predictor (share of species x N where it has the lowest bits)\n')
single = [m for m in preds if m.startswith(('ker_', 'tree_'))]
W = P[['id', 'kind', 'Nbin'] + [f'bits_{m}' for m in single]].copy()
W['winner'] = W[[f'bits_{m}' for m in single]].idxmin(1).str[5:]
wt = pd.crosstab([W.kind, W.Nbin], W.winner, normalize='index').round(2)
wt.to_csv(f'{out}/agg/winner.tsv', sep='\t')
print(md(wt, '.2f')); print()
seqv = [m for m in single if not m.endswith('jac')]
if seqv:
    W['best_seq'] = W[[f'bits_{m}' for m in seqv]].min(1)
    W['best_jac'] = W[['bits_ker_jac', 'bits_tree_jac']].min(1)
    W['jac_minus_seq_pct'] = 100 * (W.best_seq - W.best_jac) / W.best_seq
    print('Median % by which the best gene-content relatedness beats the best sequence relatedness '
          '(positive = gene content better):\n')
    print(md(W.groupby(['kind', 'Nbin']).jac_minus_seq_pct.median().unstack(), '.1f')); print()

print('## 4. Disagreement between relatedness measures\n')
C = load('*.matcor.tsv')
if len(C):
    C = C.merge(sp, on='id', how='left')
    c1 = C[C.a != 'absGCdev'].groupby(['a', 'b']).spearman.describe()[['count', '25%', '50%', '75%']]
    print('Spearman correlation of pairwise distances (per species, pool of <= 400 genomes):\n')
    print(md(c1, '.3f')); print()
    c2 = C[C.a == 'absGCdev'].groupby('b').spearman.describe()[['count', '25%', '50%', '75%']]
    print('Spearman of |GC - species median| with each genome\'s mean-distance z-score:\n')
    print(md(c2, '.3f')); print()
    C.to_csv(f'{out}/agg/matcor.tsv', sep='\t', index=False)

Dg = load('*.diag.tsv')
G = load('*.genomes.tsv.gz')
if len(Dg) and len(G):
    # each genome's bits at the largest subsample size it was scored in, averaged over subsamples
    G['Nmax'] = G.groupby('id').N.transform('max')
    Gm = G[G.N == G.Nmax].groupby(['id', 'genome'])[[c for c in G.columns if c.startswith('bits_')]].mean().reset_index()
    X = Dg.merge(Gm, on=['id', 'genome'], how='left').merge(sp[['id', 'species', 'phylum']], on='id', how='left')
    seqm = [m for m in ('ani', 'aai') if f'bits_tree_{m}' in X]
    for m in seqm + ['jac']:
        X[f'best_{m}'] = X[[f'bits_ker_{m}', f'bits_tree_{m}']].min(1)
    X.to_csv(f'{out}/agg/genome_diag.tsv.gz', sep='\t', index=False)
    flags = [c for c in X.columns if c.startswith('flag_')]
    rows = []
    for f in flags:
        g = X[X[f] == True]; rest = X[X[f] != True]
        r = dict(flag=f, genomes=len(g), species=g.id.nunique(),
                 median_abs_gc_dev=g.gc_dev.abs().median(), rest_median_abs_gc_dev=rest.gc_dev.abs().median())
        for m in seqm + ['jac']:
            r[f'bits_best_{m}'] = g[f'best_{m}'].mean(); r[f'rest_bits_best_{m}'] = rest[f'best_{m}'].mean()
        if 'bits_lr' in g: r['bits_lr'] = g.bits_lr.mean(); r['bits_prev'] = g.bits_prev.mean()
        rows.append(r)
    F = pd.DataFrame(rows)
    print('Discordant genomes (flag = outlying under one measure, typical under the other) and which relatedness '
          'predicts them better (mean bits per target, lower is better; "rest" = all other genomes):\n')
    print(md(F.set_index('flag').T, '.3f')); print()
    F.to_csv(f'{out}/agg/discordant_summary.tsv', sep='\t', index=False)
    L = X[X[flags].any(axis=1)].copy()
    keep = ['species', 'phylum', 'genome', 'gc', 'gc_dev'] + flags + \
           [c for c in X.columns if c.startswith(('meanz_', 'nnpct_'))] + [f'best_{m}' for m in seqm + ['jac']] + \
           [c for c in ('bits_prev', 'bits_lr') if c in X]
    L[keep].sort_values('gc_dev', key=np.abs, ascending=False).to_csv(f'{out}/agg/discordant_genomes.tsv',
                                                                        sep='\t', index=False)
    print(f'{len(L)} discordant genomes listed in agg/discordant_genomes.tsv; top 15 by |GC deviation|:\n')
    print(md(L[keep].sort_values('gc_dev', key=np.abs, ascending=False).head(15).set_index('genome'), '.3f'))
    # over all genomes: does the sequence-vs-content advantage track GC deviation or ANI-vs-content outlyingness?
    from scipy.stats import spearmanr
    print()
    for m in seqm:
        d = X.dropna(subset=[f'best_{m}', 'best_jac'])
        adv = d[f'best_{m}'] - d.best_jac          # > 0: content-based relatedness predicts this genome better
        print(f'all genomes: Spearman(best_{m} - best_jac, |GC dev|) = {spearmanr(adv, d.gc_dev.abs()).statistic:.3f}; '
              f'with meanz_{m} - meanz_jac = {spearmanr(adv, d[f"meanz_{m}"] - d.meanz_jac).statistic:.3f} '
              f'(n = {len(d)})')
