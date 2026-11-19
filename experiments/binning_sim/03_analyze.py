"""Summarise simulated-bin scores: per-condition means, separability (AUC) and the merge rule.
usage: 03_analyze.py OUTDIR sim_A_B.tsv [sim_B_A.tsv ...]   (ppanggolin env: pandas + matplotlib)"""
import sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out = sys.argv[1]
df = pd.concat([pd.read_csv(f, sep='\t') for f in sys.argv[2:]])
df['S_perMb'] = df['S_tot'] / df['mb']
df['Sdup_perMb'] = df['S_totdup'] / df['mb']
df['Smod_perMb'] = df['S_mod'] / df['mb']
df['Sabs1_perMb'] = df['S_abs1'] / df['mb']
# merge rule: the added piece costs more bits per Mb than the bin it joins already does
for key in ('S_tot', 'S_totdup', 'S_mod'):
    df[f'excess_{key}'] = df[f'd_{key}_perMbAdded'] - df[f'base_{key}_perMb']
df['reject_merge'] = df['excess_S_tot'] > 0
df['perMb_up'] = df['S_perMb'] > df['base_S_tot_perMb']          # merged bin more surprising per Mb
df['content_perGene'] = (df['S_pres'] + df['S_for']) / df['genes']  # absence term left out
# ANI between the base genome and the added same-species genome (skani, written by 01_build.sh)
import glob, os
ani = {}
for fn in glob.glob(os.path.join(os.path.dirname(sys.argv[2]), 'ani_*.tsv')):
    for line in open(fn):
        a, b, v = line.split()[:3]
        ani[(a.replace('.fna', ''), b.replace('.fna', ''))] = ani[(b.replace('.fna', ''), a.replace('.fna', ''))] = float(v)
df['ani'] = [ani.get((g, o), np.nan) for g, o in zip(df['genome'], df['other'].fillna(''))]


def md(t, index=True):
    """markdown table without the tabulate dependency"""
    t = t.reset_index() if index else t
    head = '| ' + ' | '.join(map(str, t.columns)) + ' |'
    sep = '|' + '---|' * len(t.columns)
    body = ['| ' + ' | '.join(str(v) for v in r) + ' |' for r in t.itertuples(index=False)]
    return '\n'.join([head, sep] + body)


def auc(pos, neg):
    """P(score of pos > score of neg), ties count half (Mann-Whitney)"""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    gt = (pos[:, None] > neg[None, :]).mean(); eq = (pos[:, None] == neg[None, :]).mean()
    return gt + eq / 2


lines = []
cols = ['comp', 'dup', 'dup_rel', 'foreign', 'content_perGene', 'S_perMb', 'perMb_up', 'Sabs1_perMb', 'Smod_perMb', 'd_S_tot_perMbAdded',
        'excess_S_tot', 'reject_merge']
for (a, b), d in df.groupby(['A', 'B']):
    lines.append(f'\n## Pangenome {a}, foreign species {b}\n')
    pc = d[d.kind.str.startswith('piece')].copy()
    d = d[~d.kind.str.startswith('piece')]
    t = d.groupby(['kind', 'frac'])[cols].mean().round(3)
    t.insert(0, 'n', d.groupby(['kind', 'frac']).size())
    lines.append(md(t))
    lines.append('\nAUC, contaminated by B vs same-species merge at equal added length '
                 '(1 = score always higher for B):\n')
    rows = []
    for base_kind, strain_kind, label in (('contam_B', 'strain', 'complete base'),
                                          ('contam_B_inc', 'strain_inc', 'base at 70%')):
        for f in sorted(d.frac[d.kind == base_kind].unique()):
            pb = d[(d.kind == base_kind) & (d.frac == f)]; ps = d[(d.kind == strain_kind) & (d.frac == f)]
            rows.append(dict(base=label, frac=f, **{c: round(auc(pb[c], ps[c]), 3) for c in
                        ('dup', 'dup_rel', 'foreign', 'content_perGene', 'S_perMb', 'Sdup_perMb', 'excess_S_tot', 'excess_S_totdup')}))
    lines.append(md(pd.DataFrame(rows), index=False))
    lines.append('\nAUC, contaminated by B vs uncontaminated bins (clean + incomplete):\n')
    neg = d[d.kind.isin(['clean', 'incomplete'])]
    rows = []
    for f in sorted(d.frac[d.kind == 'contam_B'].unique()):
        pb = d[(d.kind == 'contam_B') & (d.frac == f)]
        rows.append(dict(frac=f, **{c: round(auc(pb[c], neg[c]), 3) for c in
                    ('dup', 'dup_rel', 'foreign', 'content_perGene', 'S_perMb', 'Smod_perMb')}))
    lines.append(md(pd.DataFrame(rows), index=False))
    st = d[d.kind.isin(['strain', 'strain_inc'])]
    if st['ani'].notna().any():
        st = st.assign(ani_bin=pd.cut(st['ani'], [95, 97, 98.5, 100]))
        t = st.groupby(['ani_bin'], observed=True)[['dup', 'excess_S_tot', 'reject_merge']].mean().round(3)
        t.insert(0, 'n', st.groupby('ani_bin', observed=True).size())
        lines.append('\nSame-species merges by ANI between the two genomes:\n')
        lines.append(md(t))
    if len(pc):
        pc['len_kb'] = pd.cut(pc['frac'], [0, 10, 30, 50.1])
        pc['piece_foreign'] = np.nan
        t = pc.groupby(['kind', 'len_kb'], observed=True)[['excess_S_tot', 'reject_merge']].agg(['median', 'mean']).round(3)
        t.columns = ['_'.join(c) for c in t.columns]
        t.insert(0, 'n', pc.groupby(['kind', 'len_kb'], observed=True).size())
        lines.append('\nSingle-contig moves into the bin at 70% (frac = contig kb); reject = added bits/Mb exceed '
                     'the bin\'s own bits/Mb:\n')
        lines.append(md(t))
        rows = []
        for lb, g in pc.groupby('len_kb', observed=True):
            pos = g[g.kind == 'piece_B']; neg = g[g.kind != 'piece_B']
            rows.append(dict(len_kb=str(lb), n_B=len(pos), n_A=len(neg),
                             auc_excess=round(auc(pos['excess_S_tot'], neg['excess_S_tot']), 3),
                             auc_dup=round(auc(pos['dup'] - pos['base_dup'], neg['dup'] - neg['base_dup']), 3)))
        lines.append('\nAUC, B contig vs same-species contig (own or other strain):\n')
        lines.append(md(pd.DataFrame(rows), index=False))
    # how well does comp / undiscounted absence track true completeness?
    inc = d[d.kind.isin(['clean', 'incomplete'])].copy(); inc['true'] = 1 - inc['frac']
    r = inc[['true', 'comp', 'Sabs1_perMb', 'S_perMb', 'foreign']].corr(method='spearman')['true'].round(3)
    lines.append(f'\nSpearman rho with true completeness (clean + incomplete bins): {r.drop("true").to_dict()}\n')

    fig, ax = plt.subplots(1, 5, figsize=(21, 3.8))
    for kind, col in (('contam_B', '#c0392b'), ('strain', '#2471a3'), ('contam_B_inc', '#e59866'),
                      ('strain_inc', '#7fb3d5')):
        s = d[d.kind == kind].groupby('frac')
        for i, c in enumerate(('dup', 'dup_rel', 'foreign', 'S_perMb', 'excess_S_tot')):
            m, sd = s[c].mean(), s[c].std()
            ax[i].errorbar(m.index * 100, m.values, yerr=sd.values, label=kind, color=col, marker='o', capsize=2)
    clean = d[d.kind == 'clean']
    for i, c in enumerate(('dup', 'dup_rel', 'foreign', 'S_perMb')):
        ax[i].axhline(clean[c].mean(), color='grey', ls=':', label='clean')
    ax[4].axhline(0, color='grey', ls=':')
    for i, t_ in enumerate(('marker duplication, strict (80% id)', 'marker duplication, relaxed (30% id)',
                            'foreign gene fraction', 'surprise, bits per Mb',
                            'added bits/Mb minus base bits/Mb')):
        ax[i].set_title(t_, fontsize=10); ax[i].set_xlabel('added length, % of base')
    ax[0].legend(fontsize=7)
    fig.suptitle(f'pangenome {a}, foreign species {b}', fontsize=10)
    fig.tight_layout(); fig.savefig(f'{out}/merge_{a}_{b}.png', dpi=110); plt.close(fig)

    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    inc = d[d.kind.isin(['clean', 'incomplete'])].copy(); inc['true'] = 100 * (1 - inc['frac'])
    for i, c in enumerate(('comp', 'Sabs1_perMb', 'S_perMb')):
        g = inc.groupby('true')[c]
        ax[i].errorbar(g.mean().index, g.mean().values, yerr=g.std().values, marker='o', capsize=2)
        ax[i].set_xlabel('true completeness, % of length'); ax[i].set_title(c, fontsize=10)
    fig.suptitle(f'incompleteness, pangenome {a}', fontsize=10)
    fig.tight_layout(); fig.savefig(f'{out}/incomplete_{a}.png', dpi=110); plt.close(fig)

open(f'{out}/summary.md', 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
