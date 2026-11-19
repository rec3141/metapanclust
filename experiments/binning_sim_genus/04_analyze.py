"""Merge-rule performance per (focal species, contaminant species) against contaminant ANI.
usage: 04_analyze.py WORKDIR RESULTS_DIR   (reads RESULTS_DIR/sim_*.tsv, WORKDIR/ani.tsv, genomes.tsv)
writes RESULTS_DIR/pairs.tsv, pairs_by_frac.tsv, pieces.tsv and prints markdown tables
"""
import glob, sys
import numpy as np, pandas as pd

W, R = sys.argv[1:3]


def md(t, index=True):
    """markdown table without the optional tabulate dependency"""
    t = t.reset_index() if index else t
    cols = [' / '.join(map(str, c)) if isinstance(c, tuple) else str(c) for c in t.columns]
    out = ['| ' + ' | '.join(cols) + ' |', '|' + '---|' * len(cols)]
    out += ['| ' + ' | '.join('' if pd.isna(v) else str(v) for v in r) + ' |' for r in t.itertuples(index=False)]
    return '\n'.join(out)


def auc(pos, neg):
    """P(pos > neg) + 0.5 P(tie)"""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    gt = (pos[:, None] > neg[None, :]).mean(); eq = (pos[:, None] == neg[None, :]).mean()
    return gt + 0.5 * eq


G = pd.read_csv(f'{W}/genomes.tsv', sep='\t', header=None, names=['acc', 'species', 'genus', 'role', 'tag'], dtype=str).fillna('')
G['ncbi'] = G.acc.str[3:]
sp_of = dict(zip(G.ncbi, G.species))
tag_of = dict(zip(G.ncbi, G.tag))
ani = pd.read_csv(f'{W}/ani.tsv', sep='\t', header=None, names=['ref', 'qry', 'ani', 'af_ref', 'af_q'])
ani['focal'] = ani.ref.map(tag_of); ani['qsp'] = ani.qry.map(sp_of)
ani_pair = ani.groupby(['focal', 'qsp']).agg(ani=('ani', 'mean'), af=('af_q', 'mean')).reset_index()

df = pd.concat([pd.read_csv(f, sep='\t') for f in sorted(glob.glob(f'{R}/sim_*.tsv')) if 'strictfrac' not in f])
sf = pd.concat([pd.read_csv(f, sep='\t').assign(A=f.split('sim_')[1].split('_strictfrac')[0])
                for f in sorted(glob.glob(f'{R}/sim_*_strictfrac.tsv'))])
foc_sp = G[G.tag != ''].groupby('tag').species.first()
genus_of = G[G.tag != ''].groupby('tag').genus.first()

merge = df[~df.kind.isin(['clean']) & ~df.kind.str.startswith('piece_')]
rows, rows_f = [], []
for (A, kind), d in merge.groupby(['A', 'kind']):
    if kind == 'strain':
        continue
    qsp = foc_sp[A].split(' ')[0] + ' ' + kind.split(':', 1)[1] if kind.startswith('genus:') else 'Flavobacterium psychrophilum'
    a = ani_pair[(ani_pair.focal == A) & (ani_pair.qsp == qsp)]
    s = merge[(merge.A == A) & (merge.kind == 'strain')]
    row = dict(genus=genus_of[A], focal=foc_sp[A], contaminant=qsp,
               ani=a.ani.iloc[0] if len(a) else np.nan, af=a.af.iloc[0] if len(a) else np.nan,
               strict_frac=sf[(sf.A == A) & (sf.src == kind)].strict_frac.mean(),
               strain_strict_frac=sf[(sf.A == A) & (sf.src == 'self')].strict_frac.mean())
    for metric in ('excess', 'excess_content', 'd_foreign', 'd_dup_rel', 'd_dup'):
        row[f'auc_{metric}'] = np.nanmean([auc(d[(d.frac == f) & (d.base == b)][metric], s[(s.frac == f) & (s.base == b)][metric])
                                           for f in d.frac.unique() for b in d.base.unique()])
    row['rej_contam'] = (d.excess > 0).mean(); row['rej_strain'] = (s.excess > 0).mean()
    row['rejC_contam'] = (d.excess_content > 0).mean(); row['rejC_strain'] = (s.excess_content > 0).mean()
    row['contam_bits_perMb'] = d.d_S_perMbAdded.median(); row['strain_bits_perMb'] = s.d_S_perMbAdded.median()
    row['base_bits_perMb'] = d.base_perMb.median()
    rows.append(row)
    for f in sorted(d.frac.unique()):
        df_, sf_ = d[d.frac == f], s[s.frac == f]
        rows_f.append(dict(focal=foc_sp[A], contaminant=qsp, ani=row['ani'], frac=f,
                           auc_excess=np.nanmean([auc(df_[df_.base == b].excess, sf_[sf_.base == b].excess) for b in ('full', '70')]),
                           rej_contam=(df_.excess > 0).mean(), rej_strain=(sf_.excess > 0).mean()))
P = pd.DataFrame(rows).sort_values(['genus', 'focal', 'ani'], ascending=[True, True, False])
P.to_csv(f'{R}/pairs.tsv', sep='\t', index=False)
PF = pd.DataFrame(rows_f); PF.to_csv(f'{R}/pairs_by_frac.tsv', sep='\t', index=False)

# single-contig moves (>= 20 kb pieces): contaminant pieces vs same-species pieces (own + strain)
pc = df[df.kind.str.startswith('piece_')]
prow = []
for A, d in pc.groupby('A'):
    same = d[d.kind.isin(['piece_own', 'piece_strain'])]
    for kind, dk in d.groupby('kind'):
        if kind in ('piece_own', 'piece_strain'):
            continue
        src = kind[6:]
        qsp = foc_sp[A].split(' ')[0] + ' ' + src.split(':', 1)[1] if src.startswith('genus:') else 'Flavobacterium psychrophilum'
        a = ani_pair[(ani_pair.focal == A) & (ani_pair.qsp == qsp)]
        prow.append(dict(focal=foc_sp[A], contaminant=qsp, ani=a.ani.iloc[0] if len(a) else np.nan,
                         auc_excess=auc(dk.excess, same.excess), rej_contam=(dk.excess > 0).mean(),
                         rej_own=(d[d.kind == 'piece_own'].excess > 0).mean(), rej_strain=(d[d.kind == 'piece_strain'].excess > 0).mean()))
PC = pd.DataFrame(prow); PC.to_csv(f'{R}/pieces.tsv', sep='\t', index=False)

pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30); pd.set_option('display.max_rows', 200)
print('## Pairs (all fractions and both bases pooled per AUC; AUC = contaminant vs strain merge at equal added length)\n')
print(md(P.round(3), False))
bins = [0, 80, 85, 90, 93, 95, 101]
lab = ['<80 / no skani', '80-85', '85-90', '90-93', '93-95', '>=95']
for T, name in ((P, 'merges'), (PC, 'single 20-50 kb contigs')):
    T = T.copy(); T['ani_bin'] = pd.cut(T.ani.fillna(0), bins, labels=lab, right=False)
    cols = [c for c in ('auc_excess', 'auc_excess_content', 'auc_d_foreign', 'auc_d_dup_rel', 'auc_d_dup', 'rej_contam', 'rej_strain', 'rejC_contam', 'rejC_strain', 'rej_own', 'strict_frac') if c in T]
    print(f'\n## By ANI bin, {name}\n')
    print(md(T.groupby('ani_bin', observed=True)[cols].agg(['mean', 'min']).round(3).assign(n=T.groupby('ani_bin', observed=True).size())))
PF['ani_bin'] = pd.cut(PF.ani.fillna(0), bins, labels=lab, right=False)
print('\n## By ANI bin and added fraction (merges)\n')
print(md(PF.pivot_table(index='ani_bin', columns='frac', values=['auc_excess', 'rej_contam', 'rej_strain'], observed=True).round(3)))
print('\n## Pieces\n'); print(md(PC.round(3), False))
