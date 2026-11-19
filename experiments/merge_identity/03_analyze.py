"""AUC of identity LR / xbar / hard-cut counts vs surprise excess, per focal x contaminant pair.
AUC = P(stat of contaminant merge > stat of strain merge) at equal added fraction and base, averaged
over fractions and the two bases (as binning_sim_genus/04_analyze.py). rej5 = share of contaminant
merges above the 95th percentile of strain merges (same focal, fraction, base).
usage: 03_analyze.py OUTDIR BSG_RESULTS
"""
import glob, sys
import numpy as np, pandas as pd

O, R = sys.argv[1:3]
STATS = ['excess', 'lr', 'lrc', 'xbar', 'xclip', 'xmed', 'for95', 'for80']
pairs = pd.read_csv(f'{R}/pairs.tsv', sep='\t')
genus_of = {}
df = pd.concat([pd.read_csv(f, sep='\t') for f in sorted(glob.glob(f'{O}/sim_*.tsv'))])
foc_sp = {a: pairs.focal[pairs.focal.str.endswith(' ' + a)].iloc[0] for a in df.A.unique()}


def auc(c, s):
    c = np.asarray(c)[:, None]; s = np.asarray(s)[None, :]
    return float(((c > s) + 0.5 * (c == s)).mean())


out = []
for A, d in df.groupby('A'):
    genus = foc_sp[A].split()[0]
    for kind in sorted(set(d.kind) - {'strain'}):
        if kind.startswith('piece'):
            continue
        cont = 'Flavobacterium psychrophilum' if kind == 'far' else f'{genus} {kind.split(":")[1]}'
        pr = pairs[(pairs.focal == foc_sp[A]) & (pairs.contaminant == cont)]
        row = dict(focal=foc_sp[A], contaminant=cont, ani=pr.ani.iloc[0] if len(pr) else np.nan,
                   auc_excess_prev=pr.auc_excess.iloc[0] if len(pr) else np.nan)
        for st in STATS:
            a = []; a20 = []; rj = []
            for (f, b), g in d[d.kind.isin(['strain', kind])].groupby(['frac', 'base']):
                c, s = g[g.kind == kind][st], g[g.kind == 'strain'][st]
                v = auc(c, s); a.append(v)
                if f == 0.2:
                    a20.append(v)
                rj.append((c > np.quantile(s, 0.95)).mean())
            row[f'auc_{st}'] = np.mean(a); row[f'auc20_{st}'] = np.mean(a20); row[f'rej5_{st}'] = np.mean(rj)
        # single 50-kb contigs vs strain contigs
        pc, ps = d[d.kind == 'piece_' + kind], d[d.kind == 'piece_strain']
        for st in ('excess', 'lr', 'lrc', 'xbar', 'xclip', 'xmed', 'for95'):
            row[f'pieceauc_{st}'] = auc(pc[st], ps[st]) if len(pc) and len(ps) else np.nan
        row['xbar_contam'] = d[d.kind == kind].xbar.median(); row['xbar_strain'] = d[d.kind == 'strain'].xbar.median()
        out.append(row)
P = pd.DataFrame(out)
P.to_csv(f'{O}/pairs.tsv', sep='\t', index=False)

bins = [(0, 83, 'no skani (< ~83%), incl. far'), (83, 85, '83-85%'), (85, 90, '85-90%'), (90, 93, '90-93%'),
        (93, 96, '93-95%')]
P['bin'] = P.ani.apply(lambda a: next((l for lo, hi, l in bins[1:] if lo <= a < hi), bins[0][2]))
pd.set_option('display.width', 250)
print('## By ANI bin (mean; min in brackets for AUC)\n')
print('| ANI bin | pairs | AUC excess | AUC LR | AUC xbar | AUC for95 | AUC for80 | AUC20 excess | AUC20 LR | rej@5% excess | rej@5% LR | rej@5% xbar | piece AUC excess | piece AUC LR |')
print('|' + '---|' * 15)
for lo, hi, l in bins:
    g = P[P.bin == l]
    f = lambda c: f'{g[c].mean():.2f} ({g[c].min():.2f})'
    print(f'| {l} | {len(g)} | {f("auc_excess")} | {f("auc_lr")} | {f("auc_xbar")} | {f("auc_for95")} | {f("auc_for80")} | '
          f'{g.auc20_excess.mean():.2f} | {g.auc20_lr.mean():.2f} | {g.rej5_excess.mean():.2f} | {g.rej5_lr.mean():.2f} | '
          f'{g.rej5_xbar.mean():.2f} | {g.pieceauc_excess.mean():.2f} | {g.pieceauc_lr.mean():.2f} |')
print('\n## Pairs with skani ANI\n')
cols = ['focal', 'contaminant', 'ani', 'auc_excess_prev', 'auc_excess', 'auc_lr', 'auc_xbar', 'auc_for95', 'auc20_lr',
        'rej5_excess', 'rej5_lr', 'rej5_xbar', 'pieceauc_excess', 'pieceauc_lr', 'pieceauc_xbar', 'xbar_strain', 'xbar_contam']
print(P[P.ani.notna()].sort_values('ani', ascending=False)[cols].round(3).to_csv(sep='\t', index=False))

# along-contig clustering
AC = pd.concat([pd.read_csv(f, sep='\t') for f in sorted(glob.glob(f'{O}/autocorr_*.tsv'))])
AC['grp'] = np.where(AC.src == 'self', 'same species', 'congener')
ani_of = {(r.focal, r.contaminant): r.ani for r in P.itertuples()}
AC['ani'] = [ani_of.get((foc_sp[a], f'{foc_sp[a].split()[0]} {s.split(":")[1]}'), np.nan) if s != 'self' else 100
             for a, s in zip(AC.A, AC.src)]
AC['grp'] = np.where(AC.src == 'self', 'same species (held-out)',
                     np.where(AC.ani >= 90, 'congener >= 90% ANI', np.where(AC.ani >= 83, 'congener 83-90%', 'congener < 83%')))
c = ['acf_x_1', 'acf_x_2', 'acf_x_5', 'acf_x_10', 'acf_x_20', 'phi_x_10', 'phi_x_25', 'phi_x_50', 'phi_for_10', 'phi_for_25', 'phi_for_50']
print('\n## Along-contig clustering (median over genome x focal-pangenome combinations)\n')
print(AC.groupby('grp')[['n_shared'] + c].median().round(2).to_csv(sep='\t'))
AC.to_csv(f'{O}/autocorr_all.tsv', sep='\t', index=False)
