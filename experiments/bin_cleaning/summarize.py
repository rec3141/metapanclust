"""CheckM2 original vs cleaned (and removed part) per changed bin, joined with per_bin.tsv. usage: summarize.py OUTDIR"""
import sys
import pandas as pd
out = sys.argv[1]
q = pd.read_csv(f'{out}/checkm2/quality_report.tsv', sep='\t')
q[['ver', 'bin']] = q.Name.str.split('__', n=1, expand=True)
w = q.pivot(index='bin', columns='ver', values=['Completeness', 'Contamination'])
w.columns = [f'{a[:5]}_{b}' for a, b in w.columns]
p = pd.read_csv(f'{out}/per_bin.tsv', sep='\t').set_index('bin')
t = p[['species', 'Mb', 'removed_Mb']].join(w, how='inner')
t['dCompl'] = t.Compl_clean - t.Compl_orig; t['dConta'] = t.Conta_clean - t.Conta_orig
t = t.sort_values('removed_Mb', ascending=False)
t.round(2).to_csv(f'{out}/checkm2_compare.tsv', sep='\t')
pd.set_option('display.width', 250)
print(t.round(2).to_string())
d = t
print(f'\n{len(d)} changed bins: median dCompl {d.dCompl.median():.2f}, median dConta {d.dConta.median():.2f}; '
      f'sum orig contamination {d.Conta_orig.sum():.1f} -> {d.Conta_clean.sum():.1f}; '
      f'sum completeness {d.Compl_orig.sum():.1f} -> {d.Compl_clean.sum():.1f}')
for lab, s in (('compl drop > 5', d.dCompl < -5), ('compl drop > 1', d.dCompl < -1), ('conta drop > 1', d.dConta < -1),
               ('conta rise > 0.5', d.dConta > 0.5)):
    print(f'{lab}: {int(s.sum())}  {", ".join(d.index[s])}')
def mq(c, k): return (c >= 50) & (k < 10)
def hq(c, k): return (c >= 90) & (k < 5)
print('MQ+ (>=50/<10) orig -> clean:', int(mq(d.Compl_orig, d.Conta_orig).sum()), '->', int(mq(d.Compl_clean, d.Conta_clean).sum()))
print('HQ (>=90/<5) orig -> clean:', int(hq(d.Compl_orig, d.Conta_orig).sum()), '->', int(hq(d.Compl_clean, d.Conta_clean).sum()))
