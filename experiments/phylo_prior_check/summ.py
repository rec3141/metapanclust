import sys, pandas as pd
S = pd.concat([pd.read_csv(f, sep='\t') for f in sys.argv[1:]])
P = ['sp_eb', 'shrunk', 'mix', 'sp_2state', 'anc']
b = S.groupby(['stratum', 'm', 'genus', 'focal'])[[f'bits_{p}' for p in P]].mean().reset_index()
out = b.groupby(['stratum', 'm'])[[f'bits_{p}' for p in P]].mean().round(0)
out['n_focal'] = b.groupby(['stratum', 'm']).size()
for p in ('shrunk', 'anc', 'sp_2state'):
    out[f'{p}<sp_eb'] = b.assign(w=b[f'bits_{p}'] < b.bits_sp_eb).groupby(['stratum', 'm']).w.sum().astype(int)
out['anc<shrunk'] = b.assign(w=b.bits_anc < b.bits_shrunk).groupby(['stratum', 'm']).w.sum().astype(int)
print(out.to_string())
