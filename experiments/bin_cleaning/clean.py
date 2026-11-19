"""Surprise-guided bin cleaning of the 93 Chesterfield DAS Tool bins: per-contig keep/remove with a reason.

Inputs: coverage classes of every binned contig (coverage_test allbins_contigs.tsv) and per-contig gene-family
content for the 18 bins that sit in a species pangenome (content.py).

Rules, applied in order (first match decides):
  R1 off-population  -> remove, any length. >50% of depth in samples where the bin's core is absent.
                        Calibration: 1 of 763 known members (>= 2 persistent genes, bins outside TMED14/UBA9145).
  R2 no-reads        -> remove, any length. No mapped reads in any sample; 0 of 763 known members.
  R3 content+coverage-> remove. Contig >= 10 kb, <= 10% of its genes in a family any public genome of the
                        species carries, gene surprise per Mb above the bin's, and coverage does not vouch for it
                        (class not core-like, or composition outlier tnf > 1.25).
  K1 content only    -> keep (flagged). Same content test but core-like coverage and core composition.
  K2 short           -> keep. Everything < 10 kb not caught by R1/R2 (ratios, other-bin and surprise are noise there:
                        other-bin hits 5 of 97 short known members, 4.7% of all short contigs).
  K3 keep            -> everything else.
usage: clean.py OUTDIR"""
import os, sys
import numpy as np, pandas as pd
D = '/data/scratch/metapanclust_data'
out = sys.argv[1]; os.makedirs(f'{out}/cleaned', exist_ok=True)
TNF_OUT, MINLEN, PF_MAX = 1.25, 10_000, 0.10

m = pd.read_csv(f'{D}/coverage_test/allbins/allbins_contigs.tsv', sep='\t')
m['cov'] = m['class'].str.replace('+low-r', '', regex=False)
c = pd.read_csv(f'{out}/content.tsv', sep='\t')
x = m.merge(c, on=['contigName', 'bin'], how='left')
x['pf'] = x.pubfam / x.genes.clip(lower=1)
bb = x[x.species.notna()].groupby('bin').agg(bits=('bits', 'sum'), bp=('len', 'sum'))
x['bin_bits_per_Mb'] = x.bin.map(bb.bits / bb.bp * 1e6)
x['bits_per_Mb'] = x.bits / x.len * 1e6
content = x.species.notna() & (x.len >= MINLEN) & (x.pf <= PF_MAX) & (x.bits_per_Mb > x.bin_bits_per_Mb)
cov_ok = (x['cov'] == 'core-like') & (x.tnf <= TNF_OUT)

dec = np.select(
    [x['cov'] == 'off-population', x['cov'] == 'no-reads', content & ~cov_ok, content & cov_ok, x.len < MINLEN],
    ['R1_off-population', 'R2_no-reads', 'R3_content+coverage', 'K1_content-only', 'K2_short'], 'K3_keep')
x['rule'] = dec
x['action'] = np.where(pd.Series(dec).str.startswith('R'), 'remove', 'keep')
cols = ['bin', 'contigName', 'len', 'action', 'rule', 'cov', 'off', 'ratio_adj', 'r', 'tnf', 'gc', 'best', 'species',
        'n_public', 'genes', 'pubfam', 'persistent', 'pf', 'bits_per_Mb', 'bin_bits_per_Mb']
x[cols].to_csv(f'{out}/decisions.tsv', sep='\t', index=False, float_format='%.4g')

# cleaned FASTAs
rm = set(zip(x.bin[x.action == 'remove'], x.contigName[x.action == 'remove']))
per = []
for f in sorted(os.listdir(f'{D}/chesterfield/bins')):
    if not f.endswith('.fna'): continue
    b = f[:-4]; keep = True; kb = rb = 0
    with open(f'{D}/chesterfield/bins/{f}') as fi, open(f'{out}/cleaned/{f}', 'w') as fo:
        for l in fi:
            if l.startswith('>'): keep = (b, l[1:].split()[0]) not in rm
            elif keep: kb += len(l.strip())
            else: rb += len(l.strip())
            if keep: fo.write(l)
    sub = x[x.bin == b]
    row = dict(bin=b, species=sub.species.dropna().iloc[0] if sub.species.notna().any() else '',
               Mb=(kb + rb) / 1e6, removed_Mb=rb / 1e6, removed_frac=rb / max(kb + rb, 1),
               contigs=len(sub), removed_contigs=int((sub.action == 'remove').sum()))
    for r in ('R1_off-population', 'R2_no-reads', 'R3_content+coverage', 'K1_content-only'):
        row[r + '_Mb'] = sub.len[sub.rule == r].sum() / 1e6
    per.append(row)
per = pd.DataFrame(per)
per.to_csv(f'{out}/per_bin.tsv', sep='\t', index=False, float_format='%.4g')
print(x.groupby('rule').agg(contigs=('len', 'size'), Mb=('len', lambda s: s.sum() / 1e6)).round(2).to_string())
lt = pd.cut(x.len, [0, 5e3, 1e4, 2e4, 1e9], labels=['<5k', '5-10k', '10-20k', '>=20k'])
print(pd.crosstab(lt, x.rule, values=x.len / 1e6, aggfunc='sum').round(2).fillna(0).to_string())
print(f'bins changed: {(per.removed_contigs > 0).sum()} of {len(per)}; removed {per.removed_Mb.sum():.2f} of {per.Mb.sum():.2f} Mb')
print(per.sort_values('removed_frac', ascending=False).head(25).round(3).to_string(index=False))

# calibration: coverage classes of known members (>= 2 persistent genes) outside the two contaminated species
mem = x[(x.persistent >= 2) & ~x.species.isin(['tmed14', 'uba9145'])]
print('known members (>= 2 persistent genes, 16 bins outside TMED14/UBA9145), coverage class by length:')
print(pd.crosstab(lt[mem.index], mem['cov'], margins=True).to_string())
