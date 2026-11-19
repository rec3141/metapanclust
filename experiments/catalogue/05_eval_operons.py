"""Score locus boundaries against curated E. coli transcription units (adjacent same-strand pairs).
A pair is predicted 'same unit' when both genes fall in one locus; reported also with the gap rule.
usage: 05_eval_operons.py PAIRS_TSV GENES_TSV LOCI_TSV"""
import sys
import pandas as pd
pairs_tsv, genes_tsv, loci_tsv = sys.argv[1:4]
P = pd.read_csv(pairs_tsv, sep='\t')
G = pd.read_csv(genes_tsv, sep='\t', keep_default_na=False)
L = pd.read_csv(loci_tsv, sep='\t', keep_default_na=False)
pos = dict(zip(G.locus_tag, zip(G.contig, G.idx)))
loc = {}
for r in L.itertuples():
    for i in range(r.first, r.last + 1): loc[(r.contig, i)] = r.locus
same = []
for r in P.itertuples():
    a, b = pos.get(r.locus_a), pos.get(r.locus_b)
    la, lb = loc.get(a), loc.get(b)
    same.append(la is not None and la == lb)
P['in_locus'] = same
def rep(name, pred):
    y = P.same_tu == 1
    tp = (pred & y).sum(); fp = (pred & ~y).sum(); fn = (~pred & y).sum()
    print(f'{name}\tacc {(pred == y).mean():.3f}\tprec {tp / max(tp + fp, 1):.2f}\trecall {tp / (tp + fn):.2f}\tn_pred {pred.sum()}')
rep('same locus', P.in_locus)
rep('gap<=50 OR (same locus AND gap<=300)', (P.gap <= 50) | (P.in_locus & (P.gap <= 300)))
far = P.gap > 50
print(f'far pairs (gap>50) in same TU recovered: {(P.in_locus & far & (P.same_tu == 1)).sum()} of {(far & (P.same_tu == 1)).sum()}')
