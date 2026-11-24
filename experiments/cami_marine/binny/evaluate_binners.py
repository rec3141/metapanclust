"""Score binnings of the CAMI II marine pooled gold-standard assembly against the gold standard.

usage: evaluate_binners.py GOLD.binning SUBMISSION_DIR [EXTRA name=contig2bin.tsv ...] > table.txt
  GOLD.binning      gsa_pooled_mapping_short.binning (contig, genome, taxid, length)
  SUBMISSION_DIR    CAMI II submissions on the short-read pooled gold-standard assembly (*.binning)
  EXTRA             more binnings as contig<TAB>bin tables, e.g. binny=binny_contig2bin.tsv

Per bin: majority genome by bp; completeness = its bp in the bin / its bp in the assembly; contamination = other bp /
bin bp (AMBER's purity complement). Circular elements (plasmids) count as contamination. Reported per binner:
bins >= 50 kb, bins at >= 90% complete & <= 5% contaminated, >= 70/10 and >= 50/10, distinct genomes recovered at
>= 90/5 and >= 50/10, the share of assembly bp in bins, and the share assigned to the bin's majority genome.
"""
import glob, os, sys
import pandas as pd

gold, subdir = sys.argv[1:3]
extra = dict(a.split('=', 1) for a in sys.argv[3:])
gt = pd.read_csv(gold, sep='\t', comment='@', header=None, names=['contig', 'genome', 'taxid', 'len']).set_index('contig')
glen = gt.groupby('genome').len.sum()
total = gt.len.sum()


def load(path, cami):
    if cami:
        return pd.read_csv(path, sep='\t', comment='@', header=None, usecols=[0, 1], names=['contig', 'bin'])
    return pd.read_csv(path, sep='\t', header=None, usecols=[0, 1], names=['contig', 'bin'])


def score(b):
    b = b.drop_duplicates('contig').join(gt, on='contig').dropna(subset=['genome'])
    rows = []
    for bn, d in b.groupby('bin'):
        per = d.groupby('genome').len.sum().sort_values(ascending=False)
        bp = per.sum(); g = per.index[0]
        rows.append(dict(bin=bn, bp=bp, genome=g, compl=per.iloc[0] / glen[g], contam=1 - per.iloc[0] / bp))
    R = pd.DataFrame(rows)
    R = R[R.bp >= 50000]
    hq = R[(R.compl >= .9) & (R.contam <= .05)]; mq = R[(R.compl >= .5) & (R.contam <= .1)]
    return dict(bins=len(R), hq_90_5=len(hq), q_70_10=int(((R.compl >= .7) & (R.contam <= .1)).sum()),
                mq_50_10=len(mq), genomes_hq=hq.genome.nunique(), genomes_mq=mq.genome.nunique(),
                binned_bp=round(R.bp.sum() / total, 3),
                correct_bp=round((R.bp * (1 - R.contam)).sum() / total, 3))


rows = []
for f in sorted(glob.glob(f'{subdir}/*.binning')):
    rows.append(dict(binner=os.path.basename(f).replace('.binning', ''), **score(load(f, True))))
for name, f in extra.items():
    rows.append(dict(binner=name, **score(load(f, False))))
T = pd.DataFrame(rows).sort_values('hq_90_5', ascending=False)
pd.set_option('display.width', 200)
print(T.to_string(index=False))
