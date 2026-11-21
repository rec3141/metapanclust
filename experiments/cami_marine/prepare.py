"""CAMI II marine: true quality of every submitted bin on the pooled gold-standard assembly, and bin FASTAs.

For each binning (CAMI II submissions on the short-read pooled gold-standard assembly) and each bin of at least
MIN_BP: the majority source genome (by bp), true completeness = bp of that genome in the bin / bp of that genome
in the whole gold-standard assembly, true contamination = bp from other genomes / bin bp (AMBER's purity
definitions). Circular elements (plasmids, unknown) and viruses count as contamination of a genome bin; bins whose
majority source is one of them are skipped.

usage: prepare.py CAMI_DIR OUTDIR [MIN_BP=200000]
writes OUTDIR/bins.tsv and OUTDIR/bins/<binner>__<bin>.fna
"""
import gzip, os, sys, glob
import pandas as pd

cami, out = sys.argv[1:3]
MIN_BP = int(sys.argv[3]) if len(sys.argv) > 3 else 200000
os.makedirs(f'{out}/bins', exist_ok=True)
D = f'{cami}/eval_repo/binning/genome_binning/marine_dataset/data'
gt = pd.read_csv(f'{D}/ground_truth/gsa_pooled_mapping_short.binning', sep='\t', comment='@', header=None,
                 names=['contig', 'genome', 'taxid', 'len'])
meta = pd.read_csv(f'{cami}/simulation_short_read/metadata.tsv', sep='\t').set_index('genome_ID')
gt['cat'] = gt.genome.map(meta.novelty_category)
glen = gt.groupby('genome').len.sum()
c2g = gt.set_index('contig')
rows, members = [], {}
for f in sorted(glob.glob(f'{D}/short_read_pooled_gold_standard_assembly/*.binning')):
    binner = os.path.basename(f).replace('.binning', '')
    b = pd.read_csv(f, sep='\t', comment='@', header=None, usecols=[0, 1], names=['contig', 'bin'])
    b = b.join(c2g, on='contig')
    for bn, d in b.groupby('bin'):
        bp = d.len.sum()
        if bp < MIN_BP: continue
        per = d.groupby('genome').len.sum().sort_values(ascending=False)
        g = per.index[0]
        if meta.novelty_category.get(g) in ('plasmid', 'unknown', 'virus'): continue
        name = f'{binner}__{bn}'
        rows.append(dict(bin=name, binner=binner, bp=int(bp), contigs=len(d), genome=g,
                         category=meta.novelty_category.get(g), true_compl=round(per.iloc[0] / glen[g], 4),
                         true_contam=round(1 - per.iloc[0] / bp, 4), n_genomes=int((per > 0).sum()),
                         second_genome=per.index[1] if len(per) > 1 else '',
                         second_frac=round(per.iloc[1] / bp, 4) if len(per) > 1 else 0.0))
        members[name] = set(d.contig)
T = pd.DataFrame(rows)
T.to_csv(f'{out}/bins.tsv', sep='\t', index=False)
print(T.groupby('binner').agg(bins=('bin', 'size'), compl=('true_compl', 'median'),
                               contam=('true_contam', 'median')).to_string())
# write bin FASTAs in one pass over the assembly
want = {}
for name, cs in members.items():
    for c in cs: want.setdefault(c, []).append(name)
handles = {}
with gzip.open(f'{cami}/simulation_short_read/anonymous_gsa_pooled.fasta.gz', 'rt') as fi:
    cur = None; buf = []

    def flush():
        if cur and cur in want:
            s = ''.join(buf)
            for name in want[cur]:
                with open(f'{out}/bins/{name}.fna', 'a') as fo: fo.write(f'>{cur}\n{s}\n')
    for line in fi:
        if line.startswith('>'):
            flush(); cur = line[1:].split()[0]; buf = []
        else:
            buf.append(line.strip())
    flush()
print(f'{len(T)} bins written')
