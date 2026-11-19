"""Build inputs for catalogue surprise: bin x family observation matrix (I/P/A), family spans,
per-bin contig lengths, completeness, contamination, lineage. Writes npz + tsv to OUT."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/data/dev/metapanclust/experiments/fragmentation')
from common import fasta_lengths

D = '/data/scratch/metapanclust_data'
CAT = f'{D}/catalogue/final'
OUT = f'{D}/catalogue_surprise'

tax = pd.read_csv(f'{D}/catalogue/bin_taxonomy.tsv', sep='\t', header=None, names=['bin', 'lineage'])
bins = sorted(tax.bin)
bi = {b: i for i, b in enumerate(bins)}
cm = pd.read_csv(f'{D}/chesterfield/mag_chesterfield/binning/checkm2/quality_report.tsv', sep='\t').set_index('Name')
mem = pd.read_csv(f'{CAT}/chesterfield_members.tsv', sep='\t')
genes = pd.read_csv(f'{D}/catalogue/fam30/genes.tsv', sep='\t', usecols=['bin', 'contig', 'idx', 'start', 'end'])
g = genes.set_index(['bin', 'contig', 'idx'])
s0 = g.loc[list(zip(mem.bin, mem.contig, mem['first']))]
s1 = g.loc[list(zip(mem.bin, mem.contig, mem['last']))]
mem['span'] = np.maximum(s0.end.values, s1.end.values) - np.minimum(s0.start.values, s1.start.values) + 1
fams = sorted(mem.family.unique())
fi = {f: i for i, f in enumerate(fams)}
# span per family: median over loci not touching a contig end; else the longest truncated locus
sp = {}
for f, d in mem.groupby('family'):
    ok = d[~d.at_contig_end]
    sp[f] = float(ok.span.median()) if len(ok) else float(d.span.max())
span = np.array([sp[f] for f in fams])
# observation: 2 = intact (some locus not at contig end), 1 = partial (all loci at contig ends), 0 = absent
obs = np.zeros((len(bins), len(fams)), dtype=np.int8)
nloc = np.zeros((len(bins), len(fams)), dtype=np.int16)
for (b, f), d in mem.groupby(['bin', 'family']):
    obs[bi[b], fi[f]] = 2 if (~d.at_contig_end).any() else 1
    nloc[bi[b], fi[f]] = len(d)
lens = {b: np.array(list(fasta_lengths(f'{D}/chesterfield/bins/{b}.fna').values())) for b in bins}
meta = pd.DataFrame({'bin': bins, 'lineage': tax.set_index('bin').loc[bins, 'lineage'].values,
                     'compl': cm.loc[bins, 'Completeness'].values / 100,
                     'contam': cm.loc[bins, 'Contamination'].values / 100,
                     'size_bp': [lens[b].sum() for b in bins], 'n_contigs': [len(lens[b]) for b in bins]})
meta.to_csv(f'{OUT}/bins.tsv', sep='\t', index=False)
np.savez_compressed(f'{OUT}/inputs.npz', obs=obs, nloc=nloc, span=span, fams=np.array(fams), bins=np.array(bins),
                    **{f'L_{i}': lens[b] for i, b in enumerate(bins)})
mem[['locus', 'bin', 'contig', 'family', 'span', 'at_contig_end']].to_csv(f'{OUT}/loci_span.tsv', sep='\t', index=False)
print('bins', len(bins), 'families', len(fams), 'present', (obs > 0).sum(), 'intact', (obs == 2).sum(),
      'partial-only', (obs == 1).sum(), 'span median', np.median(span))
