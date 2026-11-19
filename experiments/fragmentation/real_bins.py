"""(3) Real MAGs: how often do Spacedust conserved-cluster instances in the 93 Chesterfield bins touch
a contig end (include the contig's first or last CDS), against
  null  - clusters placed uniformly within their contig, never truncated: 2 / (n - w + 1)
  model - random assembly breaks (simulate.py: six reference genomes cut with this bin's contigs)."""
import numpy as np, pandas as pd
from collections import defaultdict
import common as C

cds, c2b, ncds = {}, {}, {}
for b in C.bin_names():
    a, _ = C.gff_cds(f'{C.BINS}/{b}.gff3')
    for c, x in a.items():
        cds[c] = {(int(s), int(e)): i for i, (s, e) in enumerate(x)}
        c2b[c] = b; ncds[c] = len(x)
L = {b: C.fasta_lengths(f'{C.BINS}/{b}.fna') for b in C.bin_names()}
inst = C.sd_instances(f'{C.D}/chesterfield/sd_bins_mmseqs/clusters.tsv')
rows = []
for c, mem in inst:
    if c not in cds:
        continue
    ix = sorted(cds[c][m] for m in mem if m in cds[c])
    if len(ix) < 2:
        continue
    n, w = ncds[c], ix[-1] - ix[0] + 1
    touch = ix[0] == 0 or ix[-1] == n - 1
    null = 1.0 if n == w else min(1.0, 2 / (n - w + 1))
    span = max(e for s, e in mem) - min(s for s, e in mem) + 1
    rows.append((c2b[c], c, L[c2b[c]][c], n, len(ix), w, span, touch, null))
r = pd.DataFrame(rows, columns=['bin', 'contig', 'contig_len', 'contig_cds', 'k', 'width', 'span', 'touch', 'null'])
q = C.checkm2()
binfo = pd.DataFrame([(b, C.n50(list(L[b].values())), q[b][0], len(L[b])) for b in L], columns=['bin', 'N50', 'completeness', 'contigs'])
s = pd.read_csv(f'{C.OUT}/sim_clusters.tsv.gz', sep='\t')
s = s[s.set == 'sd6'].groupby('bin').agg(sim_touch=('end_touch', 'sum'), sim_det=('detected', 'sum'))
binfo['model'] = (s.sim_touch / s.sim_det).reindex(binfo.bin).values
per = r.groupby('bin').agg(instances=('touch', 'size'), observed=('touch', 'mean'), null=('null', 'mean'), median_span=('span', 'median'))
binfo = binfo.merge(per, left_on='bin', right_index=True)
binfo['n50b'] = pd.cut(binfo.N50, [0, 10e3, 20e3, 50e3, 100e3, 250e3, 1e7], labels=['<10k', '10-20k', '20-50k', '50-100k', '100-250k', '>250k'])
pd.set_option('display.width', 200)
print(f'{len(r)} cluster instances in {r.bin.nunique()} bins; median span {r.span.median():.0f} bp, median genes {r.k.median():.0f}')
print(f'overall touching a contig end: observed {r.touch.mean():.3f}, null {r.null.mean():.3f}')
w = binfo.instances
agg = binfo.groupby('n50b', observed=True).apply(lambda g: pd.Series({
    'bins': len(g), 'instances': g.instances.sum(),
    'observed': np.average(g.observed, weights=g.instances), 'null_no_trunc': np.average(g.null, weights=g.instances),
    'model_random_breaks': np.average(g.model, weights=g.instances)}))
print(agg.round(3))
print('\nper-bin correlation observed vs model r=%.3f, observed vs null r=%.3f' % (
    np.corrcoef(binfo.observed, binfo.model)[0, 1], np.corrcoef(binfo.observed, binfo.null)[0, 1]))
print('weighted ratio observed / model = %.3f; (observed - null) / (model - null) = %.3f' % (
    np.average(binfo.observed, weights=w) / np.average(binfo.model, weights=w),
    (np.average(binfo.observed, weights=w) - np.average(binfo.null, weights=w)) / (np.average(binfo.model, weights=w) - np.average(binfo.null, weights=w))))
# by contig length class: where are the excess end-touches?
r['clen'] = pd.cut(r.contig_len, [0, 5e3, 10e3, 20e3, 50e3, 100e3, 1e8], labels=['<5k', '5-10k', '10-20k', '20-50k', '50-100k', '>100k'])
print('\nby contig length:')
print(r.groupby('clen', observed=True).agg(instances=('touch', 'size'), observed=('touch', 'mean'), null=('null', 'mean')).round(3))
r.to_csv(f'{C.OUT}/real_instances.tsv.gz', sep='\t', index=False)
binfo.to_csv(f'{C.OUT}/real_bins.tsv', sep='\t', index=False)
