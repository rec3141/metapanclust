"""Per-contig gene-family content of the main-assembly Chesterfield bins that sit in a species pangenome
(pangenome survey + the two earlier species).

For each gene of our bin (PPanGGOLiN gff): family, partition, and k = number of public genomes carrying the
family (gene_presence_absence.Rtab). Gene surprise = -log2 (k+1)/(n+2) (binning_sim 'present' term, per gene).
Per contig: genes, genes in a public family (k >= 1), persistent genes, bits.
usage: content.py OUT.tsv"""
import csv, math, os, sys
from collections import defaultdict
S = '/data/scratch/metapanclust_data/pangenome_survey'
rows = []
for lst in ('survey_species.tsv', 'survey_small.tsv'):
    for l in open(f'{S}/{lst}'):
        short, spname, d = l.rstrip('\n').split('\t')[:3]
        hdr, *body = [x.rstrip('\n').split('\t') for x in open(f'{d}/gene_presence_absence.Rtab')]
        pub = [i for i, g in enumerate(hdr) if i > 0 and not g.startswith('chesterfield')]
        n = len(pub)
        k = {r[0]: sum(r[i] != '0' for i in pub) for r in body}
        ours = [g for g in hdr[1:] if g.startswith('chesterfield_') and not g.startswith(('chesterfield_marine', 'chesterfield_freshwater'))]
        for g in ours:
            c = defaultdict(lambda: dict(genes=0, pubfam=0, persistent=0, bits=0.0, rare_mod_genes=0))
            for x in open(f'{d}/gff/{g}.gff'):
                if x.startswith('#') or '\tCDS\t' not in x: continue
                f = x.rstrip('\n').split('\t'); a = dict(y.split('=', 1) for y in f[8].split(';') if '=' in y)
                kk = k.get(a.get('family'), 0); e = c[f[0]]
                e['genes'] += 1; e['pubfam'] += kk > 0; e['persistent'] += a.get('partition') == 'persistent'
                e['bits'] += -math.log2((kk + 1) / (n + 2))
            for contig, e in c.items():
                rows.append(dict(contigName=contig, bin='dastool-' + g[len('chesterfield_'):], species=short, n_public=n, **e))
with open(sys.argv[1], 'w') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)
print(len(rows), 'contigs,', len({r['bin'] for r in rows}), 'bins:', sorted({(r['species'], r['bin']) for r in rows}))
