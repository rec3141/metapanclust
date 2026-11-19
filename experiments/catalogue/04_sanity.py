"""Sanity checks on a catalogue: do known multi-gene systems come out as single families, are hubs split?
For each system (Bakta product regex) and each bin with >= MIN_GENES matching genes inside loci:
  units      = number of loci the bin's system genes fall into (1 = one locus)
  top family = the cluster family holding most system genes over all bins; its bin coverage and its
               purity (fraction of its genes that match the system regex)
usage: 04_sanity.py MEMBERS_TSV CATALOGUE_TSV GENES_TSV"""
import re, sys
from collections import Counter, defaultdict
import pandas as pd

mem_tsv, cat_tsv, genes_tsv = sys.argv[1:4]
M = pd.read_csv(mem_tsv, sep='\t', keep_default_na=False)
C = pd.read_csv(cat_tsv, sep='\t', keep_default_na=False)
G = pd.read_csv(genes_tsv, sep='\t', keep_default_na=False, dtype={'idx': int})
SYSTEMS = {
    'ATP synthase F0F1': r'^ATP synthase (subunit|gamma|epsilon|F0|F1)',
    'NADH-quinone oxidoreductase (nuo)': r'NADH-quinone oxidoreductase subunit',
    'ribosomal proteins': r'ribosomal (subunit )?protein',
    'flagellum': r'[Ff]lagell',
    'RNA polymerase core': r'DNA-directed RNA polymerase subunit (alpha|beta|beta\')$',
    'Sec translocon': r'Protein translocase subunit Sec|Sec translocon',
    'ABC transporter (any)': r'ABC transporter|ABC-type',
    'TRAP transporter': r'TRAP transporter',
    'histidine biosynthesis': r'^(Imidazole|Histidinol|Histidine biosynthesis|ATP phosphoribosyltransferase|Phosphoribosyl-(AMP|ATP)|1-\(5-phosphoribosyl\))',
}
loc_of = {}
for r in M.itertuples():
    for i in range(r.first, r.last + 1): loc_of[(r.contig, i)] = (r.locus, r.family)
G['in'] = [(c, i) in loc_of for c, i in zip(G.contig, G.idx)]
fam_genes = defaultdict(list)
for (c, i), (l, f) in loc_of.items(): fam_genes[f].append((c, i))
prod = dict(zip(zip(G.contig, G.idx), G['product']))
fam_nb = dict(zip(C.family, C.n_bins))
print('system\tgenes\tfrac_in_loci\tbins(>=3 genes in loci)\tmedian_loci_per_bin\tbins_one_locus\tfamilies\ttop_family\ttop_bins\ttop_purity\tfamilies_for_80pct')
for name, rx in SYSTEMS.items():
    p = re.compile(rx)
    S = G[G['product'].map(lambda x: bool(p.search(x)))]
    Si = S[S['in']]
    per_bin = Si.groupby('bin').apply(lambda d: len({loc_of[(c, i)][0] for c, i in zip(d.contig, d.idx)}), include_groups=False)
    nbin = Si.groupby('bin').size(); keep = nbin[nbin >= 3].index
    per_bin = per_bin[keep]
    fams = Counter(loc_of[(c, i)][1] for c, i in zip(Si.contig, Si.idx))
    tot = sum(fams.values()); acc = 0; k80 = 0
    for f, n in fams.most_common():
        acc += n; k80 += 1
        if acc >= 0.8 * tot: break
    top, ntop = fams.most_common(1)[0] if fams else ('-', 0)
    pur = sum(bool(p.search(prod[g])) for g in fam_genes[top]) / max(len(fam_genes[top]), 1)
    print(f'{name}\t{len(S)}\t{len(Si) / max(len(S), 1):.2f}\t{len(keep)}\t{per_bin.median() if len(per_bin) else 0:.0f}\t'
          f'{(per_bin == 1).sum()}\t{len(fams)}\t{top}\t{fam_nb.get(top, 0)}\t{pur:.2f}\t{k80}')
