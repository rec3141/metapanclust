"""Genus references for CAMI II marine bins.

  genus_only  bins whose majority genome's species has no PanGBank pangenome but its GTDB genus does: scored
              against the pooled genus reference (all PanGBank species of the genus, up to MAX_SP)
  loso        bins whose species does have a pangenome: scored against the genus reference with that species
              left out (leave-one-species-out), to measure what a genus reference loses against the species one

usage: genus_refs.py GENOME_SPECIES.tsv BINS.tsv CATALOG.tsv HAVE_IDS.txt OUTDIR [MAX_SP=15]
writes OUTDIR/refs.tsv (ref, ids, names), OUTDIR/bin_refs.tsv (bin, ref, mode), OUTDIR/downloads.tsv
"""
import os, sys
import pandas as pd

gs_f, bins_f, cat_f, have_f, out = sys.argv[1:6]
MAX_SP = int(sys.argv[6]) if len(sys.argv) > 6 else 15
os.makedirs(out, exist_ok=True)
g = pd.read_csv(gs_f, sep='\t')
g['genus'] = g.classification.str.extract(r'g__([^;]*)')[0]
b = pd.read_csv(bins_f, sep='\t')
cat = pd.read_csv(cat_f, sep='\t')
have = set(int(x) for x in open(have_f).read().split())
refs, binrefs, ids_all = {}, [], set()
gi = g.set_index('genome')
for genome, d in b.groupby('genome'):
    if genome not in gi.index or pd.isna(gi.genus[genome]): continue
    gen = gi.genus[genome]; sid = gi.pangbank_id[genome]
    members = cat[cat.genus == gen].sort_values('genomes')
    if pd.notna(sid):
        members = members[members.id != int(sid)]
        name, mode = f'{gen}__minus_{int(sid)}', 'loso'
    else:
        name, mode = f'{gen}__all', 'genus_only'
    members = members.head(MAX_SP)
    if len(members) == 0: continue
    refs[name] = (','.join(map(str, members.id)), ','.join(members.species.str.replace(' ', '_')))
    ids_all |= set(members.id)
    for bn in d.bin: binrefs.append((bn, name, mode))
pd.DataFrame([(k, *v) for k, v in refs.items()], columns=['ref', 'ids', 'names']).to_csv(
    f'{out}/refs.tsv', sep='\t', index=False)
B = pd.DataFrame(binrefs, columns=['bin', 'ref', 'mode'])
B.to_csv(f'{out}/bin_refs.tsv', sep='\t', index=False)
need = cat[cat.id.isin(ids_all) & ~cat.id.isin(have)]
need[['id', 'name', 'genomes', 'md5']].to_csv(f'{out}/downloads.tsv', sep='\t', index=False, header=False)
print(f'{len(refs)} genus references ({sum(r.endswith("__all") for r in refs)} full, '
      f'{sum("__minus_" in r for r in refs)} leave-one-species-out); bins: {B["mode"].value_counts().to_dict()}; '
      f'{len(ids_all)} species pangenomes, {len(need)} to download ({need.genomes.sum()} genomes)')
