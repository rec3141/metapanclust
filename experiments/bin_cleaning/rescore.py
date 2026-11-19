"""Module surprise (pangenome survey model) of our main-assembly bins before and after cleaning, without rebuilding
the pangenomes: each bin's genes are taken from its PPanGGOLiN gff, genes on removed contigs are dropped, and module
completion / persistent completeness are recomputed from the remaining families. Public prevalence is unchanged.
Survey model: present at completion >= 0.5; p = (k+1)/(n+2); present -log2 p, absent -log2(1 - p c);
rare-present = present and k/n <= 0.1.
usage: rescore.py OUTDIR   (reads OUTDIR/decisions.tsv)"""
import csv, math, sys
from collections import defaultdict
import pandas as pd
S = '/data/scratch/metapanclust_data/pangenome_survey'
out = sys.argv[1]
dec = pd.read_csv(f'{out}/decisions.tsv', sep='\t')
removed = defaultdict(set)
for b, c in zip(dec.bin[dec.action == 'remove'], dec.contigName[dec.action == 'remove']):
    removed[b].add(c)
bins_out, mods_out = [], []
for lst in ('survey_species.tsv', 'survey_small.tsv'):
    for l in open(f'{S}/{lst}'):
        short, spname, d = l.rstrip('\n').split('\t')[:3]
        hdr = open(f'{d}/gene_presence_absence.Rtab').readline().rstrip('\n').split('\t')[1:]
        pub = [g for g in hdr if not g.startswith('chesterfield')]
        present = defaultdict(set)
        for r in csv.DictReader(open(f'{d}/modules_in_genomes.tsv'), delimiter='\t'):
            if float(r['completion']) >= 0.5: present[r['module_id']].add(r['genome'])
        mfams = defaultdict(set)
        for r in csv.DictReader(open(f'{d}/functional_modules.tsv'), delimiter='\t'):
            mfams[r['module_id']].add(r['family_id'])
        persistent = {x.strip() for x in open(f'{d}/partitions/persistent.txt')}
        ours = [g for g in hdr if g.startswith('chesterfield_') and not g.startswith(('chesterfield_marine', 'chesterfield_freshwater'))]
        for g in ours:
            b = 'dastool-' + g[len('chesterfield_'):]
            fam_all, fam_kept = set(), set()
            for x in open(f'{d}/gff/{g}.gff'):
                if x.startswith('#') or '\tCDS\t' not in x: continue
                f = x.rstrip('\n').split('\t'); a = dict(y.split('=', 1) for y in f[8].split(';') if '=' in y)
                fam_all.add(a['family'])
                if f[0] not in removed[b]: fam_kept.add(a['family'])
            res = {}
            for tag, fams in (('orig', fam_all), ('clean', fam_kept)):
                c = len(fams & persistent) / len(persistent)
                tot = pb = ab = 0.0; npres = nrare = 0; rare = set()
                for m, mf in mfams.items():
                    k = len(present[m] & set(pub)); p = (k + 1) / (len(pub) + 2)
                    if len(mf & fams) / len(mf) >= 0.5:
                        s = -math.log2(p); pb += s; npres += 1
                        if k / len(pub) <= 0.1: nrare += 1; rare.add(m)
                    else:
                        s = -math.log2(1 - p * c); ab += s
                    tot += s
                res[tag] = dict(c=c, tot=tot, pb=pb, ab=ab, npres=npres, nrare=nrare, rare=rare)
            check = sum(1 for m in mfams if g in present[m])
            bins_out.append(dict(species=short, bin=b, public=len(pub), modules=len(mfams), present_ppanggolin=check,
                                 present_orig=res['orig']['npres'], present_clean=res['clean']['npres'],
                                 compl_orig=round(res['orig']['c'], 3), compl_clean=round(res['clean']['c'], 3),
                                 rare_orig=res['orig']['nrare'], rare_clean=res['clean']['nrare'],
                                 bits_orig=round(res['orig']['tot'], 1), bits_clean=round(res['clean']['tot'], 1),
                                 bits_present_orig=round(res['orig']['pb'], 1), bits_present_clean=round(res['clean']['pb'], 1)))
            for m in sorted(res['orig']['rare'] | res['clean']['rare']):
                mods_out.append(dict(species=short, bin=b, module=m, families=len(mfams[m]),
                                     rare_orig=int(m in res['orig']['rare']), rare_clean=int(m in res['clean']['rare'])))
pd.DataFrame(bins_out).to_csv(f'{out}/module_surprise_bins.tsv', sep='\t', index=False)
pd.DataFrame(mods_out).to_csv(f'{out}/module_surprise_modules.tsv', sep='\t', index=False)
print(pd.DataFrame(bins_out).to_string(index=False))
