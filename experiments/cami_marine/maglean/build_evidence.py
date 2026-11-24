"""Build MAGLean evidence JSONs for a sample of CAMI II marine bins.

usage: build_evidence.py CAMI_DIR OUTDIR KRAKEN_DB INSPECT_TXT CHECKM_DATA [N_PER_BINNER=20] [THREADS=12]
  (kraken2 must be on PATH; pyrodigal and pyhmmer importable)

Per bin JSON (MAGLean Kernel format):
  marker_sets   CheckM v1 domain-level collocated marker sets (Bacteria: 58 sets / 104 markers; Archaea: 108 / 150),
                domain chosen by the bp-majority Kraken2 domain of the bin's contigs (Bacteria if none)
  contigs       id, length, markers (CheckM marker accessions, version stripped, found by pyhmmer at the HMMs' own
                trusted cutoffs (--cut_tc), as CheckM does), taxonomy (Kraken2 standard DB lineage, ranks
                domain..species; unclassified contigs have no ranks)
  rna           not available for these bins: all false / 0 (so the MIMAG tier never reaches HQ; read tiers from
                completeness and contamination instead)
"""
import glob, json, os, random, subprocess, sys


def _s(x):
    return x.decode() if isinstance(x, bytes) else x
import pandas as pd

cami, out, kdb, inspect_txt, checkm = sys.argv[1:6]
N = int(sys.argv[6]) if len(sys.argv) > 6 else 20
T = int(sys.argv[7]) if len(sys.argv) > 7 else 12
os.makedirs(f'{out}/evidence', exist_ok=True)
bins = pd.read_csv(f'{cami}/prep/bins.tsv', sep='\t')
sel = bins.groupby('binner', group_keys=False).apply(lambda d: d.sample(min(N, len(d)), random_state=1))
sel.to_csv(f'{out}/selected_bins.tsv', sep='\t', index=False)

# contigs of the selected bins
members = {}
seqs = {}
for b in sel.bin:
    cur = None
    members[b] = []
    for line in open(f'{cami}/prep/bins/{b}.fna'):
        if line.startswith('>'):
            cur = line[1:].split()[0]; members[b].append(cur)
            if cur in seqs: cur = None
            else: seqs[cur] = []
        elif cur is not None:
            seqs[cur].append(line.strip())
seqs = {k: ''.join(v) for k, v in seqs.items()}
fa = f'{out}/contigs.fna'
if not os.path.exists(fa):
    with open(fa, 'w') as fo:
        for k, v in seqs.items(): fo.write(f'>{k}\n{v}\n')
print(f'{len(sel)} bins, {len(seqs)} contigs, {sum(map(len, seqs.values())) / 1e6:.0f} Mb', flush=True)

# Kraken2 lineage per contig
kout = f'{out}/kraken.out'
if not os.path.exists(kout):
    subprocess.run(f'kraken2 --db {kdb} --threads {T} --output {kout}.tmp {fa} > /dev/null 2> {out}/kraken.log',
                   shell=True, check=True)
    os.replace(f'{kout}.tmp', kout)
RANKS = {'D': 'domain', 'P': 'phylum', 'C': 'class', 'O': 'order', 'F': 'family', 'G': 'genus', 'S': 'species'}
lineage, stack = {}, []
for line in open(inspect_txt):
    if line.startswith('#'): continue
    p = line.rstrip('\n').split('\t')
    if len(p) < 6: continue
    code, taxid, name = p[3], int(p[4]), p[5]
    depth = (len(name) - len(name.lstrip(' '))) // 2
    name = name.strip()
    if name in ('Bacteria', 'Archaea', 'Eukaryota'):   # domains carry code D or R2 depending on the taxonomy version
        code = 'D'
    stack = stack[:depth] + [(code, name)]
    lineage[taxid] = {RANKS[c]: n for c, n in stack if c in RANKS}
ctax = {}
for line in open(kout):
    p = line.split('\t')
    if p[0] == 'C': ctax[p[1]] = lineage.get(int(p[2]), {})

# CheckM domain-level marker sets and HMMs
ms = pd.read_csv(f'{checkm}/taxon_marker_sets.tsv', sep='\t', header=None)
dom = {}
for _, r in ms[ms[0] == 'domain'].iterrows():
    sets = [sorted({m.split('.')[0] for m in s}) for s in eval(r[6])]
    dom[r[1]] = sets
need = {m for s in dom['Bacteria'] + dom['Archaea'] for m in s}
mk = f'{out}/contig_markers.json'
if not os.path.exists(mk):
    import pyhmmer, pyrodigal
    from concurrent.futures import ThreadPoolExecutor
    with pyhmmer.plan7.HMMFile(f'{checkm}/hmms/checkm.hmm') as fh:
        hmms = [h for h in fh if _s(h.accession or h.name).split('.')[0] in need]
    finder = pyrodigal.GeneFinder(meta=True)

    def genes(item):
        c, s = item
        return [(f'{c}|{k}', g.translate().rstrip('*')) for k, g in enumerate(finder.find_genes(s.encode()))]
    prots = []
    with ThreadPoolExecutor(T) as ex:
        for r in ex.map(genes, seqs.items(), chunksize=64): prots.extend(r)
    print(f'{len(prots)} genes; {len(hmms)} marker HMMs', flush=True)
    alph = pyhmmer.easel.Alphabet.amino()
    block = pyhmmer.easel.TextSequenceBlock(pyhmmer.easel.TextSequence(name=n.encode(), sequence=p)
                                            for n, p in prots if p).digitize(alph)
    best = {}
    for top in pyhmmer.hmmsearch(hmms, block, cpus=T, bit_cutoffs='trusted'):
        acc = _s(top.query.accession or top.query.name).split('.')[0]
        for h in top.reported:
            if not h.included: continue
            g = _s(h.name)
            if g not in best or h.score > best[g][1]: best[g] = (acc, h.score)
    cm = {}
    for g, (acc, _) in best.items():
        cm.setdefault(g.split('|')[0], []).append(acc)
    json.dump(cm, open(mk, 'w'))
cm = json.load(open(mk))

for _, r in sel.iterrows():
    cs = members[r.bin]
    bp = {}
    for c in cs:
        d = ctax.get(c, {}).get('domain')
        if d: bp[d] = bp.get(d, 0) + len(seqs[c])
    domain = max(bp, key=bp.get) if bp else 'Bacteria'
    domain = domain if domain in dom else 'Bacteria'
    ev = dict(mag_id=r.bin, lineage=domain, marker_sets=dom[domain],
              rna=dict(has_23S=False, has_16S=False, has_5S=False, distinct_tRNAs=0),
              contigs=[dict(id=c, length=len(seqs[c]), markers=cm.get(c, []), taxonomy=ctax.get(c, {})) for c in cs])
    json.dump(ev, open(f'{out}/evidence/{r.bin}.json', 'w'))
print('evidence written', flush=True)
