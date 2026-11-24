"""Query genes: Pyrodigal gene calls, then MMseqs2 assignment to reference families.

A gene joins a family if its best hit among the reference proteins has >= 80% identity over >= 80% of both
sequences, the PPanGGOLiN/PanGBank clustering thresholds. Other genes are novel to the species.
"""
import os, subprocess
import numpy as np, pandas as pd


def read_fasta(path):
    import gzip
    op = gzip.open if path.endswith('.gz') else open
    name, seq = None, []
    with op(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None: yield name, ''.join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None: yield name, ''.join(seq)


def call_genes(fasta, faa_out):
    import pyrodigal
    ctgs = list(read_fasta(fasta))
    total = sum(len(s) for _, s in ctgs)
    if total >= 100_000:
        gf = pyrodigal.GeneFinder()
        gf.train(*[s.encode() for _, s in ctgs if len(s) >= 1000] or [s.encode() for _, s in ctgs])
    else:
        gf = pyrodigal.GeneFinder(meta=True)
    rows = []
    with open(faa_out, 'w') as fo:
        for ci, (name, s) in enumerate(ctgs):
            for k, g in enumerate(gf.find_genes(s.encode())):
                gid = f'{ci}|{k}'
                fo.write(f'>{gid}\n{g.translate().rstrip("*")}\n')
                rows.append((gid, name, g.begin, g.end))
    genes = pd.DataFrame(rows, columns=['gid', 'contig', 'start', 'stop'])
    lengths = {n: len(s) for n, s in ctgs}
    return genes, lengths


def ref_db(refdir, threads=4):
    db = f'{refdir}/mmdb/ref'
    if not os.path.exists(db + '.dbtype'):
        os.makedirs(f'{refdir}/mmdb', exist_ok=True)
        subprocess.run(f'mmseqs createdb {refdir}/ref_prot.faa {db} -v 1', shell=True, check=True)
    return db


def map_to_ref(faa, refdir, work, threads=4, min_id=0.8, cov=0.8):
    db = ref_db(refdir, threads)
    os.makedirs(work, exist_ok=True)
    out = f'{work}/hits.m8'
    subprocess.run(f'mmseqs easy-search {faa} {db} {out} {work}/tmp --min-seq-id {min_id} -c {cov} --cov-mode 0 '
                   f'--format-output query,target,pident,bits --threads {threads} -v 1', shell=True, check=True)
    h = pd.read_csv(out, sep='\t', header=None, names=['q', 't', 'pident', 'bits'])
    h = h.sort_values('bits', ascending=False).drop_duplicates('q')
    h['fam'] = h.t.str.split('|').str[1].astype(int)
    return h.set_index('q')


def annotate(fasta, refdir, work, threads=4):
    import json
    os.makedirs(work, exist_ok=True)
    genes, lengths = call_genes(fasta, f'{work}/query.faa')
    info = json.load(open(f'{refdir}/ref.json')) if os.path.exists(f'{refdir}/ref.json') else {}
    h = map_to_ref(f'{work}/query.faa', refdir, work, threads, min_id=float(info.get('map_id', 0.8)))
    genes['fam'] = genes.gid.map(h.fam).fillna(-1).astype(int)
    genes['pident'] = genes.gid.map(h.pident)
    return genes, lengths
