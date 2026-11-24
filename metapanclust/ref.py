"""Species reference from a PPanGGOLiN / PanGBank pangenome (.h5).

A reference holds, for one species:
  genomes      genome names
  comp         per-genome completeness (CheckM2 from PanGBank metadata if present, else the fraction of
               persistent families carried, capped at 0.99)
  fams         family names; part: 'P' / 'S' / 'C' (persistent, shell, cloud)
  R            families x genomes presence (bool)
  units        scoring units: every PPanGGOLiN module is one unit (scored as a block), every family outside a
               module is its own unit. unit_members: family indices; unit_span: median span in bp of the unit
               where it lies on one contig (max stop - min start), used by the detection model
  U            units x genomes presence (a module counts as present if >= half its families are)
and writes the reference proteins (every distinct sequence, translated, labelled with its family) for MMseqs2.
"""
import os, json
import numpy as np

CODE = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'   # table 11, TCAG order
_NT = np.full(256, 4, np.int64)
for _c, _v in zip(b'TCAGtcag', (0, 1, 2, 3, 0, 1, 2, 3)): _NT[_c] = _v


def translate(dna):
    s = _NT[np.frombuffer(dna, np.uint8)]
    L = (len(s) // 3) * 3
    c = s[:L].reshape(-1, 3)
    bad = (c == 4).any(1)
    idx = c[:, 0] * 16 + c[:, 1] * 4 + c[:, 2]
    aa = np.frombuffer(CODE.encode(), np.uint8)[np.where(bad, 0, idx)].copy()
    aa[bad] = ord('X')
    p = aa.tobytes().decode()
    if p and p[0] in 'LV': p = 'M' + p[1:]          # alternative start codons
    return p.rstrip('*')


def build(h5, outdir, min_module_frac=0.5, max_genomes=400, seed=0):
    import tables
    import pandas as pd
    os.makedirs(outdir, exist_ok=True)
    H = tables.open_file(h5)
    A = H.root.annotations
    genomes = np.array([x.decode() for x in A.genomes.col('name')])
    n_all = len(genomes)
    if n_all > max_genomes:     # large species: a random subset keeps the tree and context models tractable
        genomes = np.sort(np.random.default_rng(seed).choice(genomes, max_genomes, replace=False))
    gidx = {g: i for i, g in enumerate(genomes)}
    ctg_g = dict(zip(A.contigs.col('ID').tolist(), [gidx.get(x.decode(), -1) for x in A.contigs.col('genome')]))
    genes = pd.DataFrame({'gene': A.genes.col('ID'), 'contig': A.genes.col('contig'),
                          'gd': A.genes.col('genedata_id')})
    genes['g'] = genes.contig.map(ctg_g)
    genes = genes[genes.g >= 0]
    gd = A.genedata
    gdid = gd.col('genedata_id')
    st = np.zeros(gdid.max() + 1, np.int64); sp = np.zeros_like(st)
    st[gdid] = gd.col('start'); sp[gdid] = gd.col('stop')
    genes['start'] = st[genes.gd.values]; genes['stop'] = sp[genes.gd.values]
    gf = pd.DataFrame({'gene': H.root.geneFamilies.col('gene'), 'fam': H.root.geneFamilies.col('geneFam')})
    gs = pd.DataFrame({'gene': A.geneSequences.col('gene'), 'seqid': A.geneSequences.col('seqid')})
    genes = genes.merge(gf, on='gene').merge(gs, on='gene')
    FI = H.root.geneFamiliesInfo
    fams = np.array([x.decode() for x in FI.col('name')])
    part = np.array([x.decode()[:1] for x in FI.col('partition')])
    fpos = {f.encode(): i for i, f in enumerate(fams)}
    genes['f'] = genes.fam.map(fpos).astype(np.int64)
    R = np.zeros((len(fams), len(genomes)), bool)
    R[genes.f.values, genes.g.values] = True
    comp = np.full(len(genomes), np.nan)
    try:
        M = H.root.metadata.genomes.pangbank_wf
        for g, c in zip(M.col('ID'), M.col('checkm2_completeness')):
            if g.decode() in gidx: comp[gidx[g.decode()]] = c / 100
    except Exception:
        pass
    pers = part == 'P'
    pf = R[pers].mean(0) if pers.any() else np.ones(len(genomes))
    comp = np.clip(np.where(np.isnan(comp), pf, comp), 0.3, 0.99)
    # units: modules as blocks, other families alone
    mods = {}
    if '/modules' in H:
        for f, m in zip(H.root.modules.col('geneFam'), H.root.modules.col('module')):
            mods.setdefault(int(m), []).append(fpos[f])
    in_mod = set(x for v in mods.values() for x in v)
    members = [np.array(sorted(v)) for v in mods.values()] + [np.array([f]) for f in range(len(fams))
                                                                if f not in in_mod]
    kind = np.array(['module'] * len(mods) + ['family'] * (len(members) - len(mods)))
    names = np.array([f'module_{m}' for m in mods] + [fams[m[0]] for m in members[len(mods):]])
    U = np.array([R[m].mean(0) >= min_module_frac if len(m) > 1 else R[m[0]] for m in members])
    # spans: per (unit, genome) on one contig; median over genomes; single genes: median gene length
    flen = genes.assign(L=(genes.stop - genes.start + 1)).groupby('f').L.median()
    span = np.zeros(len(members))
    f2u = {}
    for u, m in enumerate(members):
        for f in m: f2u[int(f)] = u
    genes['u'] = genes.f.map(f2u)
    sp_ = genes.groupby(['u', 'g', 'contig']).agg(a=('start', 'min'), b=('stop', 'max'), n=('f', 'nunique'))
    sp_ = sp_.reset_index()
    for u, d in sp_.groupby('u'):
        m = members[u]
        if len(m) == 1:
            span[u] = flen.get(int(m[0]), 1000)
        else:
            full = d[d.n >= max(2, int(np.ceil(min_module_frac * len(m))))]
            span[u] = np.median(full.b - full.a + 1) if len(full) else flen.reindex(m).sum()
    # reference proteins: every distinct sequence, labelled with its family
    sf = genes.drop_duplicates('seqid')[['seqid', 'f']]
    need = dict(zip(sf.seqid.values, sf.f.values))
    S = A.sequences
    with open(f'{outdir}/ref_prot.faa', 'w') as fo:
        for st_ in range(0, S.nrows, 50000):
            blk = S.read(st_, min(st_ + 50000, S.nrows))
            for dna, sid in zip(blk['dna'], blk['seqid']):
                f = need.get(int(sid))
                if f is not None:
                    fo.write(f'>{sid}|{f}\n{translate(dna)}\n')
    ctg_len = dict(zip(A.contigs.col('ID').tolist(), A.contigs.col('length').tolist())) if 'length' in \
        A.contigs.colnames else {}
    gt = genes.sort_values(['g', 'contig', 'start'])
    np.savez_compressed(f'{outdir}/genes.npz', g=gt.g.values.astype(np.int32), contig=gt.contig.values.astype(np.int64),
                        start=gt.start.values, stop=gt.stop.values, f=gt.f.values.astype(np.int32),
                        ctg_id=np.array(list(ctg_len), np.int64), ctg_len=np.array(list(ctg_len.values()), np.int64))
    np.savez_compressed(f'{outdir}/ref.npz', genomes=genomes, comp=comp, fams=fams, part=part, R=R,
                        unit_names=names, unit_kind=kind, unit_span=span, U=U,
                        unit_members=np.array([','.join(map(str, m)) for m in members], dtype=object))
    H.close()
    info = dict(kind='species', map_id=0.8, h5=os.path.abspath(h5), genomes_total=n_all, genomes=len(genomes), families=len(fams), units=len(members),
                modules=len(mods), proteins=len(need))
    json.dump(info, open(f'{outdir}/ref.json', 'w'), indent=1)
    return info


class Ref:
    def __init__(self, d):
        z = np.load(f'{d}/ref.npz', allow_pickle=True)
        self.dir = d
        self.genomes = z['genomes']; self.comp = z['comp']; self.fams = z['fams']; self.part = z['part']
        self.R = z['R']; self.U = z['U']; self.unit_names = z['unit_names']; self.unit_kind = z['unit_kind']
        self.species = z['species'] if 'species' in z else None      # per genome, clade references only
        self.info = json.load(open(f'{d}/ref.json')) if os.path.exists(f'{d}/ref.json') else {}
        self.map_id = float(self.info.get('map_id', 0.8))
        self.unit_span = z['unit_span']
        self.unit_members = [np.array([int(x) for x in s.split(',')]) for s in z['unit_members']]
        self.fam2unit = np.zeros(len(self.fams), np.int64)
        for u, m in enumerate(self.unit_members): self.fam2unit[m] = u

    def drop(self, keep):
        """copy of the reference restricted to genomes[keep] (bool or index array), for held-out tests"""
        import copy
        r = copy.copy(self)
        r.genomes = self.genomes[keep]; r.comp = self.comp[keep]; r.R = self.R[:, keep]; r.U = self.U[:, keep]
        if self.species is not None: r.species = self.species[keep]
        return r
