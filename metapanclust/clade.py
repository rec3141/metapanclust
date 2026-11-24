"""Clade (genus) reference: pooled pangenomes of the PanGBank species of a genus, for bins whose own species
has no pangenome.

Each member species contributes up to max(min_per_species, max_genomes / n_species) genomes. Family
representative proteins of all species are clustered into orthogroups (OGs; MMseqs2, 50% identity, 80% coverage
both ways), and a genome carries an OG if it carries any member family. Query genes are later mapped to the
representatives at the same 50% / 80% cut, so "novel" means no homolog anywhere in the genus at that level.

The resulting reference has the same layout as a species reference (ref.npz, ref_prot.faa, ref.json), with:
  units    one per OG (PanGBank modules are species-specific, so they are not carried over)
  part     'P' if the OG is in >= 90% of the sampled genomes of >= 80% of member species (genus core), else 'S';
           completeness is estimated from these
  species  member species per genome; the novelty baseline is the number of species-private OGs per genome,
           which is what a new species of the genus is expected to bring
"""
import os, json, subprocess
import numpy as np, pandas as pd


def read_species(h5, n_max, rng):
    import tables
    H = tables.open_file(h5)
    A = H.root.annotations
    genomes = np.array([x.decode() for x in A.genomes.col('name')])
    keep = np.sort(rng.choice(len(genomes), min(n_max, len(genomes)), replace=False))
    gidx = {genomes[k]: i for i, k in enumerate(keep)}
    ctg_g = dict(zip(A.contigs.col('ID').tolist(), [gidx.get(x.decode(), -1) for x in A.contigs.col('genome')]))
    genes = pd.DataFrame({'gene': A.genes.col('ID'), 'contig': A.genes.col('contig')})
    genes['g'] = genes.contig.map(ctg_g)
    genes = genes[genes.g >= 0]
    gf = pd.DataFrame({'gene': H.root.geneFamilies.col('gene'), 'fam': H.root.geneFamilies.col('geneFam')})
    genes = genes.merge(gf, on='gene')
    FI = H.root.geneFamiliesInfo
    fams = FI.col('name'); prot = FI.col('protein')
    fpos = {f: i for i, f in enumerate(fams)}
    genes['f'] = genes.fam.map(fpos)
    pres = genes[['f', 'g']].drop_duplicates()
    comp = np.full(len(keep), np.nan)
    try:
        M = H.root.metadata.genomes.pangbank_wf
        for g, c in zip(M.col('ID'), M.col('checkm2_completeness')):
            i = gidx.get(g.decode())
            if i is not None: comp[i] = c / 100
    except Exception:
        pass
    H.close()
    used = np.unique(pres.f.values)
    return dict(genomes=genomes[keep], comp=comp, pres=pres, used=used,
                prot={int(i): prot[i].decode().rstrip('*') for i in used})


def build_clade(h5s, species_names, outdir, max_genomes=400, min_per_species=5, map_id=0.5, threads=4, seed=0):
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(seed)
    n_per = max(min_per_species, max_genomes // max(len(h5s), 1))
    S = [read_species(h, n_per, rng) for h in h5s]
    # orthogroups over all member families
    with open(f'{outdir}/fam_reps.faa', 'w') as fo:
        for si, s in enumerate(S):
            for f, p in s['prot'].items(): fo.write(f'>{si}_{f}\n{p}\n')
    subprocess.run(f'mmseqs easy-cluster {outdir}/fam_reps.faa {outdir}/og {outdir}/tmp --min-seq-id {map_id} '
                   f'-c 0.8 --cov-mode 0 --threads {threads} -v 1', shell=True, check=True)
    cl = pd.read_csv(f'{outdir}/og_cluster.tsv', sep='\t', header=None, names=['rep', 'mem'])
    reps = pd.unique(cl.rep); rid = {r: k for k, r in enumerate(reps)}
    og_of = dict(zip(cl.mem, cl.rep.map(rid)))
    nog = len(reps)
    offs = np.cumsum([0] + [len(s['genomes']) for s in S])
    N = offs[-1]
    R = np.zeros((nog, N), bool)
    species = np.empty(N, object); genomes = np.empty(N, object); comp = np.zeros(N)
    plen = np.zeros(nog); pcount = np.zeros(nog)
    for si, s in enumerate(S):
        og = np.array([og_of[f'{si}_{f}'] for f in s['pres'].f.values])
        R[og, offs[si] + s['pres'].g.values] = True
        species[offs[si]:offs[si + 1]] = species_names[si]
        genomes[offs[si]:offs[si + 1]] = s['genomes']
        for f, p in s['prot'].items():
            o = og_of[f'{si}_{f}']; plen[o] += len(p); pcount[o] += 1
    # completeness per genome: CheckM2 if present, else share of its own species' persistent-like families
    for si, s in enumerate(S):
        sl = slice(offs[si], offs[si + 1])
        c = s['comp']
        core = R[:, sl].mean(1) >= 0.9
        pf = R[core, sl].mean(0) if core.any() else np.ones(len(c))
        comp[sl] = np.clip(np.where(np.isnan(c), pf, c), 0.3, 0.99)
    cores = np.stack([R[:, offs[i]:offs[i + 1]].mean(1) >= 0.9 for i in range(len(S))], 1)
    part = np.where(cores.mean(1) >= 0.8, 'P', 'S')        # core in >= 80% of member species
    keep = R.any(1)
    R = R[keep]; part = part[keep]
    newid = -np.ones(nog, np.int64); newid[keep] = np.arange(keep.sum())
    span = (plen[keep] / np.maximum(pcount[keep], 1)) * 3
    with open(f'{outdir}/ref_prot.faa', 'w') as fo:
        for si, s in enumerate(S):
            for f, p in s['prot'].items():
                o = newid[og_of[f'{si}_{f}']]
                if o >= 0: fo.write(f'>{si}_{f}|{o}\n{p}\n')
    names = np.array([f'OG{i}' for i in range(R.shape[0])])
    np.savez_compressed(f'{outdir}/ref.npz', genomes=genomes.astype(str), comp=comp, fams=names, part=part, R=R,
                        unit_names=names, unit_kind=np.array(['og'] * len(names)), unit_span=span, U=R,
                        unit_members=np.array([str(i) for i in range(len(names))], dtype=object),
                        species=species.astype(str))
    for f in ('fam_reps.faa', 'og_all_seqs.fasta', 'og_rep_seq.fasta'):
        if os.path.exists(f'{outdir}/{f}'): os.remove(f'{outdir}/{f}')
    subprocess.run(f'rm -rf {outdir}/tmp', shell=True)
    info = dict(kind='clade', map_id=map_id, species=list(species_names), h5=[os.path.abspath(h) for h in h5s],
                genomes=int(N), ogs=int(R.shape[0]), core_ogs=int((part == 'P').sum()))
    json.dump(info, open(f'{outdir}/ref.json', 'w'), indent=1)
    return info
