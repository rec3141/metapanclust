"""Per-sample coverage test of accessory content in co-assembly bins (contig-level depths only).

For each binned contig: depth profile across samples vs. its bin's core profile.
  ratio   sum(contig depth) / sum(core depth) over the bin's informative samples (core depth >= MIN_CORE);
          ~1 core-like, <1 carried by a fraction of cells, >1 multi-copy or shared.
  r       Pearson correlation of log10(depth + 0.05) with the core profile over all samples.
  off     share of the contig's depth that falls in samples where the core is essentially absent
          (core depth < 0.05); reads from another organism show up here.
  tnf     Euclidean distance of the contig's canonical tetranucleotide frequencies to the core's
          length-weighted mean, divided by the 95th percentile of core contigs' distances
          (> 1 = composition outside the core).
  best    the bin whose core profile the contig correlates with best.

usage:
  coverage_test.py focal DEPTHS CONTIG2BIN FASTA BAKTA_TSV GENOMAD_DIR BIN PPANGGOLIN_DIR GENOME OUTDIR
      core = contigs carrying PPanGGOLiN persistent families; module contigs = contigs carrying genes of
      modules present in GENOME and in no public genome (genome names not starting with 'chesterfield').
  coverage_test.py allbins DEPTHS CONTIG2BIN FASTA OUTDIR [FOCAL_CONTIGS_TSV ...]
      core of each bin = contigs >= 10 kb (all contigs if fewer than 5); classifies every binned contig.
      Focal-mode contig tables, if given, calibrate the classes on persistent-family contigs (known members).
  ratio_adj = ratio / median ratio of trusted contigs of the same length class (short contigs have lower depth).
"""
import csv, glob, itertools, os, sys
from collections import defaultdict
import numpy as np
import pandas as pd

MIN_CORE = 0.2      # core depth for a sample to count as informative
ABSENT = 0.05       # core depth below this: population essentially absent from the sample
PSEUDO = 0.05


def load_depths(path):
    d = pd.read_csv(path, sep='\t', index_col=0)
    cols = [c for c in d.columns if c.endswith('.bam')]
    x = d[cols].copy()
    x.columns = [c.split('.')[0] for c in cols]
    return x, d['contigLen']


def load_c2b(path):
    return pd.read_csv(path, sep='\t', header=None, index_col=0)[1]


def core_profile(x, contigs):
    return x.loc[list(contigs)].median()


def kmer_index():
    comp = str.maketrans('ACGT', 'TGCA')
    idx, n = {}, 0
    for k in map(''.join, itertools.product('ACGT', repeat=4)):
        rc = k.translate(comp)[::-1]
        if rc in idx: idx[k] = idx[rc]
        else: idx[k] = n; n += 1
    return idx, n


def tnf(fasta, keep):
    """Canonical tetranucleotide frequencies and GC of contigs in keep (vectorised)."""
    idx, n = kmer_index()
    lut = np.full(256, 4, dtype=np.int64)
    for i, b in enumerate('ACGT'): lut[ord(b)] = i; lut[ord(b.lower())] = i
    canon = np.array([idx[''.join(k)] for k in itertools.product('ACGT', repeat=4)])
    out, gc, name, seq = {}, {}, None, []
    def flush():
        if name in keep:
            a = lut[np.frombuffer(''.join(seq).encode(), dtype=np.uint8)]
            ok = (a[:-3] < 4) & (a[1:-2] < 4) & (a[2:-1] < 4) & (a[3:] < 4)
            code = (a[:-3] * 64 + a[1:-2] * 16 + a[2:-1] * 4 + a[3:])[ok]
            v = np.bincount(canon[code], minlength=n).astype(float)
            out[name] = v / max(v.sum(), 1)
            gc[name] = ((a == 1) | (a == 2)).sum() / max((a < 4).sum(), 1)
    for l in open(fasta):
        if l.startswith('>'):
            if name: flush()
            name, seq = l[1:].split()[0], []
        else: seq.append(l.strip())
    if name: flush()
    return pd.DataFrame(out).T, pd.Series(gc)


def metrics(x, contigs, core, tnfm, core_contigs, lens):
    """Per-contig ratio, r, off, tnf distance against one core profile."""
    inf = core[core >= MIN_CORE].index
    absent = core[core < ABSENT].index
    lc = np.log10(core + PSEUDO)
    xs = x.loc[list(contigs)]
    lx = np.log10(xs + PSEUDO)
    r = lx.apply(lambda row: np.corrcoef(row, lc)[0, 1] if row.std() > 0 else np.nan, axis=1)
    ratio = xs[inf].sum(axis=1) / core[inf].sum() if len(inf) else pd.Series(np.nan, index=xs.index)
    off = xs[absent].sum(axis=1) / xs.sum(axis=1).replace(0, np.nan)
    cc = [c for c in core_contigs if c in tnfm.index]
    w = lens.loc[cc]
    cen = (tnfm.loc[cc].mul(w, axis=0)).sum() / w.sum()
    dist = np.sqrt(((tnfm - cen) ** 2).sum(axis=1))
    scale = np.percentile(dist.loc[cc], 95) if len(cc) >= 3 else np.nan
    t = (dist / scale).reindex(xs.index)
    return pd.DataFrame({'len': lens.loc[xs.index], 'ratio': ratio, 'r': r, 'off': off, 'tnf': t,
                         'depth_main': xs[core.idxmax()]}), inf


def best_bin(x, contigs, cores):
    lc = np.log10(cores + PSEUDO)                      # bins x samples
    lc = lc.sub(lc.mean(axis=1), axis=0).div(lc.std(axis=1), axis=0)
    lx = np.log10(x.loc[list(contigs)] + PSEUDO)
    lx = lx.sub(lx.mean(axis=1), axis=0).div(lx.std(axis=1).replace(0, np.nan), axis=0)
    R = lx.values @ lc.values.T / (lx.shape[1] - 1)
    R = pd.DataFrame(R, index=lx.index, columns=cores.index)
    return R


def tnf_best(tnfm, lens, c2b, contigs):
    """Bin whose core composition (length-weighted TNF centroid of contigs >= 10 kb) is closest, in units of
    that bin's 95th-percentile core distance."""
    cen, sc = {}, {}
    for b, cs in c2b.groupby(c2b).groups.items():
        big = [c for c in cs if lens[c] >= 10000 and c in tnfm.index]
        cc = big if len(big) >= 5 else [c for c in cs if c in tnfm.index]
        w = lens.loc[cc]; v = tnfm.loc[cc].mul(w, axis=0).sum() / w.sum()
        cen[b] = v; sc[b] = np.percentile(np.sqrt(((tnfm.loc[cc] - v) ** 2).sum(axis=1)), 95)
    C = pd.DataFrame(cen).T; S = pd.Series(sc)
    T = tnfm.loc[[c for c in contigs if c in tnfm.index]]
    D = np.sqrt(np.maximum((T ** 2).sum(axis=1).values[:, None] + (C ** 2).sum(axis=1).values[None, :]
                           - 2 * T.values @ C.values.T, 0)) / S.values[None, :]
    D = pd.DataFrame(D, index=T.index, columns=C.index)
    return D.idxmin(axis=1), D.min(axis=1)


def summarise(df, label):
    q = lambda s: f'{s.median():.2f} [{s.quantile(.1):.2f}-{s.quantile(.9):.2f}]' if s.notna().any() else 'na'
    return (f'{label:38s} n={len(df):5d} kb={df.len.sum() / 1e3:7.0f}  ratio_adj {q(df.ratio_adj)}  r {q(df.r)}  '
            f'off {q(df.off)}  tnf {q(df.tnf)}')


LEN_EDGES = [0, 3000, 5000, 7500, 10000, 15000, 20000, 50000, 1e12]
TNF_OUT = 1.25      # composition outlier: tnf distance > 1.25 x the core's 95th percentile


def len_class(lens, edges=LEN_EDGES):
    return pd.cut(lens, edges, labels=[f'{int(a / 1000)}-{int(b / 1000) if b < 1e12 else "inf"}k'
                                       for a, b in zip(edges[:-1], edges[1:])])


def length_expect(m, mask, edges=LEN_EDGES):
    """Median ratio of trusted contigs per length class. Contig depths fall with contig length (read ends are
    clipped at contig ends), so short core contigs sit below ratio 1; ratio_adj divides this out."""
    lc = len_class(m.len, edges)
    e = m[mask].groupby(lc[mask], observed=False).ratio.median()
    return lc.map(e).astype(float), e


def classify(m, core_r10, acc_max=0.6, hi_min=2.0):
    """Coverage class from ratio_adj, r and off:
    no-reads; off-population (most depth in samples where the core is absent); other-bin (another bin's core
    profile fits clearly better and composition is outside the own core); otherwise accessory-like
    (ratio_adj < acc_max), core-like, high-copy (> hi_min); '+low-r' when r is below the 10th percentile of the
    bin's core contigs (noisy profile)."""
    lab = pd.Series('core-like', index=m.index)
    lab[m.ratio_adj < acc_max] = 'accessory-like'
    lab[m.ratio_adj > hi_min] = 'high-copy'
    lowr = m.r < m.bin.map(core_r10).clip(upper=0.8)
    lab[lowr] += '+low-r'
    lab[(m.best != m.bin) & (m.r_best >= 0.9) & (m.r_best - m.r > 0.3) & (m.tnf > TNF_OUT)] = 'other-bin'
    lab[m.off > 0.5] = 'off-population'
    lab[m.r.isna() | (m.ratio.isna())] = 'no-reads'
    return lab


def focal(depths, c2b_path, fasta, bakta, gm, bin_name, pdir, genome, out):
    os.makedirs(out, exist_ok=True)
    x, lens = load_depths(depths); c2b = load_c2b(c2b_path)
    binc = set(c2b[c2b == bin_name].index)
    # PPanGGOLiN gene table for our genome
    genes = []
    for l in open(f'{pdir}/gff/{genome}.gff'):
        if l.startswith('#') or '\tCDS\t' not in l: continue
        f = l.rstrip('\n').split('\t'); a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
        genes.append((f[0], int(f[3]), int(f[4]), a['ID'], a.get('family'), a.get('partition'), a.get('module', '')))
    g = pd.DataFrame(genes, columns=['contig', 'start', 'end', 'gene', 'family', 'partition', 'module'])
    present = defaultdict(set)
    for r in csv.DictReader(open(f'{pdir}/modules_in_genomes.tsv'), delimiter='\t'):
        if float(r['completion']) >= 0.5: present[r['module_id']].add(r['genome'])
    surp = {m for m, gs in present.items() if genome in gs and all(s.startswith('chesterfield') for s in gs)}
    fam2mod = {}
    for r in csv.DictReader(open(f'{pdir}/functional_modules.tsv'), delimiter='\t'):
        if r['module_id'] in surp: fam2mod[r['family_id']] = r['module_id']
    g['smod'] = g.family.map(fam2mod).fillna('')
    npub = defaultdict(set)                          # family -> public genomes carrying it
    for l in open(f'{pdir}/gene_families.tsv'):
        f = l.split('\t')
        if not f[1].startswith('chesterfield'): npub[f[0]].add(f[1].rsplit('_CDS_', 1)[0])
    g['n_pub_genomes'] = g.family.map(lambda f: len(npub.get(f, ())))
    prod = {}
    for l in open(bakta):
        if l.startswith('#'): continue
        f = l.rstrip('\n').split('\t')
        if len(f) > 7 and f[1] == 'cds': prod[(f[0], int(f[2]), int(f[3]))] = f[7]
    g['product'] = [prod.get((c, s, e), '') for c, s, e in zip(g.contig, g.start, g.end)]
    mge = {}
    for kind, fn in (('plasmid', 'plasmid_summary.tsv'), ('virus', 'virus_summary.tsv')):
        p = f'{gm}/{fn}'
        if os.path.exists(p):
            for r in csv.DictReader(open(p), delimiter='\t'): mge[r['seq_name'].split('|')[0]] = kind
    g.to_csv(f'{out}/{genome}_genes.tsv', sep='\t', index=False)
    pc = g.groupby('contig').agg(n_genes=('gene', 'size'), n_persist=('partition', lambda s: (s == 'persistent').sum()),
                                 n_pubfam=('n_pub_genomes', lambda s: (s > 0).sum()),
                                 n_smod=('smod', lambda s: (s != '').sum()),
                                 n_flag=('product', lambda s: s.str.contains('flagell', case=False).sum()),
                                 mods=('smod', lambda s: ','.join(sorted(set(s) - {''}))))
    missing = set(pc.index) - binc
    if missing: print(f'warning: {len(missing)} PPanGGOLiN contigs not in {bin_name}', file=sys.stderr)
    core_c = set(pc[(pc.n_persist > 0) & (pc.n_smod == 0)].index)
    mod_c = set(pc[pc.n_smod > 0].index)
    acc_c = set(pc[(pc.n_persist == 0) & (pc.n_smod == 0)].index)
    core = core_profile(x, core_c)
    # cores of all bins (median of contigs >= 10 kb) for the "other bin fits better" test and nulls
    cores = all_cores(x, lens, c2b)
    cores.loc[bin_name] = core
    allc = set(c2b.index)
    tnfm, gc = tnf(fasta, allc)
    m, inf = metrics(x, allc, core, tnfm, core_c, lens)
    m['bin'] = c2b.reindex(m.index); m['gc'] = gc
    R = best_bin(x, allc, cores)
    m['best'] = R.fillna(-2).idxmax(axis=1).where(R.notna().any(axis=1)); m['r_best'] = R.max(axis=1)
    m['tnf_best'], m['tnf_best_d'] = tnf_best(tnfm, lens, c2b, m.index)
    edges = [0, 5000, 10000, 20000, 1e12]
    m['ratio_adj'] = m.ratio / length_expect(m, m.index.isin(core_c), edges)[0]
    core_r10 = pd.Series(m.loc[list(core_c)].r.quantile(0.1), index=[bin_name])
    own = m.loc[sorted(binc)].join(pc)
    own['set'] = np.where(own.index.isin(mod_c), 'module', np.where(own.index.isin(core_c), 'core',
                          np.where(own.index.isin(acc_c), 'other_accessory', 'no_genes')))
    own['mge'] = [mge.get(c, '') for c in own.index]
    lc = np.log10(cores + PSEUDO); lf = lc.loc[bin_name]
    rbin = lc.apply(lambda row: np.corrcoef(row, lf)[0, 1] if row.std() > 0 else np.nan, axis=1).drop(bin_name)
    co = rbin[rbin >= 0.8].index
    others = m[~m.index.isin(binc)]
    contam = others[others.bin.isin(co)]
    own['class'] = classify(own, core_r10)
    own['tnf_out'] = own.tnf > TNF_OUT
    own.to_csv(f'{out}/{genome}_contigs.tsv', sep='\t')
    ct = contam.assign(bin=bin_name)                # co-occurring bins' contigs scored as if in this bin
    ct['class'] = classify(ct, core_r10)
    print(f'== {genome} = {bin_name}: {len(binc)} contigs, {lens.loc[list(binc)].sum() / 1e6:.2f} Mb; '
          f'{len(surp)} surprising modules on {len(mod_c)} contigs')
    print('core profile (informative samples): ' + ', '.join(f'{s} {core[s]:.2f}' for s in inf)
          + f'; samples with core < {ABSENT}: {(core < ABSENT).sum()}/{len(core)}')
    print(f'other bins whose core profile correlates r >= 0.8 with this core: {len(co)} '
          f'({", ".join(f"{b} {rbin[b]:.2f}" for b in rbin.sort_values(ascending=False).index[:5])} ...)')
    print('contig sets (median [10th-90th pct]):')
    print('length expectation (median ratio of persistent-core contigs): '
          + ', '.join(f'{k} {v:.2f}' for k, v in length_expect(m, m.index.isin(core_c), edges)[1].items()))
    for lab, df in (('core (persistent families)', own[own.set == 'core']),
                    ('surprising-module contigs', own[own.set == 'module']),
                    ('  of which flagellar', own[(own.set == 'module') & (own.n_flag > 0)]),
                    ('other accessory contigs (shell/cloud)', own[own.set == 'other_accessory']),
                    ('contigs of all other bins', others),
                    ('contigs of co-occurring bins (r>=0.8)', contam)):
        print('  ' + summarise(df, lab))
    # per-module table
    rows = []
    for mod in sorted(surp, key=lambda s: int(s.split('_')[1])):
        cs = sorted(set(g.loc[g.smod == mod, 'contig']))
        if not cs: continue
        sub = own.loc[cs]; w = sub.len
        prods = g.loc[g.smod == mod, 'product']
        top = prods[prods.str.len() > 0].value_counts().index[:3]
        rows.append(dict(module=mod, genes=(g.smod == mod).sum(), contigs=len(cs),
                         contigs_with_persistent=(sub.n_persist > 0).sum(),
                         contig_genes=int(sub.n_genes.sum()), contig_genes_in_public_fams=int(sub.n_pubfam.sum()),
                         len=sub.len.sum(), ratio=np.average(sub.ratio, weights=w),
                         ratio_adj=np.average(sub.ratio_adj, weights=w), r=np.average(sub.r.fillna(0), weights=w),
                         off=np.average(sub.off.fillna(0), weights=w), tnf=np.average(sub.tnf, weights=w),
                         flagellar=prods.str.contains('flagell', case=False).sum(),
                         cls='/'.join(sub['class']) + ('+tnf-out' if sub.tnf_out.any() else ''), products='; '.join(top)))
    mt = pd.DataFrame(rows)
    mt.to_csv(f'{out}/{genome}_modules.tsv', sep='\t', index=False)
    with pd.option_context('display.width', 250, 'display.max_colwidth', 70):
        print(mt.round(2).to_string(index=False))
    def tab(df):
        k = df['class'] + np.where(df.tnf > TNF_OUT, ' tnf-out', '')
        return ', '.join(f'{c} {n}' for c, n in k.value_counts().items())
    print('classes (calibration and modules):')
    print(f'  persistent-core contigs (true members, n={len((cr := own[own.set == "core"]))}): {tab(cr)}')
    print(f'  co-occurring bins scored as members (contaminant-like, n={len(ct)}): {tab(ct)}')
    print(f'  other accessory contigs (n={(own.set == "other_accessory").sum()}): {tab(own[own.set == "other_accessory"])}')
    print(f'  surprising-module contigs (n={(own.set == "module").sum()}): {tab(own[own.set == "module"])}')
    ok = lambda df: df['class'].isin(['core-like', 'accessory-like']) & (df.tnf <= TNF_OUT)
    print(f'pass (core-like or accessory-like, tnf <= {TNF_OUT}): persistent-core {ok(cr).mean():.2f}; '
          f'co-occurring-bin contigs {ok(ct).mean():.3f} ({ok(ct).sum()}/{len(ct)}); '
          f'surprising-module contigs {ok(own[own.set == "module"]).sum()}/{(own.set == "module").sum()}')
    print(f'persistent-core contigs with ratio_adj < 0.6: {(cr.ratio_adj < 0.6).mean():.3f}')
    print('genes in families shared with public genomes, by class of the contig: ' + ', '.join(
        f'{c} {d.n_pubfam.sum()}/{d.n_genes.sum()}' for c, d in own[own.set != 'core'].groupby('class')))
    return own


def all_cores(x, lens, c2b):
    cores = {}
    for b, cs in c2b.groupby(c2b).groups.items():
        big = [c for c in cs if lens[c] >= 10000]
        cores[b] = core_profile(x, big if len(big) >= 5 else cs)
    return pd.DataFrame(cores).T


def allbins(depths, c2b_path, fasta, out, *truth):
    os.makedirs(out, exist_ok=True)
    x, lens = load_depths(depths); c2b = load_c2b(c2b_path)
    cores = all_cores(x, lens, c2b)
    tnfm, gc = tnf(fasta, set(c2b.index))
    res = []
    for b, cs in c2b.groupby(c2b).groups.items():
        big = [c for c in cs if lens[c] >= 10000]
        m, inf = metrics(x, cs, cores.loc[b], tnfm, big if len(big) >= 5 else list(cs), lens)
        m['bin'] = b; m['n_inf'] = len(inf); res.append(m)
    m = pd.concat(res); m['gc'] = gc
    exp, curve = length_expect(m, (m.off <= 0.5) & m.ratio.notna())
    m['ratio_adj'] = m.ratio / exp
    R = best_bin(x, m.index, cores)
    m['best'] = R.fillna(-2).idxmax(axis=1).where(R.notna().any(axis=1)); m['r_best'] = R.max(axis=1)
    m['tnf_best'], m['tnf_best_d'] = tnf_best(tnfm, lens, c2b, m.index)
    # how distinct is each bin's profile: max correlation with another bin's core
    lc = np.log10(cores + PSEUDO); C = np.corrcoef(lc.values); np.fill_diagonal(C, np.nan)
    m['bin_max_r_other'] = m.bin.map(pd.Series(np.nanmax(C, axis=1), index=cores.index))
    big = m[m.len >= 10000]
    m['class'] = classify(m, big.groupby('bin').r.quantile(0.1).reindex(cores.index).fillna(0.8))
    m['len_class'] = len_class(m.len)
    m['tnf_out'] = m.tnf > TNF_OUT
    m.to_csv(f'{out}/allbins_contigs.tsv', sep='\t')
    tab = m.groupby('class').agg(contigs=('len', 'size'), Mb=('len', lambda s: s.sum() / 1e6))
    tab['share_contigs'] = tab.contigs / len(m); tab['share_bp'] = tab.Mb / tab.Mb.sum()
    print(f'== all bins: {m.bin.nunique()} bins, {len(m)} contigs, {lens.loc[m.index].sum() / 1e6:.1f} Mb')
    print(f'bins with >= 1 informative sample: {(cores.max(axis=1) >= MIN_CORE).sum()}; informative samples per bin: '
          f'{cores.ge(MIN_CORE).sum(axis=1).value_counts().sort_index().to_dict()}')
    print(f'bins whose core correlates r >= 0.8 with another bin core: {(np.nanmax(C, axis=1) >= 0.8).sum()}')
    print('length expectation (median ratio of contigs with off <= 0.5): ' + ', '.join(f'{k} {v:.2f}' for k, v in curve.items()))
    print(tab.round(3).to_string())
    print('composition outliers (tnf > %.2f) by class: ' % TNF_OUT + ', '.join(
        f'{c} {v:.2f}' for c, v in m.groupby('class').tnf_out.mean().items()))
    per = m.groupby(['bin', 'class']).len.sum().unstack(fill_value=0)
    per = per.div(per.sum(axis=1), axis=0)
    per['n_inf'] = m.groupby('bin').n_inf.first()
    per.round(3).to_csv(f'{out}/allbins_per_bin.tsv', sep='\t')
    lt = pd.crosstab(m.len_class, m['class'], values=m.len, aggfunc='sum', normalize='index', dropna=False)
    lt['contigs'] = m.len_class.value_counts()
    print('share of bp by contig length:'); print(lt.round(3).to_string())
    acc = m[m['class'].str.startswith('accessory-like')]
    cl = m[m['class'].str.startswith('core-like')]
    print(f'accessory-like contigs: ratio_adj median {acc.ratio_adj.median():.2f}; tnf > 1: {(acc.tnf > 1).mean():.2f} '
          f'(core-like {(cl.tnf > 1).mean():.2f}); length median {acc.len.median():.0f} (core-like {cl.len.median():.0f})')
    for c in ('accessory-like', 'core-like', 'high-copy', 'off-population', 'other-bin', 'no-reads'):
        sub = m[m['class'].str.startswith(c)]
        print(f'  {c:15s} gc-gap to bin median {(sub.gc - sub.bin.map(m.groupby("bin").gc.median())).abs().median():.3f}  '
              f'tnf median {sub.tnf.median():.2f}')
    if truth:                                       # focal-mode contig tables: persistent contigs = known members
        p = pd.concat([pd.read_csv(t, sep='\t', index_col=0) for t in truth])
        p = [c for c in p[p.set == 'core'].index if c in m.index]
        base = m['class'].str.replace('+low-r', '', regex=False)
        ok = ~base.isin(['off-population', 'no-reads', 'other-bin'])
        fp = (base[p] == 'accessory-like').groupby(m.len_class[p], observed=False).mean()
        obs = (base[ok] == 'accessory-like').groupby(m.len_class[ok], observed=False).mean()
        n = m.len_class[ok].value_counts()
        print(f'calibration on {len(p)} persistent-family contigs of the focal bins (known members): class shares '
              + str(base[p].value_counts(normalize=True).round(3).to_dict()))
        print('accessory-like share by length, all bins (coverage-consistent contigs) vs known members:')
        print(pd.DataFrame({'all_bins': obs, 'known_members': fp, 'contigs': n, 'members': m.len_class[p].value_counts()}).round(3).to_string())
        print(f'overall accessory-like {obs.mul(n).sum() / n.sum():.3f} vs length-matched noise expectation '
              f'{fp.fillna(0).mul(n).sum() / n.sum():.3f}')
    print('per-bin share of bp (median over bins): ' + ', '.join(f'{c} {per[c].median():.3f}' for c in per.columns if c != 'n_inf'))
    return m


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'focal': focal(*sys.argv[2:11])
    elif mode == 'allbins': allbins(*sys.argv[2:])
    else: sys.exit(__doc__)
