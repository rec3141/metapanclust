"""What is at the contig end where a BGC is cut? BGC Atlas R2.0 regions, fetched one by one with HTTP range requests.

usage: breakage.py ATLAS_DIR OUTDIR [N_PER_GROUP=2000] [THREADS=4]

Groups (families in the seed tier that have both a complete isolate member and a contig-edge MAG member):
  mag_edge      MAG regions at a contig edge
  iso_edge      isolate regions at a contig edge
  iso_complete  isolate regions not at an edge (their region ends are antiSMASH's flank cut-offs, the baseline)
For every region end that is a contig end (and for both ends of complete regions), looks at features within WIN bp
of that end:
  transposase   a PFAM_domain from a transposase / IS / integrase / recombinase Pfam family (TNP set below), or a
                Pfam description matching transpos|insertion|integrase|recombinase
  module        an aSDomain of a modular NRPS/PKS type (MOD set below): breaks inside repeated modules
  partial_cds   the terminal CDS is truncated at the edge (location with < or >)
Also records region length and the contig length (bgc_features). Writes OUTDIR/regions/*.gbk (cache),
OUTDIR/ends.tsv (one row per examined end) and a summary on stdout.
"""
import gzip, io, os, re, sys, tarfile, random, urllib.request
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd

A, out = sys.argv[1:3]
NPG = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
TH = int(sys.argv[4]) if len(sys.argv) > 4 else 4
WIN = 3000
URL = 'https://bgc-atlas.ziemertlab.com/downloads/all_regions_gbk.tar.gz'
TNP = {'PF01609', 'PF01527', 'PF00665', 'PF13683', 'PF13333', 'PF02371', 'PF01548', 'PF03050', 'PF05717',
       'PF04754', 'PF12762', 'PF13610', 'PF01385', 'PF12323', 'PF07282', 'PF13546', 'PF13586', 'PF13751',
       'PF00872', 'PF01710', 'PF13358', 'PF00589', 'PF02899', 'PF13700', 'PF04986', 'PF13701', 'PF05598',
       'PF14690', 'PF01526', 'PF13737', 'PF13817', 'PF03400', 'PF02914', 'PF00239', 'PF07508', 'PF12784',
       'PF13560', 'PF13007', 'PF03184', 'PF01797', 'PF13102', 'PF14659', 'PF13495', 'PF13356', 'PF09588'}
TNP_DESC = re.compile(r'transpos|insertion element|integrase|recombinase|resolvase', re.I)
MOD = {'AMP-binding', 'Condensation', 'Condensation_LCL', 'Condensation_DCL', 'Condensation_Starter',
       'Condensation_Dual', 'Epimerization', 'PCP', 'ACP', 'PKS_KS', 'PKS_AT', 'PKS_KR', 'PKS_DH', 'PKS_ER',
       'PKS_DH2', 'PKS_DHt', 'Thioesterase', 'cMT', 'nMT', 'oMT', 'Heterocyclization', 'ACP_beta', 'PP-binding'}
os.makedirs(f'{out}/regions', exist_ok=True)

X = pd.read_parquet(f'{A}/derived/genome_gcf.parquet',
                    columns=['bgc_id', 'gcf_id', 'source_type', 'contig_edge', 'is_complete', 'contig_length',
                             'start_pos', 'end_pos', 'region_len', 'product_category'])
X['edge'] = X.contig_edge.astype(str).isin(['t', 'True', 'true'])
gcfs = pd.read_csv(f'{A}/db/gcfs.tsv.gz', sep='\t', usecols=['id', 'family_class'])
seed = set(gcfs.id[gcfs.family_class != 'fragmented'])
X = X[X.gcf_id.isin(seed)]
ok_f = set(X[(X.source_type == 'isolate') & ~X.edge].gcf_id) & set(X[(X.source_type == 'MAG') & X.edge].gcf_id)
X = X[X.gcf_id.isin(ok_f)]
rng = random.Random(1)


def sample(d, n):
    # at most 3 per family so a few huge families do not dominate
    d = d.groupby('gcf_id', group_keys=False).apply(lambda g: g.sample(min(3, len(g)), random_state=1))
    return d.sample(min(n, len(d)), random_state=1)


S = pd.concat([sample(X[(X.source_type == 'MAG') & X.edge], NPG).assign(group='mag_edge'),
               sample(X[(X.source_type == 'isolate') & X.edge], NPG).assign(group='iso_edge'),
               sample(X[(X.source_type == 'isolate') & ~X.edge], NPG).assign(group='iso_complete')])
print(f'{len(ok_f):,} families eligible; sampled {S.group.value_counts().to_dict()}', flush=True)

need = set(S.bgc_id)
idx = {}
with gzip.open(f'{A}/all_regions_gbk.index.tsv.gz', 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        name, off, ln, _ = line.rstrip('\n').split('\t')
        b = name[:-4]
        if b in need: idx[b] = (int(off), int(ln))


def fetch(b):
    p = f'{out}/regions/{b}.gbk'
    if os.path.exists(p): return True
    if b not in idx: return False
    off, ln = idx[b]
    for attempt in range(4):
        try:
            req = urllib.request.Request(URL, headers={'Range': f'bytes={off}-{off + ln - 1}'})
            raw = urllib.request.urlopen(req, timeout=120).read()
            with tarfile.open(fileobj=io.BytesIO(gzip.decompress(raw))) as t:
                m = t.getmembers()[0]
                data = t.extractfile(m).read()
            open(p + '.tmp', 'wb').write(data); os.replace(p + '.tmp', p)
            return True
        except Exception:
            continue
    return False


with ThreadPoolExecutor(TH) as ex:
    got = list(ex.map(fetch, S.bgc_id))
print(f'fetched {sum(got)}/{len(S)}', flush=True)

LOC = re.compile(r'(complement\()?(?:join\()?<?(\d+)\.\.>?(\d+)')


def features(path):
    """minimal GenBank feature parser: (type, start, end, partial_left, partial_right, qualifiers text)"""
    feats, cur = [], None
    in_feat = False
    for line in open(path, errors='replace'):
        if line.startswith('FEATURES'):
            in_feat = True; continue
        if line.startswith('ORIGIN'):
            break
        if not in_feat: continue
        if len(line) > 5 and line[5] != ' ':
            if cur: feats.append(cur)
            typ = line[5:21].strip(); loc = line[21:].strip()
            m = LOC.search(loc)
            if not m: cur = None; continue
            cur = [typ, int(m.group(2)), int(m.group(3)), '<' in loc, '>' in loc, '']
        elif cur is not None:
            cur[5] += line.strip() + ' '
    if cur: feats.append(cur)
    return feats


rows = []
for r, ok in zip(S.itertuples(), got):
    if not ok: continue
    F = features(f'{out}/regions/{r.bgc_id}.gbk')
    reglen = max([f[2] for f in F], default=0)
    left_edge = r.edge and r.start_pos <= 1
    right_edge = r.edge and r.contig_length and r.end_pos >= r.contig_length - 1
    ends = []
    if r.group == 'iso_complete':
        ends = [('left', 'region_end'), ('right', 'region_end')]
    else:
        if left_edge: ends.append(('left', 'contig_end'))
        if right_edge: ends.append(('right', 'contig_end'))
    for side, kind in ends:
        near = [f for f in F if (f[1] <= WIN if side == 'left' else f[2] >= reglen - WIN)]
        cds = sorted([f for f in F if f[0] == 'CDS'], key=lambda f: f[1])
        term = (cds[0] if side == 'left' else cds[-1]) if cds else None
        pf = [f for f in near if f[0] == 'PFAM_domain']
        tnp = any(any(acc in f[5] for acc in TNP) or bool(TNP_DESC.search(f[5])) for f in pf)
        dom = [re.search(r'/aSDomain="([^"]+)"', f[5]) for f in near if f[0] == 'aSDomain']
        mod = any(d and d.group(1) in MOD for d in dom)
        rows.append(dict(bgc_id=r.bgc_id, group=r.group, kind=kind, side=side, product=r.product_category,
                         region_len=reglen, contig_len=r.contig_length, transposase=tnp, module=mod,
                         partial_cds=bool(term and (term[3] if side == 'left' else term[4])),
                         n_pfam_near=len(pf)))
E = pd.DataFrame(rows)
E.to_csv(f'{out}/ends.tsv', sep='\t', index=False)
# interior baseline: transposase Pfams anywhere in the region (per region)
pd.set_option('display.width', 200)
E['short_contig'] = E.contig_len < 20000
print('\nper examined end (within 3 kb of the end):')
print(E.groupby(['group', 'kind']).agg(ends=('bgc_id', 'size'), transposase=('transposase', 'mean'),
      module=('module', 'mean'), partial_cds=('partial_cds', 'mean'), contig_lt_20kb=('short_contig', 'mean'),
      median_contig_kb=('contig_len', lambda v: np.nanmedian(v) / 1e3)).round(3).to_string())
E['prod'] = E['product'].astype(str).str.strip('{}').str.split(',').str[0]
print('\nby product class, contig-end ends only:')
print(E[E.kind == 'contig_end'].groupby(['group', 'prod']).agg(ends=('bgc_id', 'size'),
      transposase=('transposase', 'mean'), module=('module', 'mean')).query('ends >= 30').round(3).to_string())
