"""Download BGC Atlas v2 release files listed in the site's manifest.json, verifying SHA-256.

usage: download.py OUTDIR [--include REGEX ...] [--max-gb N] [--list]
  default includes: every db/*.tsv.gz table (not the dump/tar bundles that repeat them), clustering-v3/ (BiG-SLiCE
  runs, for local query mode), bgc_uuid_source.parquet (BGC -> source assembly/genome), READMEs.
  Region sequence files (regions_nuc.fasta, regions_prot.faa, all_regions_gbk.tar.gz) are left out unless included
  explicitly. Files larger than --max-gb (default 60) are skipped and reported. Resumable (curl -C -), and files whose
  checksum already matches are skipped.
Data licence: CC BY 4.0 (BGC Atlas annotations); underlying sequences keep their source licences.
"""
import argparse, hashlib, json, os, re, subprocess, sys

BASE = os.environ.get('BGC_ATLAS_URL', 'https://bgc-atlas.ziemertlab.com/downloads')  # mirror of bgc-atlas.cs.uni-tuebingen.de
ap = argparse.ArgumentParser()
ap.add_argument('outdir')
ap.add_argument('--include', nargs='*', default=[r'^db/.*\.tsv\.gz$', r'^clustering-v3/', r'^bgc_uuid_source\.parquet$',
                                                    r'(?i)readme'])   # db/ bundles (dump, tar) repeat the tables
ap.add_argument('--max-gb', type=float, default=60)
ap.add_argument('--list', action='store_true')
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
man_f = os.path.join(a.outdir, 'manifest.json')
subprocess.run(['curl', '-sSfL', '-m', '120', '-o', man_f, f'{BASE}/manifest.json'], check=True)
man = json.load(open(man_f))
files = man.get('files', man)
if isinstance(files, dict):
    files = [dict(path=k, **(v if isinstance(v, dict) else {'size': v})) for k, v in files.items()]
rows = []
for f in files:
    p = f.get('path') or f.get('name'); size = f.get('size') or f.get('bytes') or 0
    sha = f.get('sha256') or f.get('checksum')
    rows.append((p, int(size), sha))
print(f"release {man.get('release') or man.get('id')}: {len(rows)} files, {sum(r[1] for r in rows) / 1e9:.1f} GB")
for p, s, _ in sorted(rows):
    print(f'  {s / 1e9:9.3f} GB  {p}')
if a.list:
    sys.exit(0)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 24), b''): h.update(b)
    return h.hexdigest()


for p, s, sha in rows:
    if not any(re.search(rx, p) for rx in a.include):
        continue
    if s > a.max_gb * 1e9:
        print(f'SKIP (>{a.max_gb} GB): {p} {s / 1e9:.1f} GB'); continue
    dst = os.path.join(a.outdir, p)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.exists(dst) and sha and sha256(dst) == sha:
        continue
    print(f'GET {p} ({s / 1e9:.2f} GB)', flush=True)
    subprocess.run(['curl', '-sSfL', '--retry', '5', '-C', '-', '-o', dst, f'{BASE}/{p}'], check=True)
    if sha and sha256(dst) != sha:
        print(f'CHECKSUM MISMATCH {p}', flush=True)
print('DONE')
