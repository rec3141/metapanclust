"""Download the PanGBank .h5 files listed in downloads.tsv (from select_species.py), sequentially and politely.
Run on a login node (compute nodes have no internet). Skips files already present with the right md5;
resumable. Files are written as H5DIR/<id>.h5.

usage: download.py downloads.tsv H5DIR [PAUSE_SECONDS=1]
"""
import hashlib, os, sys, time, urllib.request

API = "https://pangbank-api.genoscope.cns.fr"
lst, h5dir = sys.argv[1], sys.argv[2]
pause = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
os.makedirs(h5dir, exist_ok=True)


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''): h.update(b)
    return h.hexdigest()


rows = [l.rstrip('\n').split('\t') for l in open(lst) if l.strip()]
tot, fails = 0, []
t0 = time.time()
for i, (pid, name, ng, want) in enumerate(rows):
    dst = f'{h5dir}/{pid}.h5'
    if os.path.exists(dst) and md5(dst) == want: continue
    ok = False
    for attempt in range(3):
        try:
            tmp = dst + '.part'
            with urllib.request.urlopen(f'{API}/pangenomes/{pid}/file', timeout=600) as r, open(tmp, 'wb') as f:
                while True:
                    b = r.read(1 << 22)
                    if not b: break
                    f.write(b)
            got = md5(tmp)
            if got != want: raise IOError(f'md5 {got} != {want}')
            os.replace(tmp, dst); ok = True; break
        except Exception as e:
            print(f'{pid} {name}: attempt {attempt + 1} failed: {e}', file=sys.stderr, flush=True)
            time.sleep(30 * (attempt + 1))
    if not ok: fails.append(pid); continue
    tot += os.path.getsize(dst)
    print(f'{i + 1}/{len(rows)} {pid} {name} {ng} genomes {os.path.getsize(dst) / 1e6:.0f} MB '
          f'(total {tot / 1e9:.1f} GB, {time.time() - t0:.0f}s)', flush=True)
    time.sleep(pause)
print(f'done: {len(rows) - len(fails)}/{len(rows)} present; failed: {" ".join(fails) or "none"}', flush=True)
sys.exit(1 if fails else 0)
