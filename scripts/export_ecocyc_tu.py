"""Export EcoCyc flat files (transunits.dat, genes.dat, ...) from a Pathway Tools SIF, for curated E. coli
operons. Reuses danaSeq's metacyc_scorer.install helpers. Licensed data: written owner-only, never commit.
usage: export_ecocyc_tu.py DANASEQ_DIR SIF WORKDIR OUTDIR"""
import glob, os, shutil, sys
dana, sif, work, out = sys.argv[1:5]
sys.path.insert(0, os.path.join(dana, 'mag_analysis', 'bin'))
from metacyc_scorer import install as I
os.makedirs(f'{work}/home', exist_ok=True)
probe = I._exec(sif, [work], 'ls -d /opt/pathway-tools/aic-export/pgdbs/biocyc/ecocyc/*/data; '
                'ls -d /opt/ptools-local-template /data/ptools-local 2>/dev/null', work)
dirs = [x for x in probe.stdout.split() if x.startswith('/')]
data = next((d for d in dirs if d.endswith('/data') and '/ecocyc/' in d), None)
tmpl = next((d for d in dirs if 'ptools-local' in d), None)
if not data:
    sys.exit(f'no EcoCyc PGDB in {sif}: {probe.stdout[-500:]}')
ec, pl = f'{work}/ecdata', f'{work}/ptools-local'
os.makedirs(ec, exist_ok=True); os.makedirs(pl, exist_ok=True)
if tmpl:
    I._exec(sif, [work], f'cp -a {tmpl}/. {pl}/', work)
I._export(sif, [work, f'{ec}:{data}', f'{pl}:/data/ptools-local'], '/opt/pathway-tools/pathway-tools', pl, work,
          f'{work}/export.log', org='ecoli')
os.makedirs(out, exist_ok=True); os.chmod(out, 0o700)
for name in ('transunits.dat', 'genes.dat', 'promoters.dat', 'terminators.dat'):
    p = os.path.join(ec, name)
    if os.path.exists(p):
        shutil.copyfile(p, os.path.join(out, name)); os.chmod(os.path.join(out, name), 0o600)
print('exported:', sorted(os.listdir(out)), file=sys.stderr)
