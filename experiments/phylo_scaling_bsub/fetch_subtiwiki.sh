#!/bin/bash
# SubtiWiki v5 API: operon list (gene ids) + per-gene genomic annotation (locus tag, coordinates, strand) for strain 168.
# usage: fetch_subtiwiki.sh OUTDIR   (writes operon.json, operon_v4.csv, genes/<id>.json, gene_coords.tsv)
set -euo pipefail
out=$1; mkdir -p "$out/genes"; cd "$out"
api=https://subtiwiki.uni-goettingen.de/v5/api
[ -s operon.json ] || curl -sfL -m 120 -o operon.json $api/operon/
[ -s operon_v4.csv ] || curl -sfL -m 120 -o operon_v4.csv https://subtiwiki.uni-goettingen.de/v4/operon/exporter
python3 -I -c "import json;d=json.load(open('operon.json'))['data'];print('\n'.join(map(str,sorted({g['id'] for o in d for g in o['genes']}))))" > gene_ids.txt
xargs -P 3 -I{} sh -c '[ -s genes/{}.json ] || curl -sfL -m 60 -o genes/{}.json '$api'/gene/{} || echo fail {}' < gene_ids.txt
python3 -I - <<'P'
import json, glob
with open('gene_coords.tsv', 'w') as o:
    o.write('id\tname\tlocus_tag\tstart\tend\tstrand\n')
    for f in sorted(glob.glob('genes/*.json')):
        d = json.load(open(f))['data']
        for a in d.get('genomic_annotations') or []:
            if a.get('strain_id') == 1:
                o.write(f"{a['gene']['id']}\t{a['gene']['name']}\t{a.get('locus_tag')}\t{a['start']}\t{a['end']}\t{a['orientation']}\n")
P
wc -l gene_coords.tsv
