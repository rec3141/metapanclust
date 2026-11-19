"""Gene table + MMseqs2 gene families for the bins Spacedust was run on.
Gene IDs follow Spacedust: <contig>_<index>, index = 0-based CDS order in the bin GFF3 (counted across contigs).
usage: 01_genes_families.py BINS_DIR BAKTA_FAA OUTDIR [min_id=0.3] [cov=0.8] [threads=2]
-> OUTDIR/genes.tsv (gene, bin, contig, idx, start, end, strand, locus_tag, product, family)"""
import glob, os, subprocess, sys
bins_dir, faa, out = sys.argv[1:4]
min_id = sys.argv[4] if len(sys.argv) > 4 else '0.3'
cov = sys.argv[5] if len(sys.argv) > 5 else '0.8'
threads = sys.argv[6] if len(sys.argv) > 6 else '2'
os.makedirs(out, exist_ok=True)
rows, tag2gene = [], {}
for g in sorted(glob.glob(f'{bins_dir}/*.gff3')):
    b = os.path.basename(g)[:-5]
    n = 0
    for l in open(g):
        if l.startswith('#'): continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'CDS': continue
        a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
        i = n; n += 1
        gid = f'{f[0]}_{i}'
        tag = a.get('locus_tag', a.get('ID'))
        rows.append([gid, b, f[0], str(i), f[3], f[4], f[6], tag, a.get('product', '').replace('%2C', ',')])
        tag2gene[tag] = gid
sel = f'{out}/bin_proteins.faa'
with open(sel, 'w') as o:
    keep = False
    for l in open(faa):
        if l.startswith('>'):
            t = l[1:].split()[0]; keep = t in tag2gene
            if keep: o.write(f'>{tag2gene[t]}\n')
        elif keep: o.write(l)
tmp = f'{out}/mmseqs_tmp'
if not os.path.exists(f'{out}/fam_cluster.tsv'):              # reuse an existing clustering
  subprocess.run(['mmseqs', 'easy-cluster', sel, f'{out}/fam', tmp, '--min-seq-id', min_id, '-c', cov,
                '--cov-mode', '0', '--threads', threads, '-s', '7.5', '--cluster-mode', '0'], check=True,
               stdout=subprocess.DEVNULL)
  subprocess.run(['rm', '-rf', tmp, f'{out}/fam_all_seqs.fasta'])
fam = {}
for l in open(f'{out}/fam_cluster.tsv'):
    r, m = l.split(); fam[m] = r
ids = {r: f'GF{k:06d}' for k, r in enumerate(sorted(set(fam.values())))}
with open(f'{out}/genes.tsv', 'w') as o:
    o.write('gene\tbin\tcontig\tidx\tstart\tend\tstrand\tlocus_tag\tproduct\tfamily\n')
    for r in rows:
        o.write('\t'.join(r + [ids.get(fam.get(r[0]), 'NA')]) + '\n')
print(f'{len(rows)} genes, {len(ids)} families ({sum(1 for r in rows if r[0] not in fam)} without protein)', file=sys.stderr)
