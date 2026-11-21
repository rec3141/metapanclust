"""metapanclust command line.

  python -m metapanclust ref   PANGENOME.h5 REFDIR
  python -m metapanclust score REFDIR BIN.fna OUTPREFIX [--comp 0.93] [--threads 4] [--no-context]
  python -m metapanclust eval  REFDIR FOREIGN_REFDIR OUTDIR [--n 20] [--threads 4]

score writes OUTPREFIX.summary.tsv (one row), OUTPREFIX.units.tsv (every reference unit: p, detection,
observed, bits, calibrated excess) and OUTPREFIX.contigs.tsv (per contig: novel genes, persistent genes,
presence-direction excess, median identity to the reference).
"""
import argparse, json, os, sys


def main():
    ap = argparse.ArgumentParser(prog='metapanclust')
    sub = ap.add_subparsers(dest='cmd', required=True)
    a = sub.add_parser('ref'); a.add_argument('h5'); a.add_argument('refdir')
    a = sub.add_parser('score'); a.add_argument('refdir'); a.add_argument('fasta'); a.add_argument('out')
    a.add_argument('--comp', type=float); a.add_argument('--threads', type=int, default=4)
    a.add_argument('--no-context', action='store_true')
    a = sub.add_parser('eval'); a.add_argument('refdir'); a.add_argument('frefdir'); a.add_argument('outdir')
    a.add_argument('--n', type=int, default=20); a.add_argument('--threads', type=int, default=4)
    a.add_argument('--no-context', action='store_true'); a.add_argument('--seed', type=int, default=1)
    x = ap.parse_args()
    if x.cmd == 'ref':
        from .ref import build
        print(json.dumps(build(x.h5, x.refdir)))
    elif x.cmd == 'score':
        import pandas as pd
        from .ref import Ref
        from .annotate import annotate
        from .score import score
        genes, lengths = annotate(x.fasta, x.refdir, x.out + '.work', x.threads)
        s, U, C = score(Ref(x.refdir), genes, lengths.values(), comp=x.comp, use_context=not x.no_context)
        s['bin'] = os.path.basename(x.fasta)
        pd.DataFrame([s]).to_csv(x.out + '.summary.tsv', sep='\t', index=False)
        U.to_csv(x.out + '.units.tsv', sep='\t', index=False)
        C.to_csv(x.out + '.contigs.tsv', sep='\t', index=False)
        print(json.dumps(s))
    elif x.cmd == 'eval':
        from .simulate import run
        os.makedirs(x.outdir, exist_ok=True)
        run(x.refdir, x.frefdir, x.outdir, n=x.n, threads=x.threads, use_context=not x.no_context, seed=x.seed)


if __name__ == '__main__':
    main()
