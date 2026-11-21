"""Prune the GTDB r226 species trees (bac120, ar53) to each clade's PanGBank species.

Every GTDB tree tip is a species representative genome, so a species maps to the tip whose taxonomy names it.
Output per clade: OUTDIR/<clade>.nwk, a rooted binary Newick tree with tips named by PanGBank id; unary nodes left
by pruning are merged (branch lengths summed), polytomies are resolved with zero-length branches. Also writes
OUTDIR/tips.tsv (id, species, accession, clade, in_tree).

usage: gtdb_tree.py CLADES.tsv GTDB_DIR OUTDIR
"""
import gzip, os, re, sys
import pandas as pd


def parse_newick(s):
    """iterative Newick parser (the bac120 tree is too deep for recursion). Returns parent, bl, name lists;
    node 0 is the root. Quoted labels ('100.0:p__X') are kept as names of internal nodes."""
    parent, bl, name = [-1], [0.0], ['']
    cur = 0; i = 0; n = len(s)
    while i < n:
        ch = s[i]
        if ch == '(':
            parent.append(cur); bl.append(0.0); name.append(''); cur = len(parent) - 1; i += 1
        elif ch == ',':
            p = parent[cur]
            parent.append(p); bl.append(0.0); name.append(''); cur = len(parent) - 1; i += 1
        elif ch == ')':
            cur = parent[cur]; i += 1
        elif ch == ';':
            break
        elif ch == ':':
            j = i + 1
            while j < n and s[j] not in ',();': j += 1
            bl[cur] = float(s[i + 1:j]); i = j
        elif ch == "'":
            j = s.index("'", i + 1); name[cur] = s[i + 1:j]; i = j + 1
        elif ch.isspace():
            i += 1
        else:
            j = i
            while j < n and s[j] not in ':,();': j += 1
            name[cur] = s[i:j]; i = j
    return parent, bl, name


def prune(parent, bl, name, keep_names, rename):
    """subtree spanning the kept tips; returns Newick string"""
    N = len(parent)
    children = [[] for _ in range(N)]
    for v in range(1, N): children[parent[v]].append(v)
    keep = [False] * N
    tipset = set(keep_names)
    for v in range(N):
        if not children[v] and name[v] in tipset: keep[v] = True
    order = []                                  # post-order
    stack = [(0, False)]
    while stack:
        v, done = stack.pop()
        if done: order.append(v); continue
        stack.append((v, True))
        for c in children[v]: stack.append((c, False))
    for v in order:
        if children[v] and any(keep[c] for c in children[v]): keep[v] = True
    # build pruned structure: each kept node -> list of (kept child subtree string, length)
    sub = {}
    for v in order:
        if not keep[v]: continue
        kc = [c for c in children[v] if keep[c]]
        if not kc:
            sub[v] = (rename[name[v]], bl[v])
        elif len(kc) == 1:                       # unary: merge into the child
            s_, l_ = sub.pop(kc[0]); sub[v] = (s_, l_ + bl[v])
        else:
            parts = [sub.pop(c) for c in kc]
            while len(parts) > 2:                # resolve polytomy with zero-length branches
                a, b = parts.pop(), parts.pop()
                parts.append((f'({a[0]}:{a[1]:.6g},{b[0]}:{b[1]:.6g})', 0.0))
            sub[v] = (f'({parts[0][0]}:{parts[0][1]:.6g},{parts[1][0]}:{parts[1][1]:.6g})', bl[v])
    return sub[0][0] + ';'


clades_f, gdir, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)
C = pd.read_csv(clades_f, sep='\t')
tax = {}
for f in ('bac120_taxonomy_r226.tsv.gz', 'ar53_taxonomy_r226.tsv.gz'):
    for line in gzip.open(f'{gdir}/{f}', 'rt'):
        acc, t = line.rstrip('\n').split('\t')
        tax[acc] = t.split(';s__')[-1]
trees = {}
for dom, f in (('bac', 'bac120_r226.tree.gz'), ('ar', 'ar53_r226.tree.gz')):
    trees[dom] = parse_newick(gzip.open(f'{gdir}/{f}', 'rt').read())
sp2tip = {}
for dom, (parent, bl, name) in trees.items():
    childless = set(range(len(parent))) - set(parent)
    for v in childless:
        if name[v] in tax: sp2tip[tax[name[v]]] = (dom, name[v])
rows = []
for clade, sub in C.groupby('clade'):
    m = {r.id: sp2tip.get(r.species) for r in sub.itertuples()}
    ok = {k: v for k, v in m.items() if v}
    doms = {v[0] for v in ok.values()}
    assert len(doms) == 1, (clade, doms)
    dom = doms.pop()
    rename = {v[1]: str(k) for k, v in ok.items()}
    nwk = prune(*trees[dom], list(rename), rename)
    open(f'{out}/{clade}.nwk', 'w').write(nwk + '\n')
    for r in sub.itertuples():
        rows.append(dict(id=r.id, species=r.species, accession=m[r.id][1] if m[r.id] else '', clade=clade,
                         in_tree=bool(m[r.id])))
    print(f'{clade}: {len(ok)}/{len(sub)} species in the {dom} tree')
pd.DataFrame(rows).to_csv(f'{out}/tips.tsv', sep='\t', index=False)
