"""Summarise context_model.py outputs for both species: log-loss (bits) per model, per-module gains with
co-occurrence partners, overfitting gap, learning curve, and Chesterfield-bin surprise.
usage: analyse.py OUTROOT"""
import sys, os, numpy as np, pandas as pd
root = sys.argv[1]
SPECIES = {'uba9145': 'chesterfield_maxbin_040', 'mag120531': 'chesterfield_semibin_019'}
BAKTA = '/data/scratch/metapanclust_data/chesterfield/mag_chesterfield/annotation/bakta/basic/annotation.tsv'
EPS = 1e-6

def bits(p, y):
    p = np.clip(p, EPS, 1 - EPS)
    return -np.where(y > 0, np.log2(p), np.log2(1 - p))

def models(df):
    return [c for c in df.columns if c.startswith(('prev', 'knn', 'lr'))]

def annotate(sp, bin_):
    d = f'/data/scratch/metapanclust_data/pangenomes/{sp}/ppanggolin'
    gff = pd.read_csv(f'{d}/gff/{bin_}.gff', sep='\t', comment='#', header=None)
    gff = gff[gff[2] == 'CDS']
    gff['fam'] = gff[8].str.extract(r'family=([^;]+)')[0]
    ann = pd.read_csv(BAKTA, sep='\t', comment='#', header=None, usecols=[0, 2, 3, 7], names=['c', 's', 'e', 'prod'])
    ann = ann.set_index(['c', 's', 'e'])['prod']
    prod = {}
    for c, s, e, f in zip(gff[0], gff[3], gff[4], gff.fam):
        p = ann.get((c, s, e))
        if isinstance(p, str) and p != 'hypothetical protein': prod.setdefault(f, p)
    fm = pd.read_csv(f'{d}/functional_modules.tsv', sep='\t')
    out = {}
    for m, g in fm.groupby('module_id'):
        ps = [prod[f] for f in g.family_id if f in prod]
        out[m] = '; '.join(sorted(set(ps), key=ps.index)[:3]) if ps else ''
    return out

def module_matrix(sp):
    d = f'/data/scratch/metapanclust_data/pangenomes/{sp}/ppanggolin'
    mig = pd.read_csv(f'{d}/modules_in_genomes.tsv', sep='\t'); mig = mig[mig.completion >= 0.5]
    M = pd.crosstab(mig.module_id, mig.genome) > 0
    spots = pd.read_csv(f'{d}/modules_spots.tsv', sep='\t').groupby('module_id').spot_id.apply(lambda s: ','.join(sorted(s))).to_dict()
    return M, spots

summary = {}
for sp, bin_ in SPECIES.items():
    o = f'{root}/{sp}'
    if not os.path.exists(f'{o}/chesterfield_predictions.tsv'): print(f'{sp}: not finished'); continue
    loo = pd.read_csv(f'{o}/loo_predictions.tsv', sep='\t'); ins = pd.read_csv(f'{o}/insample_predictions.tsv', sep='\t')
    ms = models(loo)
    for m in ms: loo['b_' + m] = bits(loo[m], loo.y); ins['b_' + m] = bits(ins[m], ins.y)
    print(f'\n######## {sp}: LOO log-loss, mean bits per (genome, target); in-sample in brackets')
    tab = []
    for kind in ('module', 'family'):
        L = loo[loo.kind == kind]; I = ins[ins.kind == kind]
        if L.empty: continue
        per_t = L.groupby('target')[['b_' + m for m in ms]].mean()
        base = per_t['b_prev']
        for m in ms:
            g = base - per_t['b_' + m]
            tab.append(dict(kind=kind, model=m, loo=L['b_' + m].mean(), insample=I['b_' + m].mean(),
                            gain=g.mean(), frac_targets_better=(g > 0.05).mean(), frac_worse=(g < -0.05).mean(),
                            n_targets=len(per_t)))
    tab = pd.DataFrame(tab); summary[sp] = tab
    print(tab.to_string(index=False, float_format='%.3f'))

    # per-module detail for the fixed model choices
    M, spots = module_matrix(sp)
    ann = annotate(sp, bin_)
    pubg = [g for g in M.columns if not g.startswith('chesterfield')]
    Mp = M.reindex(columns=sorted(set(pubg) | set(loo.genome.unique())), fill_value=False).astype(float)
    L = loo[loo.kind == 'module']
    per_t = L.groupby('target')[['b_' + m for m in ms]].mean()
    corr = np.corrcoef(Mp.values)
    idx = {m: i for i, m in enumerate(Mp.index)}
    rows = []
    for t, r in per_t.iterrows():
        i = idx[t]; c = corr[i].copy(); c[i] = 0; c = np.nan_to_num(c)
        j = int(np.argmax(np.abs(c))); partner = Mp.index[j]
        rows.append(dict(module=t, n_pub=int(Mp.loc[t].sum()), prev=r['b_prev'], knn8=r['b_knn8'], lrmod1=r['b_lrmod1'],
                         lrfam=r['b_lrfam0.1'], gain_lrmod=r['b_prev'] - r['b_lrmod1'], gain_knn=r['b_prev'] - r['b_knn8'],
                         partner=partner, phi=c[j], same_spot=bool(spots.get(t) and spots.get(t) == spots.get(partner)),
                         spot=spots.get(t, ''), ann=ann.get(t, '')[:70]))
    det = pd.DataFrame(rows).sort_values('gain_lrmod', ascending=False)
    det.to_csv(f'{o}/module_gains.tsv', sep='\t', index=False)
    print(f'\n{sp}: modules by gain of lrmod1 over prevalence (bits/genome), with most correlated module')
    print(det.head(15).to_string(index=False, float_format='%.2f'))
    print('worst:'); print(det.tail(5).to_string(index=False, float_format='%.2f'))
    print(f'modules with |phi| >= 0.8 partner: {(det.phi.abs() >= 0.8).sum()}/{len(det)}, negative phi <= -0.5: {(det.phi <= -0.5).sum()}')

    # learning curve
    lc = pd.read_csv(f'{o}/learning_curve.tsv', sep='\t') if os.path.exists(f'{o}/learning_curve.tsv') else None
    if lc is None: lms = []
    else: lms = models(lc)
    for m in lms: lc['b_' + m] = bits(lc[m], lc.y)
    g = None if lc is None else lc.groupby(['kind', 'n'])[['b_' + m for m in lms]].mean()
    for m in (lms if g is not None else []):
        if m != 'prev': g['gain_' + m] = g['b_prev'] - g['b_' + m]
    if g is not None: print(f'\n{sp}: learning curve, gain over prevalence (bits per genome-target) vs training genomes')
    if g is not None: print(g[[c for c in g.columns if c.startswith('gain_')] + ['b_prev']].to_string(float_format='%.3f'))
    if g is not None: g.to_csv(f'{o}/learning_curve_summary.tsv', sep='\t')

    # Chesterfield bins
    ch = pd.read_csv(f'{o}/chesterfield_predictions.tsv', sep='\t')
    for m in ms: ch['b_' + m] = bits(ch[m], ch.y)
    # README convention for the prevalence model: absent discounted by completeness
    ch['conv_prev'] = np.where(ch.y > 0, -np.log2(ch.prev), -np.log2(1 - ch.prev * ch.c))
    ch['ann'] = ch.module.map(ann).fillna('').str[:60]
    ch.to_csv(f'{o}/chesterfield_surprise.tsv', sep='\t', index=False)
    print(f'\n{sp}: Chesterfield bins, total surprise over all modules (bits)')
    cols = ['conv_prev', 'b_prev', 'b_prevc', 'b_knn8', 'b_lrmod1', 'b_lrfam0.1']
    print(ch.groupby(['setting', 'genome'])[cols].sum().to_string(float_format='%.1f'))
    novel = ch[(ch.n_pub == 0) & (ch.y > 0)]
    print(f'Chesterfield-only modules present in the bin ({novel.module.nunique()}): mean bits per module')
    print(novel.groupby(['setting', 'genome'])[cols].mean().to_string(float_format='%.2f'))
    rest = ch[~((ch.n_pub == 0) & (ch.y > 0))]
    print('all other module states: total bits')
    print(rest.groupby(['setting', 'genome'])[cols].sum().to_string(float_format='%.1f'))
    top = ch[ch.setting == 'public'].assign(d=lambda x: x.b_prev - x.b_lrmod1)
    top = top[top.n_pub > 0].sort_values('b_prev', ascending=False)
    print('most surprising states with public carriers (public training):')
    print(top.head(12)[['genome', 'module', 'y', 'n_pub', 'n_fam', 'b_prev', 'b_knn8', 'b_lrmod1', 'b_lrfam0.1', 'ann']].to_string(index=False, float_format='%.2f'))

# cross-species hyperparameter choice: pick on one species, evaluate on the other
if len(summary) == 2:
    print('\n######## hyperparameters chosen on one species, LOO gain (bits/genome-target) on the other')
    for a, b in (('uba9145', 'mag120531'), ('mag120531', 'uba9145')):
        for kind in ('module', 'family'):
            A = summary[a][summary[a].kind == kind].set_index('model'); B = summary[b][summary[b].kind == kind].set_index('model')
            for fam in ('knn', 'lrmod', 'lrfam'):
                best = A[A.index.str.startswith(fam)].gain.idxmax()
                print(f'  {kind:6s} {fam:5s}: chosen on {a} = {best:10s} -> gain on {b} {B.loc[best, "gain"]:+.3f} (best possible {B[B.index.str.startswith(fam)].gain.max():+.3f})')
