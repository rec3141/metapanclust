"""Strain similarity vs coverage needed: method limits from reports/Strain resolution in metagenomics.md.
Each point is a method's corner: the highest ANI between strains it separates, at its minimum coverage per strain.
The Pareto frontier is the staircase envelope: for each ANI, the lowest coverage any method in the group needs."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INK, INK2, MUTED, SURF = '#0b0b0b', '#52514e', '#8a8984', '#fcfcfb'
C = {'single': '#2a78d6', 'multi': '#eb6834', 'ref': '#8a8984'}
NR = 60.0   # y position of the "coverage not reported" band
pts = [  # name, ANI limit (%), coverage per strain (x), group, coverage reported, note
    ('Short-read assembly + binning', 95.0, NR, 'single', False, 'CAMI I/II'),
    ('metaFlye / metaMDBG (default)', 98.65, NR, 'single', False, 'collapse E. coli K-12/W'),
    ('myloasm (R10, author claim)', 99.0, NR, 'single', False, 'preprint'),
    ('Strainberry (long-read phasing)', 99.5, 20, 'single', True, 'plateau at 30x'),
    ('MetaStrainer (1 sample, <=3 strains)', 99.5, 8, 'single', True, ''),
    ('DESMAN (~30 samples)', 99.62, 5, 'multi', True, 'core genes + accessory'),
    ('STRONG (3-10+ samples)', 99.95, 5.7, 'multi', True, 'core genes only'),
    ('StrainGE (reference)', 99.89, 0.1, 'ref', True, ''),
    ('StrainScan (reference)', 99.89, 1.0, 'ref', True, ''),
    ('StrainEst (reference)', 99.89, 5.0, 'ref', True, ''),
    ('inStrain (same vs different strain only)', 99.999, 5.0, 'ref', True, ''),
]
d = lambda ani: 100.0 - ani                      # divergence (%), plotted on a reversed log axis

fig, ax = plt.subplots(figsize=(10, 6.2), dpi=160)
fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(8, 0.0006); ax.set_ylim(0.05, 110)
ax.axhspan(35, 110, color='#efeeea', zorder=0)
ax.text(0.0007, 85, 'coverage requirement not reported', color=INK2, fontsize=8.5, va='center', ha='right')


def stair(group_pts, color, ls, label):
    """lowest coverage available at each divergence: methods resolve all pairs less similar than their limit"""
    P = sorted([(d(a), y) for _, a, y, g, ok, _ in group_pts if ok], key=lambda t: -t[0])
    xs = np.logspace(np.log10(8), np.log10(0.0006), 600)
    ys = np.array([min([y for (dl, y) in P if dl <= x] or [np.nan]) for x in xs])
    ax.plot(xs, ys, color=color, lw=2, ls=ls, zorder=2, label=label, drawstyle='steps-post')


stair([p for p in pts if p[3] == 'single'], C['single'], '-', 'frontier: one sample, de novo')
stair([p for p in pts if p[3] in ('single', 'multi')], C['multi'], '--', 'frontier: many samples, de novo')
stair([p for p in pts if p[3] == 'ref'], C['ref'], ':', 'frontier: needs close reference genomes')
off = {'Short-read assembly + binning': (-6, 14), 'metaFlye / metaMDBG (default)': (0, -16),
       'myloasm (R10, author claim)': (10, 0), 'Strainberry (long-read phasing)': (10, 0),
       'MetaStrainer (1 sample, <=3 strains)': (-10, 12), 'DESMAN (~30 samples)': (-10, -14),
       'STRONG (3-10+ samples)': (10, 12), 'StrainGE (reference)': (10, 0), 'StrainScan (reference)': (10, 0),
       'StrainEst (reference)': (10, -14), 'inStrain (same vs different strain only)': (-10, -14)}
for name, ani, y, g, ok, note in pts:
    x = d(ani)
    ax.scatter([x], [y], s=70, zorder=4, facecolor=C[g] if ok else SURF, edgecolor=C[g], linewidth=2)
    ax.scatter([x], [y], s=130, zorder=3, color=SURF)          # surface ring
    dx, dy = off[name]
    if name.startswith('Short-read'): dx = 0
    ax.annotate(name, (x, y), xytext=(dx, dy), textcoords='offset points', fontsize=8.5, color=INK,
                ha='left' if dx > 0 else ('center' if dx == 0 else 'right'), va='center', zorder=5)
ticks = [95, 98, 99, 99.5, 99.9, 99.95, 99.99, 99.999]
ax.set_xticks([d(t) for t in ticks]); ax.set_xticklabels([f'{t:g}' for t in ticks])
ax.minorticks_off()
ax.set_yticks([0.1, 1, 5, 10, 20, 30]); ax.set_yticklabels(['0.1', '1', '5', '10', '20', '30'])
ax.set_xlabel('ANI between strains (%)  -  more similar to the right', color=INK2)
ax.set_ylabel('minimum coverage per strain (x)', color=INK2)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
for s in ('left', 'bottom'): ax.spines[s].set_color('#c9c8c2')
ax.tick_params(colors=INK2, labelsize=9)
ax.grid(True, which='major', color='#e6e5e0', lw=0.8, zorder=0)
fig.text(0.012, 0.975, 'Strain resolution: similarity limit vs coverage needed', color=INK, fontsize=13, va='top')
fig.text(0.012, 0.935, 'Each point: the most similar strains a method separates, at its minimum coverage per strain. '
        'Feasible region is up-left of each frontier.\nHollow: coverage not reported. Sources and caveats in '
        'reports/Strain resolution in metagenomics.md (Strainy omitted: no verified limits).',
        fontsize=8, color=INK2, va='top')
leg = ax.legend(loc='lower left', fontsize=8.5, frameon=False, labelcolor=INK2)
fig.tight_layout(rect=(0, 0, 1, 0.87))
fig.savefig('/data/dev/metapanclust/docs/figures/strain_frontier.png', facecolor=SURF)
fig.savefig('/data/dev/metapanclust/docs/figures/strain_frontier.svg', facecolor=SURF)
