"""Merge-rule AUC and rejection rate against contaminant ANI (one point per focal x contaminant species).
usage: 06_plot.py RESULTS_DIR FIG_DIR
"""
import os, sys
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

R, F = sys.argv[1:3]
os.makedirs(F, exist_ok=True)
P = pd.read_csv(f'{R}/pairs.tsv', sep='\t')
NOANI = 76  # x position for contaminants below skani's range (no ANI reported)
P['x'] = P.ani.fillna(NOANI)
INK, MUTED, GRID, SURF = '#2b2b2a', '#6b6b68', '#e4e4e1', '#fcfcfb'
# categorical slots in fixed order; shape is the secondary encoding (4 series on a scatter)
STYLE = {'Alteromonas': ('#2a78d6', 'o'), 'Pseudoalteromonas': ('#eb6834', 's'),
         'Shewanella': ('#1baf7a', '^'), 'Sulfitobacter': ('#eda100', 'D')}
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'xtick.color': MUTED, 'ytick.color': MUTED, 'text.color': INK})
fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), facecolor=SURF)
panels = [('auc_excess', 'Merge rule AUC\n(contaminant vs same-species merge)', 'Full surprise (with absences)'),
          ('auc_excess_content', 'Merge rule AUC', 'Content surprise only (no absence term)'),
          ('rej_contam', 'Contaminant merges rejected', 'Rejection rate, full surprise')]
for ax, (col, ylab, title) in zip(axes, panels):
    ax.set_facecolor(SURF)
    for g, d in P.groupby('genus'):
        c, m = STYLE.get(g, (MUTED, 'o'))
        ax.scatter(d.x, d[col], s=36, c=c, marker=m, edgecolors=SURF, linewidths=1.5, label=g, zorder=3)
    if col == 'rej_contam':
        ax.axhline(P.rej_strain.mean(), color=MUTED, lw=1, ls='--', zorder=2)
        ax.text(77, P.rej_strain.mean() + 0.03, 'same-species merges rejected (mean)', color=MUTED, fontsize=8)
    else:
        ax.axhline(0.5, color=MUTED, lw=1, ls='--', zorder=2)
    ax.axvline(95, color=GRID, lw=1, zorder=1)
    ax.set_ylim(-0.03, 1.05); ax.set_xlim(74, 100)
    ax.set_xticks([NOANI, 80, 85, 90, 95, 100]); ax.set_xticklabels(['none\n(<~80)', '80', '85', '90', '95', '100'])
    ax.set_xlabel('skani ANI, contaminant vs focal held-out genomes (%)')
    ax.set_ylabel(ylab); ax.set_title(title, fontsize=9, loc='left')
    ax.grid(axis='y', color=GRID, lw=0.6); ax.set_axisbelow(True)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
axes[0].legend(frameon=False, fontsize=8, loc='lower left', title='Genus', title_fontsize=8)
fig.tight_layout()
fig.savefig(f'{F}/merge_vs_ani.png', dpi=150, facecolor=SURF)
print(f'{F}/merge_vs_ani.png')
