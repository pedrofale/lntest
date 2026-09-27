"""Capsule QC: why the kidney runs drop the 35 capsules with the lowest total UMI.

Left, the distribution of capsule sizes, split at the median into the shallower and the deeper half, as
kidney/full_depth.py splits them. Right, the DEGs each method calls on that split at full depth as the
smallest capsules are removed one by one: the filter is meant to leave every method calling none.

    python -m kidney.capsule_qc      # from reproducibility/, a few minutes
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter

import plot_style
from kidney.full_depth import count_degs, load, split_by_depth
from method_colors import method_color, method_draw_order, method_label
from paths import fig_name, results_dir

N_REMOVED = 35
DEEP, SHALLOW = '0.25', '0.7'


def removal_scan(umis):
    return pd.concat([count_degs(*split_by_depth(umis, k), k) for k in range(N_REMOVED + 1)], ignore_index=True)


def plot(umis, scan, out):
    sizes = np.sort(umis.sum(1))
    median = np.median(sizes)
    deep = sizes >= median  # as full_depth.split_by_depth
    fig, (ax_hist, ax_scan) = plt.subplots(1, 2, figsize=plot_style.figsize(2))

    bins = np.logspace(np.log10(sizes.min()) - 0.05, np.log10(sizes.max()) + 0.05, 30)
    ax_hist.hist([sizes[~deep], sizes[deep]], bins=bins, stacked=True, color=[SHALLOW, DEEP], edgecolor='none')
    ax_hist.axvline(median, color='black', ls='--', lw=0.8)
    ax_hist.set_xscale('log')
    ax_hist.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:,.0f}'))
    ax_hist.set_xlabel('Total UMIs per capsule')
    ax_hist.set_ylabel('Capsules')
    top = ax_hist.get_ylim()[1] * 1.15  # headroom, so the labels clear the tallest bars
    ax_hist.set_ylim(0, top)
    ax_hist.text(median / 1.1, top * 0.98, 'median', ha='right', va='top', fontsize=plot_style.SMALL)

    bh = scan[scan['correction'] == 'BH']
    zero = plot_style.zero_row(bh['n_degs'], np.zeros(len(bh)))
    for m in method_draw_order(bh['method'].unique()):
        d = bh[bh['method'] == m].sort_values('n_removed')
        n = d['n_degs'].to_numpy()
        ax_scan.plot(d['n_removed'], np.where(n > 0, n, zero), color=method_color(m), marker='o', ms=2, lw=1.3,
                     label=method_label(m))
    plot_style.log_axis_with_zero(ax_scan, zero, bh['n_degs'].max())
    ax_scan.set_xlabel('Smallest capsules removed')
    ax_scan.set_ylabel('Number of DEGs')

    fig.tight_layout()
    plot_style.legend_outside(fig, *ax_scan.get_legend_handles_labels())
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    plot_style.use()
    umis = load()
    scan = removal_scan(umis)
    csv = results_dir(__file__) / 'capsule_removal_scan.csv'
    scan.to_csv(csv, index=False)
    out = results_dir(__file__) / fig_name('capsule_qc')
    plot(umis, scan, out)
    most = scan[scan['correction'] == 'BH'].groupby('n_removed')['n_degs'].max()
    print(f'every method calls 0 from {most[most > 0].index.max() + 1} removed (BH)')
    print(f'wrote {csv} and {out}')
