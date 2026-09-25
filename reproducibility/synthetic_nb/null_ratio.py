import argparse
import math

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

from method_colors import method_color
from paths import fig_name, require_input, results_dir

INK, MUTED = '#0b0b0b', '#52514e'
METHODS = {'log1p t-test': r'log1p $t$-test', "LN's t-test": r"LN's $t$-test"}
# Var(X) from lowest (Poisson) to highest: solid, long dashes, short dashes, dots
VAR_STYLES = ['-', (0, (6, 2)), (0, (3, 2)), (0, (1, 1.5))]
ALPHA = 0.8


def load(base_mus):
    out = results_dir(__file__)
    frames = []
    for mu in base_mus:
        stem = f'variance_vs_fpr_mu{int(mu)}_ratio_grid'
        path = require_input(out / f'{stem}.csv', what=f"ratio-grid null sweep at base_mu={mu:g}",
                             source=f'python -m synthetic_nb.null --base-mu {mu:g}')
        frames.append(pd.read_csv(path))
    df = pd.concat(frames)
    df['ratio'] = df.var_y / df.var_x
    return df


def ratio_axis(ax, sub):
    """Linear x-axis from 1, labelled at every integer, ending at the first integer
    by which every log1p curve has reached FPR 0.99 (or the last ratio, if one never does)."""
    hi = sub.ratio.max()
    lp = sub[(sub.method == 'log1p t-test') & (sub.ratio > 1)]
    reached = lp[lp.value >= 0.99].groupby('var_x').ratio.min()
    if len(reached) == lp.var_x.nunique():
        hi = min(math.ceil(reached.max()), hi)
    ax.set_xticks(range(1, int(hi) + 1))
    pad = 0.02 * (hi - 1)
    ax.set_xlim(1 - pad, hi + pad)


def plot(sub, mu):
    # method by colour, Var(X) by line style
    var_xs = sorted(sub.var_x.unique())
    fig, ax = plt.subplots(figsize=(3.6, 2.6))
    for method in METHODS:
        for vx, ls in zip(var_xs, VAR_STYLES):
            s = sub[(sub.var_x == vx) & (sub.method == method)].sort_values('ratio')
            ax.plot(s.ratio, s.value, color=method_color(method), lw=1.3, ls=ls, alpha=ALPHA)

    ax.set_title(rf'$\mu$ = {mu:g}', fontsize=8, color=INK)
    ratio_axis(ax, sub)
    ax.set_xlabel(r'Var$(Y)$ / Var$(X)$')
    ax.set_ylabel('False positive rate')
    ax.set_ylim(-0.03, 1.03)

    style = {'fontsize': 6, 'title_fontsize': 6, 'frameon': False, 'loc': 'upper left',
             'alignment': 'left'}
    methods = ax.legend(handles=[Line2D([], [], color=method_color(m), lw=1.4, alpha=ALPHA, label=label)
                                 for m, label in METHODS.items()],
                        title='Method', bbox_to_anchor=(1.01, 1.0), **style)
    ax.add_artist(methods)
    var_key = ax.legend(handles=[Line2D([], [], color=INK, lw=1.3, ls=ls, label=f'{vx:g}')
                                 for vx, ls in zip(var_xs, VAR_STYLES)],
                        title=r'Var$(X)$', bbox_to_anchor=(1.01, 0.62), handlelength=3, **style)

    out = results_dir(__file__) / fig_name(f'variance_ratio_vs_fpr_mu{int(mu)}')
    fig.savefig(out, dpi=200, bbox_inches='tight', bbox_extra_artists=[methods, var_key])
    plt.close(fig)
    print(f'wrote {out}')


def run(base_mus):
    df = load(base_mus)
    n = df.groupby('base_mu').var_x.nunique().max()
    if n > len(VAR_STYLES):
        raise SystemExit(f'{n} values of Var(X), but only {len(VAR_STYLES)} line styles')

    plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.edgecolor': MUTED, 'xtick.color': MUTED, 'ytick.color': MUTED,
                         'axes.labelcolor': INK})
    for mu in base_mus:
        plot(df[df.base_mu == mu], mu)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(
        description="null.py's FPR sweeps against Var(Y)/Var(X), one figure per base_mu."
    )
    ap.add_argument('--base-mu', type=float, nargs='+', default=[5, 50],
                    help='sweeps to plot, each from null.py --base-mu')
    run(ap.parse_args().base_mu)
