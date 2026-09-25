import argparse

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

from method_colors import LN as LN_COLOR, LOG1P_RAMP
from paths import fig_name, require_input, results_dir

INK, MUTED = '#0b0b0b', '#52514e'
STYLES = [('-', 'o', None), ((0, (3, 1.5)), 's', 'white'), ((0, (1, 1)), '^', 'white')]


def load(base_mus):
    out = results_dir(__file__)
    frames = []
    for mu in base_mus:
        stem = f'variance_vs_fpr_mu{int(mu)}'
        path = require_input(out / f'{stem}.csv', what=f"null sweep at base_mu={mu:g}",
                             source=f'python -m synthetic_nb.null --base-mu {mu:g}')
        frames.append(pd.read_csv(path))
    df = pd.concat(frames)
    df['ratio'] = df.var_y / df.var_x
    df['phi_x'] = ((df.var_x - df.base_mu) / df.base_mu ** 2).round(2)
    return df


def run(base_mus):
    df = load(base_mus)
    phis = sorted(df.phi_x.unique())
    if len(phis) > len(LOG1P_RAMP):
        raise SystemExit(f'{len(phis)} dispersions of X, but only {len(LOG1P_RAMP)} colours')
    color = dict(zip(phis, LOG1P_RAMP))

    plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.edgecolor': MUTED, 'xtick.color': MUTED, 'ytick.color': MUTED,
                         'axes.labelcolor': INK})
    fig, ax = plt.subplots(figsize=(3.6, 2.6))

    for mu, (ls, marker, mfc) in zip(base_mus, STYLES):
        sub = df[df.base_mu == mu]
        for phi in phis:
            s = sub[(sub.phi_x == phi) & (sub.method == 'log1p t-test')].sort_values('ratio')
            ax.plot(s.ratio, s.value, color=color[phi], lw=1.2, ls=ls,
                    marker=marker, ms=2.2, mfc=mfc or color[phi], mew=0.7)
            s = sub[(sub.phi_x == phi) & (sub.method == "LN's t-test")].sort_values('ratio')
            ax.plot(s.ratio, s.value, color=LN_COLOR, lw=1.4, ls=ls)

    ax.axvline(1, color=MUTED, lw=0.6, ls=':')
    ax.text(1, 1.06, 'equal variance', ha='center', va='bottom', color=MUTED, fontsize=7,
            transform=ax.get_xaxis_transform())
    ax.set_xscale('log', base=2)
    ax.set_xticks([1 / 32, 1 / 8, 1 / 2, 1, 2, 8, 32])
    ax.set_xticklabels(['1/32', '1/8', '1/2', '1', '2', '8', '32'])
    ax.set_xlabel(r'Var$(Y)$ / Var$(X)$')
    ax.set_ylabel('False positive rate')
    ax.set_ylim(-0.03, 1.03)

    handles = ([Line2D([], [], color=LN_COLOR, lw=1.4, label="LN's t-test (all settings)")]
               + [Line2D([], [], color=color[p], lw=1.4, label=f'log1p t-test, $\\phi_X$={p:.2f}')
                  for p in phis]
               + [Line2D([], [], color=MUTED, lw=1.4, ls=ls, label=f'$\\mu$={mu:g}')
                  for mu, (ls, _, _) in zip(base_mus, STYLES)])
    ax.legend(handles=handles, fontsize=6, frameon=False, loc='center left',
              bbox_to_anchor=(1.01, 0.5))

    out = results_dir(__file__) / fig_name('variance_ratio_vs_fpr')
    fig.savefig(out, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'wrote {out}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(
        description="null.py's FPR sweeps on one panel, against Var(Y)/Var(X)."
    )
    ap.add_argument('--base-mu', type=float, nargs='+', default=[50, 5],
                    help='sweeps to overlay, each from null.py --base-mu; '
                         'the first is drawn solid')
    args = ap.parse_args()
    if len(args.base_mu) > len(STYLES):
        ap.error(f'at most {len(STYLES)} values of --base-mu')
    run(args.base_mu)
