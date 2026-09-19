import argparse

import numpy as np
import pandas as pd
import statsmodels.stats.multitest as smm
from lntest import get_LN_lfcs
import matplotlib.pyplot as plt
from baselines import scanpy_sig_test, get_test_results
from paths import results_dir

NX = 1000
NY = 1000
N_GENES = 1500


def run(base_mu, metric, seed=0):
    np.random.seed(seed)
    mu1 = base_mu + np.zeros((1, N_GENES))
    mu2 = base_mu
    true_lfcs = np.log2(mu2 / mu1)

    plt.rcParams.update({'font.size': 50})
    fig, (sc_ax, LN_ax) = plt.subplots(1, 2, figsize=(30, 15), sharey=True)

    d1_list = np.linspace(0.01, 1, 4)
    d2_list = np.linspace(0.0001, 2, 20)

    var_1 = d1_list * base_mu ** 2 + base_mu
    var_2 = d2_list * base_mu ** 2 + base_mu

    rows = []
    for i, d1 in enumerate(d1_list):
        results = {"sc": [], "LN": []}
        for d2 in d2_list:
            r1 = 1 / d1
            r2 = 1 / d2
            p1 = 1 / (1 + d1 * mu1)
            p2 = 1 / (1 + d2 * mu2)

            # Generate synthetic gene expression data
            X = np.random.negative_binomial(n=r1, p=p1, size=(NX, N_GENES))
            Y = np.random.negative_binomial(n=r2, p=p2, size=(NY, N_GENES))

            sc_lfs, sc_adj_pvals = scanpy_sig_test(X, Y)
            if np.sum(sc_adj_pvals >= 0.05) == N_GENES:
                # 100% accuracy and 0% TPR --- conf_mat crashes in this setting
                sc_results = {"accuracy": 1., "fpr": 0.}
            else:
                sc_results = get_test_results(sc_adj_pvals, true_lfcs, verbose=False)

            LN_lfcs, LN_p_vals = get_LN_lfcs(Y, X, test='t')
            LN_adj_pvals = smm.multipletests(LN_p_vals, alpha=0.05, method='bonferroni')[1]
            if np.sum(LN_adj_pvals >= 0.05) == N_GENES:
                # 100% accuracy and 0% TPR --- conf_mat crashes in this setting
                ln_results = {"accuracy": 1., "fpr": 0.}
            else:
                ln_results = get_test_results(LN_adj_pvals, true_lfcs, verbose=False)

            results["sc"] += [sc_results[metric]]
            results["LN"] += [ln_results[metric]]

        for method, key in (('log1p t-test', 'sc'), ("LN's t-test", 'LN')):
            rows += [
                {'base_mu': base_mu, 'metric': metric, 'var_x': var_1[i],
                 'var_y': vy, 'method': method, 'value': v}
                for vy, v in zip(var_2, results[key])
            ]

        sc_ax.plot(var_2.astype(int), results["sc"], lw=5)
        LN_ax.plot(var_2.astype(int), results["LN"],
                   label=f"Var$(X$)={int(var_1[i])}", lw=5)

    sc_ax.set_xlabel(r'Var$(Y)$')
    LN_ax.set_xlabel(r'Var$(Y)$')
    LN_ax.set_title(r"LN's $t$-test")
    sc_ax.set_title(r'$\log1$p $t$-test')
    sc_ax.set_ylabel(metric.upper())
    if metric == 'fpr':
        sc_ax.set_yticks(np.arange(0, 11) * 1 / 10)
    LN_ax.legend(loc='best')

    out = results_dir(__file__)
    stem = f'variance_vs_{metric}_mu{int(base_mu)}'
    fig.savefig(out / f'{stem}.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    pd.DataFrame(rows).to_csv(out / f'{stem}.csv', index=False)
    print(f'wrote {out / stem}.{{png,csv}}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(
        description="Null-scenario FPR/accuracy against Var(Y), log1p t-test vs LN's t-test."
    )
    ap.add_argument('--base-mu', type=float, default=50,
                    help="shared NB mean (was hardcoded). RECOMB's var_v_fpr.png "
                         "came from 5; var_v_fpr_large_mu.png from this default")
    ap.add_argument('--metric', default='fpr', choices=['fpr', 'accuracy'])
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()
    run(args.base_mu, args.metric, args.seed)
