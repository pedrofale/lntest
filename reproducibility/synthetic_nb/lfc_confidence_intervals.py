import numpy as np
from lntest import get_LN_lfcs as get_DELN_lfcs
from baselines import scanpy_sig_test
from paths import fig_name, results_dir
import plot_style
from method_colors import method_color, method_label


def _save(name):
    out = results_dir(__file__) / fig_name(name)
    plot_style.rasterize_dense(plt.gcf())
    plt.savefig(out, dpi=plot_style.RASTER_DPI, bbox_inches='tight')
    plt.close()
    print(f'wrote {out}')
import matplotlib.pyplot as plt

"""
Scanpy function to get mean and vars used for computing the t statistic

def _get_mean_var(
    X: _SupportedArray, *, axis: Literal[0, 1] = 0
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    if isinstance(X, sparse.spmatrix):
        mean, var = sparse_mean_variance_axis(X, axis=axis)
    else:
        mean = axis_mean(X, axis=axis, dtype=np.float64)
        mean_sq = axis_mean(elem_mul(X, X), axis=axis, dtype=np.float64)
        var = mean_sq - mean**2
    # enforce R convention (unbiased estimator) for variance
    var *= X.shape[axis] / (X.shape[axis] - 1)
    return mean, var
"""

np.random.seed(10)
nx = 1000
ny = 1000
n_genes = 1500
mu1 = 10 + np.abs(np.random.normal(0, 5, (1, n_genes))) * np.random.binomial(1, 0.1, (1, n_genes))
mu2 = 10 * np.exp(1.)
d1 = 2.
d2 = 1.

r1 = 1 / d1
r2 = 1 / d2
p1 = 1 / (1 + d1 * mu1)
p2 = 1 / (1 + d2 * mu2)

# Generate synthetic gene expression data
X = np.random.negative_binomial(n=r1, p=p1, size=(nx, n_genes))
Y = np.random.negative_binomial(n=r2, p=p2, size=(ny, n_genes))

true_lfcs = np.log2(mu2 * np.exp(-1.) / mu1)
non_de_idx = (true_lfcs == 0.).reshape(-1)


sc_lfcs = scanpy_sig_test(X, Y)[0].to_numpy()  # the Series index is in ranked order
X_raw, Y_raw = X.copy(), Y.copy()
# this is how Scanpy calculates the statistics used in the t test (based on their code)
# the results here were reproduced in the debugging console within the scanpy code, too
X = np.log(1e4 * X / (X.sum(1, keepdims=True)) + 1)
Y = np.log(1e4 * Y / (Y.sum(1, keepdims=True)) + 1)

mean_Y, mean_X = np.mean(Y, axis=0), np.mean(X, axis=0)

mean_sq_Y, mean_sq_X = np.mean(Y * Y, axis=0), np.mean(X * X, axis=0)
var_Y = (mean_sq_Y - mean_Y ** 2) * ny / (ny - 1)
var_X = (mean_sq_X - mean_X ** 2) * nx / (nx - 1)

se_Y = np.sqrt(var_Y / ny)
se_X = np.sqrt(var_X / nx)
se = np.sqrt(se_Y ** 2 + se_X ** 2) / np.log(2)
# for simplicity, the confidence intervals are formed with a critical z value (not t)
# but for these sample sizes the values are equal
confidence_intervals = (mean_Y - mean_X) / np.log(2) + 1.96 * np.array([-se, se])
print(f"Estimated LFC is below CI {np.sum(sc_lfcs[~non_de_idx] < confidence_intervals[0, ~non_de_idx])} times for DEGs")
print(f"Estimated LFC is above CI {np.sum(sc_lfcs[~non_de_idx] > confidence_intervals[1, ~non_de_idx])} times for DEGs")
print()
print(f"Estimated LFC is below CI {np.sum(sc_lfcs[non_de_idx] < confidence_intervals[0, non_de_idx])} times for non-DEGs")
print(f"Estimated LFC is above CI {np.sum(sc_lfcs[non_de_idx] > confidence_intervals[1, non_de_idx])} times for non-DEGs")
print(f"Frequency of estimated LFC that are above CI {np.sum(sc_lfcs[non_de_idx] > confidence_intervals[1, non_de_idx]) / np.sum(non_de_idx)} for  non-DEGs")
print()
plot_style.use()
plt.figure(figsize=plot_style.figsize())
plt.errorbar(true_lfcs[0], (mean_Y - mean_X) / np.log(2), yerr=1.96 * se, fmt='none', label='Estimated CIs', color=method_color('t-test'), alpha=0.5, lw=0.8)
plt.plot(true_lfcs[0], sc_lfcs, 'o', ms=2, label='Estimated LFCs', color=method_color('t-test'))
plt.plot(true_lfcs[0], true_lfcs[0], 'o', ms=2, label='True LFCs', color='0.7')
plt.xlabel('LFCs')
plt.ylabel('LFCs')
plt.title(method_label('t-test'))
plt.legend()
plt.tight_layout()
_save('lfc_ci_log1p.png')

deln_lfcs, _, gamma = get_DELN_lfcs(Y_raw, X_raw, return_standard_error=True)
deln_confidence_intervals = deln_lfcs + 1.96 * np.array([-gamma, gamma])
print(f"LN's Estimated LFC is below CI {np.sum(deln_lfcs < deln_confidence_intervals[0])} times")
print(f"LN's Estimated LFC is above CI {np.sum(deln_lfcs > deln_confidence_intervals[1])} times")

plt.figure(figsize=plot_style.figsize())
plt.errorbar(true_lfcs[0], deln_lfcs, yerr=1.96 * gamma, fmt='none', label='Estimated CIs', color=method_color('LN'), alpha=0.5, lw=0.8)
plt.plot(true_lfcs[0], deln_lfcs, 'o', ms=2, label='Estimated LFCs', color=method_color('LN'))
plt.plot(true_lfcs[0], true_lfcs[0], 'o', ms=2, label='True LFCs', color='0.7')
plt.xlabel('LFCs')
plt.ylabel('LFCs')
plt.title(method_label('LN'))
plt.legend()
plt.tight_layout()
_save('lfc_ci_ln.png')
plt.show()