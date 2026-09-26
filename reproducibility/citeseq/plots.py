"""Every CITE-seq figure, from the CSVs that exp.py and the R chain write; nothing is drawn in R.

    python -m citeseq.plots      # from reproducibility/
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

import plot_style
from method_colors import method_color, method_label
from paths import fig_name, require_input, results_dir

RESULTS = results_dir(__file__)
GATE_SHADES = {'-': '0.8', '+': '0.5', '++': '0.2'}
# Scanpy's and Seurat's log1p t-tests appear side by side here, so both are named by tool
SCANPY_LOG1P = f"Scanpy {method_label('t-test')}"


def save(fig, stem, verbose=True):
    out = RESULTS / fig_name(stem)
    plot_style.rasterize_dense(fig)
    fig.savefig(out, dpi=plot_style.RASTER_DPI, bbox_inches='tight')
    plt.close(fig)
    if verbose:
        print(f'wrote {out}')


def gating():
    """CLR density of each gating marker, one fill per GMM gate."""
    df = pd.read_csv(require_input(RESULTS / 'adt_gating.csv',
                                   what='CLR values and GMM gates of the PBMC10k ADTs',
                                   source='cd citeseq/R && Rscript pbmc10k_process.R'))
    cd4_t = df[(df.CD3_cluster == '+') & (df.CD4_cluster == '++')]  # CD45RA is gated within these
    for marker, sub in (('CD3', df), ('CD4', df), ('CD45RA', cd4_t)):
        fig, ax = plt.subplots(figsize=plot_style.figsize())
        grid = np.linspace(sub[marker].min(), sub[marker].max(), 400)
        for gate in [g for g in GATE_SHADES if g in set(sub[f'{marker}_cluster'])]:
            x = sub.loc[sub[f'{marker}_cluster'] == gate, marker]
            ax.fill_between(grid, gaussian_kde(x, bw_method='silverman')(grid), color=GATE_SHADES[gate],
                            alpha=0.7, lw=0, label=f"{marker}{gate.replace('-', '−')}")
        ax.set_xlabel(f'{marker} (CLR)')
        ax.set_ylabel('Density')
        ax.set_title(f'{marker} CLR distribution with GMM clustering')
        fig.tight_layout()
        plot_style.legend_outside(fig)
        save(fig, f'{marker}_density_plot')


def annotate_r2(ax, true_lfc, est):
    """Squared Pearson correlation of estimated and true LFC over the plotted genes, top left."""
    true_lfc, est = np.asarray(true_lfc), np.asarray(est)
    ok = np.isfinite(true_lfc) & np.isfinite(est)
    r = np.corrcoef(true_lfc[ok], est[ok])[0, 1]
    ax.text(0.14, 0.84, f'$R^2$ = {r ** 2:.3f}', transform=ax.transAxes, ha='left', va='top',
            fontsize=plot_style.TEXT)


def scatter_panels(true_lfc, panels, stem, error=False, verbose=True):
    """Estimated LFC (or its error) against the true LFC, one panel per estimator."""
    fig, axes = plt.subplots(1, len(panels), figsize=plot_style.figsize(len(panels)), sharex=True, sharey=True)
    for ax, (title, color, est) in zip(axes, panels):
        ax.scatter(true_lfc, true_lfc - est if error else est, s=2, color=color, alpha=0.5, lw=0)
        if error:
            ax.axhline(0, color='black', ls='--', lw=0.8)
        else:
            lo, hi = np.nanmin(true_lfc), np.nanmax(true_lfc)
            ax.plot([lo, hi], [lo, hi], color='black', ls='--', lw=0.8)
            annotate_r2(ax, true_lfc, est)
        ax.set_title(title)
        ax.set_xlabel('True LFC')
    axes[0].set_ylabel('Bias' if error else 'Estimated LFC')
    fig.tight_layout()
    save(fig, stem, verbose)


def lfc():
    res = pd.read_csv(require_input(RESULTS / 'CITE_seq_results.csv', what='per-gene LFCs of every replicate',
                                    source='python -m citeseq.exp'), index_col=0)
    sig = res[res.is_signal_gene]
    panels = [(method_label('LN'), method_color('LN'), sig.ln_lfc),
              (SCANPY_LOG1P, method_color('t-test'), sig.scanpy_lfc)]
    scatter_panels(sig.true_lfc, panels, 'CITE_seq_lfc_plot')
    scatter_panels(sig.true_lfc, panels, 'CITE_seq_lfc_error_plot', error=True)

    seurat_csv = RESULTS / 'seurat_lfc_results.csv'
    if not seurat_csv.exists():
        print(f'skipped the Seurat figures: no {seurat_csv} (cd citeseq/R && Rscript pbmc10k_seurat.R)')
        return
    seurat = pd.read_csv(seurat_csv)
    fig, ax = plt.subplots(figsize=plot_style.figsize())
    scatter = seurat[seurat.true_lfcs != 0]
    ax.scatter(scatter.true_lfcs, scatter.seurat_lfc, s=2, color=method_color('seurat'), alpha=0.5, lw=0)
    lo, hi = scatter.true_lfcs.min(), scatter.true_lfcs.max()
    ax.plot([lo, hi], [lo, hi], color='black', ls='--', lw=0.8)
    annotate_r2(ax, scatter.true_lfcs, scatter.seurat_lfc)
    ax.set_xlabel('True LFC')
    ax.set_ylabel('Estimated LFC')
    ax.set_title(method_label('seurat'))
    fig.tight_layout()
    save(fig, 'CITE_seq_seurat_lfc_plot')

    merged = res.rename_axis('gene_names').reset_index().rename(columns={'replicate': 'rep'}).merge(
        seurat[['gene_names', 'seurat_lfc', 'rep']], on=['gene_names', 'rep'])
    merged = merged[merged.true_lfc != 0]
    panels = [(method_label('LN'), method_color('LN'), merged.ln_lfc),
              (SCANPY_LOG1P, method_color('t-test'), merged.scanpy_lfc),
              (method_label('seurat'), method_color('seurat'), merged.seurat_lfc)]
    scatter_panels(merged.true_lfc, panels, 'CITE_seq_lfc_plot_all_methods')
    scatter_panels(merged.true_lfc, panels, 'CITE_seq_lfc_error_plot_all', error=True)


if __name__ == '__main__':
    plot_style.use()
    gating()
    lfc()
