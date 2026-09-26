import os
import sys
import numpy as np
import pandas as pd

from paths import fig_name, require_input, results_dir
import plot_style
from method_colors import method_color, method_draw_order, method_label, method_order
import json
import matplotlib.pyplot as plt
import argparse
from matplotlib.lines import Line2D
plot_style.use()

def method_and_p_legends(fig, methods, selected_p, linestyle_map):
    """Two keys right of the axes: colour is the method, line style the split probability."""
    plot_style.legend_outside(fig, [Line2D([0], [0], color=method_color(m), lw=1.3) for m in methods],
                              [method_label(m) for m in methods], title='Method')
    plot_style.legend_outside(fig, [Line2D([0], [0], color='black', lw=1.3, linestyle=linestyle_map[p])
                                    for p in selected_p],
                              [f'{p:.3f}' for p in selected_p], title='Split probability $p$', y=0.6)


def plot_tpr_fpr_results(csv_file, output_file):
    """Plot TPR and FPR vs p for all methods."""
    df = pd.read_csv(require_input(
        csv_file,
        what="this arm's spot-splitting results",
        source="run `python -m kidney.spot_split` first",
    ))
    
    # Get unique methods and p values
    methods = method_draw_order(df['method'])
    p_values = sorted(df['p'].unique())
    
    method_to_color = {method: method_color(method) for method in methods}
    
    fig, axes = plt.subplots(1, 2, figsize=plot_style.figsize(2))
    
    # Plot TPR
    for method in methods:
        method_df = df[df['method'] == method]
        method_df = method_df.sort_values('p')
        p_vals = method_df['p'].values
        tpr_means = method_df['tpr_mean'].values
        tpr_stds = method_df['tpr_std'].values
        
        display_name = method_label(method)
        axes[0].errorbar(p_vals, tpr_means, yerr=tpr_stds,
                        label=display_name,
                        color=method_to_color[method], **plot_style.ERRORBAR)
    
    axes[0].set_xlabel(r'Split probability $p$')
    axes[0].set_ylabel('TPR')
    axes[0].set_title('TPR vs Split Probability')
    axes[0].set_ylim([0, 1])
    
    # Plot FPR
    for method in methods:
        method_df = df[df['method'] == method]
        method_df = method_df.sort_values('p')
        p_vals = method_df['p'].values
        fpr_means = method_df['fpr_mean'].values
        fpr_stds = method_df['fpr_std'].values
        
        display_name = method_label(method)
        axes[1].errorbar(p_vals, fpr_means, yerr=fpr_stds,
                        label=display_name,
                        color=method_to_color[method], **plot_style.ERRORBAR)
    
    axes[1].set_xlabel(r'Split probability $p$')
    axes[1].set_ylabel('FPR')
    axes[1].set_title('FPR vs Split Probability')
    axes[1].set_ylim([0, 1])
    
    plt.tight_layout()
    plot_style.legend_outside(fig)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_file}")
    plt.close()


def plot_ap_pr_auc_results(csv_file, output_file):
    """Plot AP and PR-AUC vs p for all methods."""
    df = pd.read_csv(require_input(
        csv_file,
        what="this arm's spot-splitting results",
        source="run `python -m kidney.spot_split` first",
    ))
    
    # Get unique methods and p values
    methods = method_draw_order(df['method'])
    p_values = sorted(df['p'].unique())
    
    method_to_color = {method: method_color(method) for method in methods}
    
    fig, axes = plt.subplots(1, 2, figsize=plot_style.figsize(2))
    
    # Plot AP
    for method in methods:
        method_df = df[df['method'] == method]
        method_df = method_df.sort_values('p')
        p_vals = method_df['p'].values
        ap_means = method_df['ap_mean'].values
        ap_stds = method_df['ap_std'].values
        
        display_name = method_label(method)
        axes[0].errorbar(p_vals, ap_means, yerr=ap_stds,
                         label=display_name,
                         color=method_to_color[method], **plot_style.ERRORBAR)
    axes[0].set_xlabel(r'Split probability $p$')
    axes[0].set_ylabel('Average Precision (AP)')
    axes[0].set_title('AP vs Split Probability')
    axes[0].set_ylim([0, 1])
    
    # Plot PR-AUC
    for method in methods:
        method_df = df[df['method'] == method]
        method_df = method_df.sort_values('p')
        p_vals = method_df['p'].values
        pr_auc_means = method_df['pr_auc_mean'].values
        pr_auc_stds = method_df['pr_auc_std'].values
        
        display_name = method_label(method)
        axes[1].errorbar(p_vals, pr_auc_means, yerr=pr_auc_stds,
                         label=display_name,
                         color=method_to_color[method], **plot_style.ERRORBAR)
    axes[1].set_xlabel(r'Split probability $p$')
    axes[1].set_ylabel('PR-AUC (Trapezoidal)')
    axes[1].set_title('PR-AUC vs Split Probability')
    axes[1].set_ylim([0, 1])
    
    plt.tight_layout()
    plot_style.legend_outside(fig)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_file}")
    plt.close()


def plot_pr_curves_all_p(json_file, output_file, q, lfc):
    """Plot Precision vs Recall curves for each p value (original style)."""
    with open(json_file, 'r') as f:
        curve_data = json.load(f)
    
    p_values = np.array(curve_data['p_values'])
    curves = curve_data['curves']
    
    # Get methods from first p value
    first_p_key = list(curves.keys())[0]
    methods = method_draw_order(curves[first_p_key])
    
    method_to_color = {method: method_color(method) for method in methods}
    
    # Create subplots: one for each p value
    n_p = len(p_values)
    n_cols = min(3, n_p)
    n_rows = (n_p + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=plot_style.figsize(n_cols, n_rows))
    if n_p == 1:
        axes = [axes]
    else:
        axes = axes.flatten() if n_rows > 1 else axes
    
    for p_idx, p in enumerate(p_values):
        ax = axes[p_idx]
        p_key = f"p_{p:.6f}"
        
        if p_key in curves:
            for method in methods:
                if method in curves[p_key]:
                    method_data = curves[p_key][method]
                    recall = np.array(method_data['recall'])
                    precision_mean = np.array(method_data['precision_mean'])
                    precision_std = method_data.get('precision_std')
                    pr_auc_curve_mean = method_data.get('pr_auc_curve_mean', np.nan)
                    
                    # Plot mean curve; its AUC is written in the panel, so one legend serves every panel
                    ax.plot(recall, precision_mean,
                           label=method_label(method),
                           color=method_to_color[method], linewidth=1.3)
                    if not np.isnan(pr_auc_curve_mean):
                        k = method_order(methods).index(method)
                        ax.text(0.03, 0.03 + 0.09 * (len(methods) - 1 - k), f'AUC {pr_auc_curve_mean:.3f}',
                                transform=ax.transAxes, color=method_to_color[method],
                                fontsize=plot_style.SMALL)
                    
                    # Plot std as shaded region (lighter shade)
                    if precision_std is not None:
                        precision_std_arr = np.array(precision_std)
                        ax.fill_between(recall,
                                      precision_mean - precision_std_arr,
                                      precision_mean + precision_std_arr,
                                      alpha=0.1, color=method_to_color[method])
        
        ax.set_xlabel('Recall')
        ax.set_ylabel('Precision')
        ax.set_title(f'p = {p:.3f}')
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
    
    # Hide unused subplots
    for idx in range(n_p, len(axes)):
        axes[idx].axis('off')
    
    fig.suptitle(f'Precision vs recall (q={q}, lfc={lfc})')
    plt.tight_layout()
    plot_style.legend_outside(fig)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_file}")
    plt.close()


def plot_pr_curves_selected_p(json_file, output_file, q, lfc):
    """Plot Precision vs Recall curves for selected p values (smallest, median, largest) in one plot."""
    with open(json_file, 'r') as f:
        curve_data = json.load(f)
    
    p_values = np.array(curve_data['p_values'])
    curves = curve_data['curves']
    
    # Get methods from first p value
    first_p_key = list(curves.keys())[0]
    methods = method_draw_order(curves[first_p_key])
    
    method_to_color = {method: method_color(method) for method in methods}
    
    # Select p values: smallest, largest, and closest to median
    p_sorted = sorted(p_values)
    p_min = p_sorted[0]
    p_max = p_sorted[-1]
    p_median_val = np.median(p_sorted)
    # Find p value closest to median
    p_median_idx = np.argmin(np.abs(p_sorted - p_median_val))
    p_median = p_sorted[p_median_idx]
    
    selected_p = [p_min, p_median, p_max]
    
    fig, ax = plt.subplots(1, 1, figsize=plot_style.figsize())
    
    # Define line styles: dot for smallest, dash for median, solid for largest
    linestyle_map = {p_min: ':', p_median: '--', p_max: '-'}
    
    for method in methods:
        for p in selected_p:
            p_key = f"p_{p:.6f}"
            
            if p_key in curves and method in curves[p_key]:
                method_data = curves[p_key][method]
                recall = np.array(method_data['recall'])
                precision_mean = np.array(method_data['precision_mean'])
                precision_std = method_data.get('precision_std')
                
                # Plot mean curve (no variance shades for selected p plot)
                # Use dot (:) for smallest, dash (--) for median, solid (-) for largest
                ax.plot(recall, precision_mean,
                       color=method_to_color[method], linewidth=1.3, linestyle=linestyle_map[p])
    
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_title(f'Precision vs Recall (q={q}, lfc={lfc})')
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    
    plt.tight_layout()
    method_and_p_legends(fig, methods, selected_p, linestyle_map)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_file}")
    plt.close()


def plot_summary_figure(csv_file, json_file, output_file, q, lfc):
    """
    Create a summary figure with:
      Top row:    [TPR vs p]  [FPR vs p]
      Bottom row: [PR-AUC vs p]  [Selected PR curves (min/median/max p)]
    Keys right of the axes: colour is the method, line style the split probability of the PR curves.
    """
    # Load CSV and JSON
    df = pd.read_csv(require_input(
        csv_file,
        what="this arm's spot-splitting results",
        source="run `python -m kidney.spot_split` first",
    ))
    with open(json_file, 'r') as f:
        curve_data = json.load(f)
    p_values_json = np.array(curve_data['p_values'])
    curves = curve_data['curves']
    # Methods from CSV (ensures alignment with DE summaries)
    methods = method_draw_order(df['method'])
    method_to_color = {method: method_color(method) for method in methods}
    # Select p values: min, median (closest), max from JSON (where PR curves exist)
    p_sorted = np.sort(p_values_json)
    p_min = p_sorted[0]
    p_max = p_sorted[-1]
    p_median_val = np.median(p_sorted)
    p_median = p_sorted[np.argmin(np.abs(p_sorted - p_median_val))]
    selected_p = [p_min, p_median, p_max]
    linestyle_map = {p_min: ':', p_median: '--', p_max: '-'}
    # Build figure
    fig, axes = plt.subplots(2, 2, figsize=plot_style.figsize(2, 2))
    ax_tpr, ax_fpr = axes[0, 0], axes[0, 1]
    ax_prauc, ax_prcurves = axes[1, 0], axes[1, 1]
    # TPR vs p
    for method in methods:
        mdf = df[df['method'] == method].sort_values('p')
        ax_tpr.errorbar(mdf['p'].values, mdf['tpr_mean'].values, yerr=mdf['tpr_std'].values,
                        color=method_to_color[method], **plot_style.ERRORBAR,
                        label=method_label(method))
    ax_tpr.set_xlabel(r'Split probability $p$')
    ax_tpr.set_ylabel('TPR')
    ax_tpr.set_ylim([0, 1])
    ax_tpr.tick_params(axis='both', which='major')
    # FPR vs p
    for method in methods:
        mdf = df[df['method'] == method].sort_values('p')
        ax_fpr.errorbar(mdf['p'].values, mdf['fpr_mean'].values, yerr=mdf['fpr_std'].values,
                        color=method_to_color[method], **plot_style.ERRORBAR)
    ax_fpr.set_xlabel(r'Split probability $p$')
    ax_fpr.set_ylabel('FPR')
    ax_fpr.set_ylim([0, 1])
    ax_fpr.tick_params(axis='both', which='major')
    # PR-AUC vs p (recomputed from JSON PR curves)
    # Use JSON p grid to ensure availability of curves
    p_values_for_auc = np.sort(p_values_json)
    for method in methods:
        pr_auc_list = []
        pr_auc_err = []
        p_valid = []
        for p in p_values_for_auc:
            p_key = f"p_{p:.6f}"
            if p_key in curves and method in curves[p_key]:
                md = curves[p_key][method]
                recall = np.array(md.get('recall', []))
                precision_mean = np.array(md.get('precision_mean', []))
                if recall.size > 1 and precision_mean.size == recall.size:
                    # Recompute AUC from mean curve
                    auc_val = np.trapz(precision_mean, recall)
                    pr_auc_list.append(auc_val)
                    # If std available, use it as yerr; else no error bar
                    pr_auc_std = md.get('pr_auc_curve_std', None)
                    if pr_auc_std is None:
                        pr_auc_err.append(np.nan)
                    else:
                        pr_auc_err.append(pr_auc_std)
                    p_valid.append(p)
        if len(p_valid) > 0:
            p_valid = np.array(p_valid)
            pr_auc_arr = np.array(pr_auc_list, dtype=float)
            pr_auc_err_arr = np.array(pr_auc_err, dtype=float)
            # Plot with error bars where available
            # If all yerr are NaN, omit yerr to avoid warnings
            if np.all(np.isnan(pr_auc_err_arr)):
                ax_prauc.errorbar(p_valid, pr_auc_arr,
                                  color=method_to_color[method], **plot_style.ERRORBAR)
            else:
                ax_prauc.errorbar(p_valid, pr_auc_arr, yerr=pr_auc_err_arr,
                                  color=method_to_color[method], **plot_style.ERRORBAR)
    ax_prauc.set_xlabel(r'Split probability $p$')
    ax_prauc.set_ylabel('PR-AUC')
    ax_prauc.set_ylim([0, 1])
    ax_prauc.tick_params(axis='both', which='major')
    # Selected PR curves (min, median, max p); the keys are drawn once, outside
    for method in methods:
        color = method_to_color[method]
        for p in selected_p:
            p_key = f"p_{p:.6f}"
            if p_key in curves and method in curves[p_key]:
                md = curves[p_key][method]
                recall = np.array(md['recall'])
                precision_mean = np.array(md['precision_mean'])
                ax_prcurves.plot(recall, precision_mean,
                                 color=color, linewidth=1.3, linestyle=linestyle_map[p])
    ax_prcurves.set_xlabel('Recall')
    ax_prcurves.set_ylabel('Precision')
    ax_prcurves.set_xlim([0, 1])
    ax_prcurves.set_ylim([0, 1])
    ax_prcurves.tick_params(axis='both', which='major')
    plt.tight_layout()
    method_and_p_legends(fig, methods, selected_p, linestyle_map)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_file}")
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Visualize results from vishd_test_shape_split_shared.py')
    parser.add_argument('--csv_file', type=str, required=True,
                        help='Input CSV file with results')
    parser.add_argument('--json_file', type=str, required=True,
                        help='Input JSON file with PR curve data')
    parser.add_argument('--output_prefix', type=str, default=str(results_dir(__file__) / 'spot_split'),
                        help='Output file prefix for plots (default: kidney/results/spot_split)')
    args = parser.parse_args()
    
    csv_file = args.csv_file
    json_file = args.json_file
    output_prefix = args.output_prefix
    
    # Read q and lfc from CSV (should be same for all rows)
    df = pd.read_csv(require_input(
        csv_file,
        what="this arm's spot-splitting results",
        source="run `python -m kidney.spot_split` first",
    ))
    q = df['q'].iloc[0]
    lfc = df['lfc'].iloc[0]
    
    print(f"Visualizing results with q={q}, lfc={lfc}")
    
    # Plot 1: TPR and FPR
    print("Plotting TPR and FPR...")
    plot_tpr_fpr_results(csv_file, fig_name(f"{output_prefix}_tpr_fpr"))
    
    # Plot 2: AP and PR-AUC
    print("Plotting AP and PR-AUC...")
    plot_ap_pr_auc_results(csv_file, fig_name(f"{output_prefix}_ap_pr_auc"))
    
    # Plot 3: PR curves for all p values
    print("Plotting PR curves for all p values...")
    plot_pr_curves_all_p(json_file, fig_name(f"{output_prefix}_pr_curves_all_p"), q, lfc)
    
    # Plot 4: PR curves for selected p values (min, median, max)
    print("Plotting PR curves for selected p values...")
    plot_pr_curves_selected_p(json_file, fig_name(f"{output_prefix}_pr_curves_selected_p"), q, lfc)
    
    # Plot 5: Summary figure (TPR/FPR/PR-AUC + selected PR curves)
    print("Plotting summary figure...")
    plot_summary_figure(csv_file, json_file, fig_name(f"{output_prefix}_summary"), q, lfc)
    
    print("Done!")

