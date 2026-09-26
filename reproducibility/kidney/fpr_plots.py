import os
import sys
import numpy as np
import pandas as pd

from paths import fig_name, require_input, results_dir
import plot_style
from method_colors import method_color, method_draw_order, method_label, method_order
import matplotlib.pyplot as plt
import argparse

plot_style.use()


def plot_fpr_results(csv_file, output_file, aspect_ratio=1.0, title_suffix=None,
                     p_label='Sampling probability'):
    """
    Plot log10(FPR) vs p from CSV results file.
    
    Parameters
    ----------
    csv_file : str
        Path to input CSV file
    output_file : str
        Path to output plot file
    aspect_ratio : float
        Aspect ratio for the plot (default: 1.0)
    title_suffix : str, optional
        Text to append to the plot title (default: None)
    p_label : str
        What p is: 'Sampling probability' for UMI downsampling, 'Split probability' for spot splitting
    """
    # Read the CSV file
    df = pd.read_csv(require_input(
        csv_file,
        what="this arm's FPR-vs-subsampling results",
        source="run `python -m kidney.umi_null` first",
    ))
    
    df['method_display'] = df['method'].map(method_label)
    methods = method_draw_order(df['method_display'])
    p_values = sorted(df['p'].unique())

    # Create the plot
    fig, ax = plt.subplots(1, 1, figsize=plot_style.figsize())
    
    # Plot each method
    # Use raw FPR values and plot on log scale to avoid epsilon capping issues
    for method in methods:
        method_data = df[df['method_display'] == method]
        method_data = method_data.sort_values('p')
        
        fpr_means = method_data['fpr_mean'].values
        fpr_stds = method_data['fpr_std'].values
        p_vals = method_data['p'].values
        
        # Handle zero FPR values by setting a minimum for plotting
        # Use asymmetric error bars: only show positive error bars to avoid going below log scale minimum
        epsilon = 1e-10
        fpr_means_plot = np.maximum(fpr_means, epsilon)  # Minimum for log scale
        
        # Compute error bar bounds in linear space
        fpr_lower = np.maximum(fpr_means - fpr_stds, epsilon)
        fpr_upper = np.maximum(fpr_means + fpr_stds, epsilon)  # Also cap upper bound
        
        # Convert to log10 for plotting
        fpr_means_log10 = np.log10(fpr_means_plot)
        fpr_lower_log10 = np.log10(fpr_lower)
        fpr_upper_log10 = np.log10(fpr_upper)
        
        # Asymmetric error bars: ensure non-negative
        # Lower error: how far down from mean (should be non-negative since fpr_lower <= fpr_means)
        yerr_lower = np.maximum(0, fpr_means_log10 - fpr_lower_log10)
        # Upper error: how far up from mean (should be non-negative since fpr_upper >= fpr_means)
        yerr_upper = np.maximum(0, fpr_upper_log10 - fpr_means_log10)
        
        ax.errorbar(p_vals, fpr_means_log10, 
                   yerr=[yerr_lower, yerr_upper], 
                   label=method, 
                   color=method_color(method), **plot_style.ERRORBAR)
    
    # Prepare title text
    title = f'FPR vs {p_label}'
    if title_suffix:
        title = f'{title} {title_suffix}'
    
    ax.set_xlabel(f'{p_label} $p$')
    ax.set_ylabel(r'$\log_{10}(\mathrm{FPR})$')
    ax.set_title(title)
    ax.tick_params(axis='both', which='major')
    ax.set_aspect(aspect_ratio, adjustable='box')

    plt.tight_layout()
    plot_style.legend_outside(fig)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_file}")
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot FPR results from CSV file')
    parser.add_argument('--csv_file', type=str, required=True,
                        help='Input CSV file with FPR results')
    parser.add_argument('--output', type=str, default=str(results_dir(__file__) / fig_name('fpr_plot')),
                        help='Output file for the plot (default: kidney/results/fpr_plot.pdf)')
    parser.add_argument('--aspect_ratio', type=float, default=1.0,
                        help='Aspect ratio for the plot (default: 1.0)')
    parser.add_argument('--title_suffix', type=str, default=None,
                        help='Text to append to the plot title (default: None)')
    parser.add_argument('--p_label', type=str, default='Sampling probability',
                        help="Name for p in the title and x-axis; 'Split probability' for spot_split (default: Sampling probability)")
    args = parser.parse_args()
    
    if not os.path.exists(args.csv_file):
        raise FileNotFoundError(f"CSV file not found: {args.csv_file}")
    
    print(f"Reading results from {args.csv_file}...")
    plot_fpr_results(args.csv_file, args.output, args.aspect_ratio, title_suffix=args.title_suffix,
                     p_label=args.p_label)
    print("Done!")

