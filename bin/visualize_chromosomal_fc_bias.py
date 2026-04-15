#!/usr/bin/env python3
"""
Visualize chromosomal log2 fold-change (FC) bias from mbarq results.

Takes mbarq barcodes results CSV and good_barcodes mapping (barcode->contig,position),
determines the main chromosomal contig using the same logic as barseq_qc, then plots
log2FC values in a smoothed window over the main chromosome.
"""

import pandas as pd
import numpy as np
import argparse
from pathlib import Path

def determine_main_contig(df_with_good_barcodes_and_contigs):
    """Determine the main contig using the same logic as barseq_qc.py.
    
    - If 1 contig: use it
    - If 2 contigs: check if one is <10% of other AND >1Mb → use larger (assumes smaller is plasmid)
    - If >2 contigs OR conditions not met: ERROR
    
    Returns: contig_id or None if error
    """
    unique_contigs = set(df_with_good_barcodes_and_contigs['contig'])
    
    if len(unique_contigs) == 1:
        return list(unique_contigs)[0]
    
    if len(unique_contigs) == 2:
        print(
            "It looks like you have more than one contig... Checking if we have 2 and one is much longer than the other, "
            "which would suggest that the shorter one is a plasmid and we can proceed with the longer one."
        )
        contig_lengths = df_with_good_barcodes_and_contigs.groupby('contig')['position'].max()
        
        if contig_lengths.min() < 0.1 * contig_lengths.max() and contig_lengths.max() > 1E6:
            main_contig = contig_lengths.idxmax()
            print(
                f"One contig is much longer than the other, and the longer one is > 1Mb. "
                f"Assuming the shorter one ({contig_lengths.idxmin()}) is a plasmid and proceeding with the longer one ({main_contig})."
            )
            return main_contig
        else:
            raise ValueError(
                "We have 2 contigs but the length difference is not large enough to assume one is a plasmid. "
                "Cannot proceed with chromosomal FC visualization."
            )
    
    if len(unique_contigs) > 2:
        raise ValueError(
            f"We have {len(unique_contigs)} contigs. Cannot proceed with chromosomal FC visualization "
            "since it relies on having a single main contig."
        )


def smooth_fc_values(positions, lfc_values, window=10000):
    """Compute smoothed FC values using sliding-window mean.
    
    Args:
        positions: list of genomic positions
        lfc_values: list of log-fold-change values at each position
        window: window size for smoothing (bp)
    
    Returns:
        (bin_centers, smoothed_lfc, observed_sites_per_window)
    """
    if not positions:
        return [], [], []
    
    positions = np.asarray(positions, dtype=float)
    lfc_values = np.asarray(lfc_values, dtype=float)
    
    # Sort by position
    sort_idx = np.argsort(positions)
    positions = positions[sort_idx]
    lfc_values = lfc_values[sort_idx]
    
    min_pos = int(positions.min())
    max_pos = int(positions.max())
    
    # Slide window along genomic coordinates
    window_starts = list(range(min_pos, max_pos + 1, window))
    
    moving_mean_array = []
    bin_centers_list = []
    observed_sites_per_window = []
    
    left_idx = 0
    right_idx = 0
    running_sum = 0.0
    running_n = 0
    
    for start in window_starts:
        end = start + window
        
        # Add positions that fall in [start, end)
        while right_idx < len(positions) and positions[right_idx] < end:
            running_sum += lfc_values[right_idx]
            running_n += 1
            right_idx += 1
        
        # Remove positions that fall before start
        while left_idx < len(positions) and positions[left_idx] < start:
            running_sum -= lfc_values[left_idx]
            running_n -= 1
            left_idx += 1
        
        if running_n > 0:
            moving_mean = running_sum / running_n
        else:
            moving_mean = 0.0
        
        moving_mean_array.append(moving_mean)
        bin_centers_list.append(start + (window / 2))
        observed_sites_per_window.append(running_n)
    
    return bin_centers_list, moving_mean_array, observed_sites_per_window


def main():
    parser = argparse.ArgumentParser(
        description="Visualize chromosomal log2 fold-change bias from mbarq results."
    )
    parser.add_argument('--mbarq-barcodes-results', type=str, required=True,
                        help="mbarq barcodes_results CSV (with barcode, LFC, Name, etc.)")
    parser.add_argument('--good-barcodes', type=str, required=True,
                        help="Good barcodes CSV (with barcode, contig, position columns)")
    parser.add_argument('--comparison-id', type=str, required=True,
                        help="Comparison identifier for labels")
    parser.add_argument('--window', type=int, default=10000,
                        help="Smoothing window size in bp (default: 10000)")
    parser.add_argument('--output-plot', type=str, required=True,
                        help="Output plot filename (PNG)")
    parser.add_argument('--output-data', type=str, required=True,
                        help="Output data CSV (bin centers, smoothed LFC, count per window)")
    
    args = parser.parse_args()
    
    # Load data
    print(f"Loading mbarq barcodes results from {args.mbarq_barcodes_results}...")
    mbarq_df = pd.read_csv(args.mbarq_barcodes_results)
    
    print(f"Loading good barcodes from {args.good_barcodes}...")
    good_bc_df = pd.read_csv(args.good_barcodes)
    
    # Merge mbarq results with good_barcodes on barcode
    print("Merging mbarq results with good barcodes...")
    merged_df = mbarq_df.merge(good_bc_df[['barcode', 'contig', 'position']], on='barcode', how='inner')
    
    if merged_df.empty:
        raise ValueError("No matching barcodes found between mbarq results and good barcodes.")
    
    print(f"Found {len(merged_df)} barcodes with LFC and contig/position info.")
    
    # Determine main contig
    print("Determining main contig...")
    main_contig = determine_main_contig(merged_df)
    
    # Filter to main contig
    merged_df = merged_df[merged_df['contig'] == main_contig].copy()
    print(f"Using contig {main_contig} with {len(merged_df)} mapped barcodes.")
    
    if merged_df.empty:
        raise ValueError(f"No barcodes found on main contig {main_contig}.")
    
    # Get position and LFC
    positions = merged_df['position'].tolist()
    lfc_values = merged_df['LFC'].tolist()
    
    # Compute smoothed FC values
    print(f"Computing smoothed FC values (window={args.window} bp)...")
    bin_centers, smoothed_lfc, observed_sites = smooth_fc_values(positions, lfc_values, window=args.window)
    
    if not bin_centers:
        raise ValueError("No windows with sufficient data for smoothing.")
    
    print(f"Total number of bins: {len(bin_centers)}")
    print(f"Bins with sufficient mapped reads: {sum(1 for x in observed_sites if x > 0)}")
    
    # Save data
    print(f"Saving data to {args.output_data}...")
    data_df = pd.DataFrame({
        'bin_center': bin_centers,
        'smoothed_lfc': smoothed_lfc,
        'observed_barcodes_in_window': observed_sites
    })
    data_df.to_csv(args.output_data, index=False)
    
    # Create plot
    print(f"Creating plot...")
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        fig, ax = plt.subplots(figsize=(14, 5))
        
        # Scatter plot with color based on count of barcodes per window
        scatter = ax.scatter(bin_centers, smoothed_lfc, s=20, c=observed_sites, 
                            cmap='viridis', alpha=0.7, edgecolors='none')
        
        # Add zero line
        ax.axhline(y=0, color='red', linestyle='--', linewidth=1, alpha=0.5, label='LFC=0')
        
        ax.set_xlabel('Genomic position (bp)')
        ax.set_ylabel('Smoothed log2 fold-change')
        ax.set_title(f'{args.comparison_id} - Chromosomal FC bias (window={args.window} bp)')
        
        cbar = fig.colorbar(scatter, ax=ax)
        cbar.set_label('Observed barcodes in window')
        
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        fig.savefig(args.output_plot, dpi=150, bbox_inches='tight')
        plt.close(fig)
        
        print(f'Saved plot to {args.output_plot}')
    except Exception as e:
        print(f'Warning: Could not create plot: {e}')
        Path(args.output_plot).touch()
    
    print("Done!")


if __name__ == '__main__':
    main()
