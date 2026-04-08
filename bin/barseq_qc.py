#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import argparse
import numpy as np


def fit_piecewise_ptr(xs, ys):
    """Fit a circular two-segment piecewise model (down then up).

    We optimize over:
    - a circular cut point (genome coordinate where we unwrap the circle), and
    - a split point in the unwrapped axis.

    This keeps the model as two linear segments on a circular genome while allowing
    apparent 2- or 3-segment shapes on the displayed linear coordinate axis.
    """
    if len(xs) < 6:
        return None

    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    if len(x) < 2:
        return None

    diffs = np.diff(x)
    diffs = diffs[diffs > 0]
    if len(diffs) == 0:
        return None

    # Approximate circular genome length from span + typical bin spacing.
    circumference = float((x[-1] - x[0]) + np.median(diffs))
    if circumference <= 0:
        return None

    n = len(x)

    best = None
    for cut_idx in range(n):
        # Unwrap circle starting at cut_idx.
        x_rot = np.concatenate([x[cut_idx:], x[:cut_idx] + circumference])
        y_rot = np.concatenate([y[cut_idx:], y[:cut_idx]])
        rot_to_sorted_idx = np.concatenate([np.arange(cut_idx, n), np.arange(0, cut_idx)])

        for split_idx in range(2, n - 2):
            split_x = x_rot[split_idx]

            z = x_rot - x_rot[0]
            split_z = float(z[split_idx])

            if split_z <= 0 or split_z >= circumference:
                continue

            # Circular-continuous piecewise model (no jump at split and no jump at wrap):
            # y(z) = a + b*z + c*max(0, z-s), with y(0)=y(C)
            # => b*C + c*(C-s) = 0  => c = -b*C/(C-s)
            # Substitute to fit only [a,b] for each split s.
            denom = circumference - split_z
            if denom <= 0:
                continue

            h = np.maximum(0.0, z - split_z)
            basis = z - (circumference / denom) * h
            X = np.column_stack([np.ones_like(z), basis])
            beta, _, _, _ = np.linalg.lstsq(X, y_rot, rcond=None)
            a, b = beta
            c = -b * circumference / denom
            m_left = float(b)
            m_right = float(b + c)

            # Constrain to PTR shape on the circle: down then up.
            if not (m_left < 0 and m_right > 0):
                continue

            yhat_rot = a + b * z + c * h

            # Map predictions back to sorted original coordinate order.
            yhat_sorted = np.empty_like(y)
            yhat_sorted[rot_to_sorted_idx] = yhat_rot

            ss_res = float(np.sum((y - yhat_sorted) ** 2))
            ss_tot = float(np.sum((y - np.mean(y)) ** 2))
            r2 = 0.0 if ss_tot == 0 else 1.0 - (ss_res / ss_tot)

            if best is None or r2 > best['r2']:
                peak_idx = int(np.argmax(yhat_sorted))
                trough_idx = int(np.argmin(yhat_sorted))
                trough_val = float(yhat_sorted[trough_idx])
                peak_val = float(yhat_sorted[peak_idx])

                if trough_val > 0:
                    ptr = peak_val / trough_val
                else:
                    ptr = None

                left_size = float(split_z)
                right_size = float(circumference - split_z)
                if left_size > 0 and right_size > 0:
                    half_size_ratio = max(left_size, right_size) / min(left_size, right_size)
                else:
                    half_size_ratio = None

                best = {
                    'x_sorted': x,
                    'yhat_sorted': yhat_sorted,
                    'r2': r2,
                    'split_x': float(split_x),
                    'cut_x': float(x[cut_idx]),
                    'm_left': float(m_left),
                    'm_right': float(m_right),
                    'peak': peak_val,
                    'trough': trough_val,
                    'peak_location': int(round(float(x[peak_idx]))),
                    'trough_location': int(round(float(x[trough_idx]))),
                    'ptr': ptr,
                    'half_size_ratio': half_size_ratio,
                }

    return best

# Code was adapted from https://github.com/alipirani88/Growth-rate-estimate/blob/master/modules/generate_PTR_dataframe.py
def smoothing_1(
    position_count_df: pd.DataFrame,
    window,
    logger, 
    Config,
    sample_id
    ):

    moving_mean_array = []
    sufficient_coverage = 0
    bin_centers = []
    observed_sites_per_window = []

    step = window
    if position_count_df.empty:
        print("Total number of bins: 0")
        print("The number of bins with sufficient mapped reads: 0")
        return [[sample_id, None, None, None, None, None], ""]

    position_count_df = position_count_df.sort_values('position').reset_index(drop=True)
    min_pos = int(position_count_df['position'].min())
    max_pos = int(position_count_df['position'].max())

    # Slide window along genomic coordinates, not over a dense per-bp array.
    window_starts = list(range(min_pos, max_pos + 1, step))

    if not window_starts:
        print("Total number of bins: 0")
        print("The number of bins with sufficient mapped reads: 0")
        return [[sample_id, None, None, None, None, None], ""]

    total_bins = len(window_starts)

    positions = position_count_df['position'].tolist()
    counts = position_count_df['count'].tolist()

    left_idx = 0
    right_idx = 0
    running_sum = 0.0
    running_n = 0

    for start in window_starts:
        end = start + window

        while right_idx < len(positions) and positions[right_idx] < end:
            running_sum += counts[right_idx]
            running_n += 1
            right_idx += 1

        while left_idx < len(positions) and positions[left_idx] < start:
            running_sum -= counts[left_idx]
            running_n -= 1
            left_idx += 1

        if running_n > 0:
            moving_mean = running_sum / running_n
            sufficient_coverage += 1
        else:
            moving_mean = 0.0

        moving_mean_array.append(moving_mean)
        bin_centers.append(start + (window / 2))
        observed_sites_per_window.append(running_n)
            
    print(f"Total number of bins: {total_bins}")
    print(f"The number of bins with sufficient mapped reads: {sufficient_coverage}")

    xs_all = [int(x) for x in bin_centers]
    ys_all = [float(y) for y in moving_mean_array]

    fit_full = fit_piecewise_ptr(xs_all, ys_all)
    if fit_full is None:
        PTR_median = None
        peak = None
        through = None
        peak_location = None
        through_location = None
    else:
        PTR_median = fit_full['ptr']
        peak = fit_full['peak']
        through = fit_full['trough']
        peak_location = fit_full['peak_location']
        through_location = fit_full['trough_location']

    a = [
        sample_id, PTR_median,peak_location, through_location,peak,through
    ]

    # --- Plotting ---
    try:
        # Use non-interactive backend for headless environments
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        xs = []
        ys = []
        cs = []
        for ii, val in enumerate(moving_mean_array):
            xs.append(int(bin_centers[ii]))
            ys.append(val)
            cs.append(int(observed_sites_per_window[ii]))

        if xs and ys:
            fig, ax = plt.subplots(figsize=(10, 4))
            scatter = ax.scatter(xs, ys, s=12, c=cs, cmap='viridis', alpha=0.85)

            fit_for_plot = fit_piecewise_ptr(xs, ys)
            r2_text = 'NA'
            half_ratio_text = 'NA'
            ptr_text = 'NA'
            if fit_for_plot is not None:
                ax.plot(
                    fit_for_plot['x_sorted'],
                    fit_for_plot['yhat_sorted'],
                    color='crimson',
                    linewidth=2,
                    label='piecewise fit'
                )
                r2_text = f"{fit_for_plot['r2']:.3f}"
                if fit_for_plot.get('half_size_ratio') is not None:
                    half_ratio_text = f"{fit_for_plot['half_size_ratio']:.3f}"
                if fit_for_plot.get('ptr') is not None:
                    ptr_text = f"{fit_for_plot['ptr']:.3f}"

            ax.set_xlabel('Genomic position (bp)')
            ax.set_ylabel('Sliding-window mean coverage')
            ax.set_title(
                f'{sample_id} sliding_window_mean\n'
                f'PTR={ptr_text} | R2={r2_text} | half_size_ratio={half_ratio_text}'
            )
            cbar = fig.colorbar(scatter, ax=ax)
            cbar.set_label('Observed genomic sites in window')
            if peak_location is not None:
                ax.axvline(peak_location, color='C1', linestyle='--', label='peak')
            if through_location is not None:
                ax.axvline(through_location, color='C2', linestyle=':', label='trough')
            if peak_location is not None or through_location is not None or fit_for_plot is not None:
                ax.legend()
            outfn = f'{sample_id}_median_sliding_window.png'
            plt.tight_layout()
            fig.savefig(outfn, dpi=150)
            plt.close(fig)
            print(f'Saved sliding-window mean plot to {outfn}')

            # 2nd version: keep only windows in the top 50th percentile of observed site counts
            cs_series = pd.Series(cs)
            site_count_threshold = float(cs_series.quantile(0.5))
            xs_top = [x for x, c in zip(xs, cs) if c >= site_count_threshold]
            ys_top = [y for y, c in zip(ys, cs) if c >= site_count_threshold]
            cs_top = [c for c in cs if c >= site_count_threshold]

            if xs_top and ys_top:
                fig2, ax2 = plt.subplots(figsize=(10, 4))
                scatter2 = ax2.scatter(xs_top, ys_top, s=12, c=cs_top, cmap='viridis', alpha=0.85)

                fit_top = fit_piecewise_ptr(xs_top, ys_top)
                r2_top_text = 'NA'
                half_ratio_top_text = 'NA'
                ptr_top_text = 'NA'
                if fit_top is not None:
                    ax2.plot(
                        fit_top['x_sorted'],
                        fit_top['yhat_sorted'],
                        color='crimson',
                        linewidth=2,
                        label='piecewise fit'
                    )
                    r2_top_text = f"{fit_top['r2']:.3f}"
                    if fit_top.get('half_size_ratio') is not None:
                        half_ratio_top_text = f"{fit_top['half_size_ratio']:.3f}"
                    if fit_top.get('ptr') is not None:
                        ptr_top_text = f"{fit_top['ptr']:.3f}"

                ax2.set_xlabel('Genomic position (bp)')
                ax2.set_ylabel('Sliding-window mean coverage')
                ax2.set_title(
                    f'{sample_id} sliding_window_mean top50pct_sites\n'
                    f'PTR={ptr_top_text} | R2={r2_top_text} | half_size_ratio={half_ratio_top_text}'
                )
                cbar2 = fig2.colorbar(scatter2, ax=ax2)
                cbar2.set_label('Observed genomic sites in window')
                if fit_top is not None:
                    ax2.legend()
                outfn2 = f'{sample_id}_median_sliding_window_top50pct_sites.png'
                plt.tight_layout()
                fig2.savefig(outfn2, dpi=150)
                plt.close(fig2)
                print(f'Saved top-50pct-sites sliding-window mean plot to {outfn2}')
            else:
                print('No windows passed top-50pct-sites filter for second plot')
        else:
            print('No valid points to plot for moving_mean_array')
    except Exception as e:
        print(f'Could not create sliding-window mean plot: {e}')

    b = []
    header = "bin,bin_genomic_center,count\n"
    count = 0
    for ii, i in enumerate(moving_mean_array):
        bin_genomic_center = int(bin_centers[ii])
        count += 1
        b.append(sample_id+','+str(count)+','+str(bin_genomic_center)+','+str(i)+'\n')
    b = ''.join(b)

    return(
        [a, b]
    )


# compute_barcode_count_sample_metrics.py \\
#     --input ${counts_path} \\
#     --output ${sample_id}.barcode_matrics.tsv

# Parse this above
parser = argparse.ArgumentParser(description="Compute sample-level metrics from barcode count files.")
parser.add_argument('--sample_id', type=str, required=True, help="the smaple iD ;)")
parser.add_argument('--input_counts', type=str, required=True, help="2-column barcode count matrix with headers")
parser.add_argument('--sample_goodbarcodes_contig_position_map', type=str, required=True, help="TODO") 
parser.add_argument('--output_barcode_metrics', type=str, required=True, help="Output TSV file for sample metrics.")
#parser.add_argument('--output_PTR', type=str, required=True, help="TODO")
parser.add_argument('--output_median_of_means_over_genomes', type=str, required=True, help="TODO")
parser.add_argument('--output_passed', type=str, required=True, help="Output file to write '1' if sample passed QC, '0' otherwise.")
parser.add_argument('--min_read_sum_for_qc', type=int, help="Minimum total read count for sample to pass QC.")
parser.add_argument('--min_median_barcode_count', type=int, help="Minimum median barcode count for sample to pass QC.")

args = parser.parse_args()
counts_file_path = args.input_counts

rows = []
counts_file_path = Path(counts_file_path)
df = pd.read_csv(counts_file_path, sep=',')
df.columns = ['barcode', 'count']

df_with_good_barcodes_and_contigs = df.join(pd.read_csv(args.sample_goodbarcodes_contig_position_map, sep=',').set_index('barcode'), on='barcode', how='inner')

if len(set(df_with_good_barcodes_and_contigs['contig'])) > 1:
    print(
        """
        It looks like you have more than one contig... Checking if we have 2 and one is much longer than the other, 
        which would suggest that the shorter one is a plasmid and we can proceed with the longer one. If not, we cannot 
        proceed with PTR calculation since it relies on having a single contig.
        """
    )
    if len(set(df_with_good_barcodes_and_contigs['contig'])) == 2:
        print(
            "We have 2 contigs. Checking if one is much longer than the other..."
        )
        contig_lengths = df_with_good_barcodes_and_contigs.groupby('contig')['position'].max()
        if contig_lengths.min() < 0.1 * contig_lengths.max() and contig_lengths.max() > 1E6:
            print(
                "One contig is much longer than the other, and the longer one is > 1Mb. Assuming the shorter one is a plasmid and proceeding with the longer one."
            )
            df_with_good_barcodes_and_contigs = df_with_good_barcodes_and_contigs[df_with_good_barcodes_and_contigs['contig'] == contig_lengths.idxmax()]
    else:
        exit(
            "We have more than 2 contigs, or the length difference is not large enough to assume one is a plasmid. Cannot proceed with PTR calculation since it relies on having a single contig."
        )


df_with_good_barcodes_and_contigs_summed_by_pos = (
    df_with_good_barcodes_and_contigs
    .groupby('position', as_index=True)['count']
    .sum()
    .to_frame('count')
    .reset_index()
)

PTR_results = smoothing_1(
    position_count_df=df_with_good_barcodes_and_contigs_summed_by_pos,
    window=10000,
    logger=None,
    Config=None,
    sample_id=args.sample_id
)

with open(args.output_median_of_means_over_genomes, 'w') as f:
    f.write(str(PTR_results[1]))


total_reads     = int(df['count'].sum())
n_barcodes_detected      = int((df['count'] > 0).sum())
max_count       = int(df['count'].max()) if len(df) else 0
mean_count      = float(df['count'].mean()) if len(df) else 0.0
median_count    = float(df['count'].median()) if len(df) else 0.0
ptr             = PTR_results[0][1]
peak_location   = PTR_results[0][2]
trough_location = PTR_results[0][3]
peak            = PTR_results[0][4]
trough          = PTR_results[0][5]

qc_passed = (
    total_reads >= int(args.min_read_sum_for_qc)
    and median_count >= int(args.min_median_barcode_count)
)

if qc_passed:
    with open(args.output_passed, 'w') as f:
        f.write('1')
else:
    with open(args.output_passed, 'w') as f:
        f.write('0')

rows.append({
    'sample_id': args.sample_id,
    'total_reads': total_reads,
    'n_barcodes_detected': n_barcodes_detected,
    'max_count': max_count,
    'mean_count': mean_count,
    'median_count': median_count,
    'ptr': ptr,
    'peak_location': peak_location,
    'trough_location': trough_location,
    'peak': peak,
    'trough': trough
})

out = pd.DataFrame(rows)
out.to_csv(args.output_barcode_metrics, sep=',', index=False)
