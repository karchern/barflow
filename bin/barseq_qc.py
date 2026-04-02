#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import argparse
import statistics

# Code was adapted from https://github.com/alipirani88/Growth-rate-estimate/blob/master/modules/generate_PTR_dataframe.py
def smoothing_1(
    read_counts: list, 
    window,
    logger, 
    Config,
    sample_id
    ):

    moving_median_array = []
    sufficient_coverage = 0

    step = 100
    median_of_medians_size = 100
    ll = list(range(0, len(read_counts) - window + 1, step))

    if not ll:
        print("Total number of bins: 0")
        print("The number of bins with no mapped reads: 0")
        return

    total_bins = len(ll)

    # initialise first window
    first_start = ll[0]
    first_end = first_start + window
    zeros_in_window = read_counts[first_start:first_end].count(0)

    for idx, i in enumerate(ll):
        start = i
        end = start + window

        # update zeros_in_window when window moves
        if idx > 0:
            # element leaving at start-step
            left_idx = start - step
            right_idx = end - 1
            # subtract zeros that left
            for j in range(left_idx, left_idx + step):
                if read_counts[j] == 0:
                    zeros_in_window -= 1
            # add zeros that entered
            for j in range(end - step, end):
                if j < len(read_counts) and read_counts[j] == 0:
                    zeros_in_window += 1

        # all zeros in this window?
        if zeros_in_window < (0.1*(window)):
            sufficient_coverage += 1
        moving_median_array.append(statistics.median(read_counts[start:end]))
            

    print(f"Total number of bins: {total_bins}")
    print(f"The number of bins with sufficient mapped reads: {sufficient_coverage}")

    ## Step 3
    median_sliding_window_array = []
    for i in range(0, len(moving_median_array), 1):
        start = i
        end = i + median_of_medians_size
        if len(moving_median_array[start:end]) > (median_of_medians_size/2):
            median_sliding_window_array.append(statistics.median(moving_median_array[start:end]))
        else:
            median_sliding_window_array.append(None)

    ## Step 4

    


    peak = max([x for x in median_sliding_window_array if x is not None])  # ignore None
    peak_index = median_sliding_window_array.index(peak)
    if all([(x == 0 or x == None) for x in median_sliding_window_array]):
        #print("All median sliding window values are 0. Cannot calculate PTR.")
        #return
        through = 1
        through_index = None
    else:
        through = min([x for x in median_sliding_window_array if x !=0 and x is not None])  # ignore None and 0
        through_index = median_sliding_window_array.index(through)
    PTR_median = peak/through

    def second_level_index_to_genomic_position(idx, ll, step, median_of_medians_size, window):
            if idx is None:
                return None

            if idx >= len(ll):
                return None

            start_bp = ll[idx]

            # approximate center of second-level smoothing region
            center_offset = ((median_of_medians_size - 1) * step) / 2

            # optional: add half the first-level window width to refer to window center
            center_bp = start_bp + center_offset + (window / 2)

            return int(center_bp)

    peak_location = second_level_index_to_genomic_position(
        peak_index, ll, step, median_of_medians_size, window
    )

    through_location = second_level_index_to_genomic_position(
        through_index, ll, step, median_of_medians_size, window
    )

    a = [
        sample_id, PTR_median,peak_location, through_location,peak,through
    ]


    #with open(out_path_base + "/" + "base.csv", 'w') as out:
    b = []
    header = "bin,bin_genomic_center,count\n"
    count = 0
    for ii, i in enumerate(median_sliding_window_array):
        bin_genomic_center = second_level_index_to_genomic_position(
            ii, ll, step, median_of_medians_size, window
        )
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
parser.add_argument('--output_median_of_medians_over_genomes', type=str, required=True, help="TODO")
parser.add_argument('--output_passed', type=str, required=True, help="Output file to write '1' if sample passed QC, '0' otherwise.")
parser.add_argument('--min_read_sum_for_qc', type=int, help="Minimum total read count for sample to pass QC.")

args = parser.parse_args()
counts_file_path = args.input_counts

rows = []
counts_file_path = Path(counts_file_path)
df = pd.read_csv(counts_file_path, sep=',')
df.columns = ['barcode', 'count']

df_with_good_barcodes_and_contigs = df.join(pd.read_csv(args.sample_goodbarcodes_contig_position_map, sep=',').set_index('barcode'), on='barcode', how='inner')
if len(set(df_with_good_barcodes_and_contigs['contig'])) > 1:
    print(
        "PTR calculation is not possible with fragmented assemblies. Skipping rest of script."
    )

# df_with_good_barcodes_and_contigs does not contain entries for each position.
# Assume that max(position) is the length of the contig, and add all MISSING positions between 0 and that value


df_with_good_barcodes_and_contigs_summed_by_pos = (
    df_with_good_barcodes_and_contigs
    .groupby('position', as_index=True)['count']
    .sum()
    .to_frame('count')
)

ally = pd.DataFrame(
    list(range(0, int(df_with_good_barcodes_and_contigs['position'].max()) + 1)), columns = ['position']
)
ally['count_dummy'] = 0
ally = ally.merge(df_with_good_barcodes_and_contigs_summed_by_pos, on='position', how='left').fillna(0)
ally['count'] = ally['count'] + ally['count_dummy']
ally.drop(columns=['count_dummy'], inplace=True)
ally.sort_values('position', inplace=True)

PTR_results = smoothing_1(
    read_counts=ally['count'].tolist(),
    window=10000,
    logger=None,
    Config=None,
    sample_id=args.sample_id
)

# with open(args.output_PTR, 'w') as f:
#     f.write(str(PTR_results[0]))

with open(args.output_median_of_medians_over_genomes, 'w') as f:
    f.write(str(PTR_results[1]))

total_reads     = int(df['count'].sum())
n_barcodes      = int((df['count'] > 0).sum())
max_count       = int(df['count'].max()) if len(df) else 0
mean_count      = float(df['count'].mean()) if len(df) else 0.0
median_count    = float(df['count'].median()) if len(df) else 0.0
ptr             = PTR_results[0][1]
peak_location   = PTR_results[0][2]
trough_location = PTR_results[0][3]
peak            = PTR_results[0][4]
trough          = PTR_results[0][5]

if total_reads >= int(args.min_read_sum_for_qc):
    with open(args.output_passed, 'w') as f:
        f.write('1')
else:
    with open(args.output_passed, 'w') as f:
        f.write('0')

rows.append({
    'sample_id': args.sample_id,
    'total_reads': total_reads,
    'n_barcodes': n_barcodes,
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
