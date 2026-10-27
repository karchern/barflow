#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(description="Compute sample-level metrics from barcode count files.")
parser.add_argument('--sample_id', type=str, required=True, help="the smaple iD ;)")
parser.add_argument('--input_counts', type=str, required=True, help="2-column barcode count matrix with headers")
parser.add_argument('--sample_goodbarcodes_contig_position_map', type=str, required=True, help="csv with header barcode,locus_tag,contig,position")
parser.add_argument('--output_barcode_metrics', type=str, required=True, help="Output TSV file for sample metrics.")
parser.add_argument('--output_passed', type=str, required=True, help="Output file to write '1' if sample passed QC, '0' otherwise.")
parser.add_argument('--min_read_sum_for_qc', type=int, help="Minimum total read count for sample to pass QC.")
parser.add_argument('--min_mean_barcode_count', type=float, help="Minimum mean barcode count for sample to pass QC.")
parser.add_argument('--output_good_barcode_counts', type=str, required=True, help="Output counts restricted to good barcodes, in 2fast2q format.")

args = parser.parse_args()

rows = []
df_raw = pd.read_csv(Path(args.input_counts), sep=',')
original_feature_col = df_raw.columns[0]  # e.g. '#Feature'
original_count_col   = df_raw.columns[1]  # e.g. the sample name
df = df_raw.copy()
df.columns = ['barcode', 'count']

df_with_good_barcodes = df.join(pd.read_csv(args.sample_goodbarcodes_contig_position_map, sep=',').set_index('barcode'), on='barcode', how='inner')

# Write the counts of good barcodes only, in the original 2fast2q column format
good_out = df_with_good_barcodes[['barcode', 'count']].copy()
good_out['count'] = good_out['count'].round().astype(int)
good_out.columns = [original_feature_col, original_count_col]
good_out.to_csv(args.output_good_barcode_counts, index=False)


total_reads     = int(df['count'].sum())
n_barcodes_detected      = int((df['count'] > 0).sum())
max_count       = int(df['count'].max()) if len(df) else 0
mean_count      = float(df['count'].mean()) if len(df) else 0.0
median_count    = float(df['count'].median()) if len(df) else 0.0
min_mean_threshold = float(args.min_mean_barcode_count)

qc_passed = (
    total_reads >= int(args.min_read_sum_for_qc)
    and mean_count >= min_mean_threshold
)

qc_failure_reasons = []
if total_reads < int(args.min_read_sum_for_qc):
    qc_failure_reasons.append('total_reads_below_threshold')
if mean_count < min_mean_threshold:
    qc_failure_reasons.append('mean_barcode_count_below_threshold')

qc_failure_reason = ';'.join(qc_failure_reasons) if qc_failure_reasons else ''

if qc_passed:
    with open(args.output_passed, 'w') as f:
        f.write('1')
else:
    with open(args.output_passed, 'w') as f:
        f.write('0')

rows.append({
    'sample_id': args.sample_id,
    'qc_passed': int(qc_passed),
    'qc_failure_reason': qc_failure_reason,
    'qc_min_read_sum_threshold': int(args.min_read_sum_for_qc),
    'qc_min_mean_barcode_count_threshold': min_mean_threshold,
    'total_reads': total_reads,
    'n_barcodes_detected': n_barcodes_detected,
    'max_count': max_count,
    'mean_count': mean_count,
    'median_count': median_count,
})

out = pd.DataFrame(rows)
out.to_csv(args.output_barcode_metrics, sep=',', index=False)
