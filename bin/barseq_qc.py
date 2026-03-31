#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import argparse

# compute_barcode_count_sample_metrics.py \\
#     --input ${counts_path} \\
#     --output ${sample_id}.barcode_matrics.tsv

# Parse this above
parser = argparse.ArgumentParser(description="Compute sample-level metrics from barcode count files.")
parser.add_argument('--sample_id', type=str, required=True, help="the smaple iD ;)")
parser.add_argument('--input', type=str, required=True, help="2-column barcode count matrix with headers")
parser.add_argument('--output', type=str, required=True, help="Output TSV file for sample metrics.")
parser.add_argument('--output_passed', type=str, required=True, help="Output file to write '1' if sample passed QC, '0' otherwise.")
parser.add_argument('--min_read_sum_for_qc', type=int, help="Minimum total read count for sample to pass QC.")

args = parser.parse_args()
counts_file_path = args.input

rows = []
counts_file_path = Path(counts_file_path)
df = pd.read_csv(counts_file_path, sep=',')
df.columns = ['barcode', 'count']

total_reads  = int(df['count'].sum())
n_barcodes   = int((df['count'] > 0).sum())
max_count    = int(df['count'].max()) if len(df) else 0
mean_count   = float(df['count'].mean()) if len(df) else 0.0
median_count = float(df['count'].median()) if len(df) else 0.0

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
    'median_count': median_count
})

out = pd.DataFrame(rows)
out.to_csv(args.output, sep=',', index=False)
