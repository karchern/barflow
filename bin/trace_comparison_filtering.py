#!/usr/bin/env python3

import argparse
import json
from collections import Counter, defaultdict

import pandas as pd

from comparison_filtering_logic import as_string_list, derive_status, resolve_selectors


def derive_drop_reason(pre_treat, pre_ctrl, post_treat, post_ctrl):
    if not pre_treat and not pre_ctrl:
        return 'no_treatments_or_controls_matched_before_barseq_qc'
    if not pre_treat:
        return 'no_treatments_matched_before_barseq_qc'
    if not pre_ctrl:
        return 'no_controls_matched_before_barseq_qc'
    if not post_treat and not post_ctrl:
        return 'all_treatments_and_controls_failed_barseq_qc'
    if not post_treat:
        return 'all_treatments_failed_barseq_qc'
    if not post_ctrl:
        return 'all_controls_failed_barseq_qc'
    if set(pre_treat) != set(post_treat) or set(pre_ctrl) != set(post_ctrl):
        return 'comparison_retained_with_partial_sample_loss_after_barseq_qc'
    return 'comparison_retained_without_sample_loss'


def format_sample_reason(sample_id, sample_metrics):
    info = sample_metrics.get(sample_id, {})
    qc_failure_reason = info.get('qc_failure_reason', '') or 'unknown'
    total_reads = info.get('total_reads', 'NA')
    median_count = info.get('median_count', 'NA')
    return f"{sample_id}[reason={qc_failure_reason};total_reads={total_reads};median_count={median_count}]"


def build_summary_text(summary_df, reason_counts_df, sample_contrib_df):
    total = len(summary_df)
    dropped = int(summary_df['dropped_after_barseq_qc'].sum())
    retained = total - dropped

    lines = []
    lines.append('============================================')
    lines.append('Comparison Filtering Trace')
    lines.append('============================================')
    lines.append(f'Total comparisons examined: {total}')
    lines.append(f'Comparisons retained after BarSeq QC: {retained}')
    lines.append(f'Comparisons dropped after BarSeq QC: {dropped}')
    lines.append('')
    lines.append('Dropped comparison reasons:')

    dropped_reason_counts = reason_counts_df[reason_counts_df['dropped_after_barseq_qc'] == 1]
    if dropped_reason_counts.empty:
        lines.append('  none')
    else:
        for _, row in dropped_reason_counts.sort_values(['n_comparisons', 'comparison_drop_reason'], ascending=[False, True]).iterrows():
            lines.append(f"  {row['comparison_drop_reason']}: {row['n_comparisons']}")

    lines.append('')
    lines.append('Samples most often responsible for dropped comparisons:')
    if sample_contrib_df.empty:
        lines.append('  none')
    else:
        for _, row in sample_contrib_df.head(10).iterrows():
            lines.append(
                f"  {row['sample_id']}: {row['n_dropped_comparisons']} dropped comparisons | "
                f"qc_failure_reason={row['qc_failure_reason']} | "
                f"total_reads={row['total_reads']} | median_count={row['median_count']}"
            )

    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description='Trace why comparisons are dropped after BarSeq QC.')
    parser.add_argument('--comparisons-json', required=True)
    parser.add_argument('--sample-metrics', required=True)
    parser.add_argument('--min-read-sum-for-qc', type=int, required=True)
    parser.add_argument('--min-median-barcode-count', type=float, required=True)
    parser.add_argument('--output-comparison-trace', required=True)
    parser.add_argument('--output-comparison-sample-trace', required=True)
    parser.add_argument('--output-drop-reason-counts', required=True)
    parser.add_argument('--output-summary-text', required=True)
    args = parser.parse_args()

    with open(args.comparisons_json) as handle:
        comparisons = json.load(handle)

    metrics = pd.read_csv(args.sample_metrics, sep=None, engine='python')
    metrics['sample_id'] = metrics['sample_id'].astype(str)

    if 'qc_passed' not in metrics.columns:
        metrics['qc_passed'] = (
            (metrics['total_reads'] >= args.min_read_sum_for_qc)
            & (metrics['median_count'] >= args.min_median_barcode_count)
        ).astype(int)

    if 'qc_failure_reason' not in metrics.columns:
        failure_reasons = []
        for _, row in metrics.iterrows():
            reasons = []
            if row['total_reads'] < args.min_read_sum_for_qc:
                reasons.append('total_reads_below_threshold')
            if row['median_count'] < args.min_median_barcode_count:
                reasons.append('median_barcode_count_below_threshold')
            failure_reasons.append(';'.join(reasons))
        metrics['qc_failure_reason'] = failure_reasons

    all_ids = sorted(metrics['sample_id'].tolist())
    passed_ids = set(metrics.loc[metrics['qc_passed'] == 1, 'sample_id'].tolist())
    sample_metrics = metrics.set_index('sample_id').to_dict(orient='index')

    comparison_rows = []
    sample_rows = []
    dropped_sample_counter = Counter()
    dropped_sample_details = {}

    for comparison in comparisons['comparisons']:
        name = comparison['name']
        treat_sel = as_string_list(comparison.get('treatments'))
        ctrl_sel = as_string_list(comparison.get('controls'))
        treat_neg_sel = as_string_list(comparison.get('treatments_negative_selection'))
        ctrl_neg_sel = as_string_list(comparison.get('controls_negative_selection'))

        pre_treat_info = resolve_selectors(treat_sel, all_ids, treat_neg_sel)
        pre_ctrl_info = resolve_selectors(ctrl_sel, all_ids, ctrl_neg_sel)

        pre_treat = pre_treat_info['resolved_ids']
        pre_ctrl = pre_ctrl_info['resolved_ids']
        post_treat = sorted([sample_id for sample_id in pre_treat if sample_id in passed_ids])
        post_ctrl = sorted([sample_id for sample_id in pre_ctrl if sample_id in passed_ids])

        removed_treat = sorted(set(pre_treat) - set(post_treat))
        removed_ctrl = sorted(set(pre_ctrl) - set(post_ctrl))

        status_before = derive_status(pre_treat, pre_ctrl)
        status_after = derive_status(post_treat, post_ctrl)
        drop_reason = derive_drop_reason(pre_treat, pre_ctrl, post_treat, post_ctrl)
        dropped_after_qc = int(status_before == 'OK' and (not post_treat or not post_ctrl))

        for arm, pre_ids, post_ids, removed_ids in [
            ('treatment', pre_treat, post_treat, removed_treat),
            ('control', pre_ctrl, post_ctrl, removed_ctrl),
        ]:
            for sample_id in pre_ids:
                info = sample_metrics.get(sample_id, {})
                dropped_by_qc = int(sample_id in removed_ids)
                sample_rows.append({
                    'comparison': name,
                    'arm': arm,
                    'sample_id': sample_id,
                    'present_before_barseq_qc': 1,
                    'present_after_barseq_qc': int(sample_id in post_ids),
                    'dropped_by_barseq_qc': dropped_by_qc,
                    'sample_qc_passed': info.get('qc_passed', ''),
                    'sample_qc_failure_reason': info.get('qc_failure_reason', ''),
                    'total_reads': info.get('total_reads', ''),
                    'median_count': info.get('median_count', ''),
                })

                if dropped_after_qc and dropped_by_qc:
                    dropped_sample_counter[sample_id] += 1
                    dropped_sample_details[sample_id] = {
                        'qc_failure_reason': info.get('qc_failure_reason', ''),
                        'total_reads': info.get('total_reads', ''),
                        'median_count': info.get('median_count', ''),
                    }

        comparison_rows.append({
            'comparison': name,
            'status_before_barseq_qc': status_before,
            'status_after_barseq_qc': status_after,
            'comparison_drop_reason': drop_reason,
            'dropped_after_barseq_qc': dropped_after_qc,
            'treatment_selectors': ';'.join(treat_sel),
            'control_selectors': ';'.join(ctrl_sel),
            'treatment_negative_selectors': ';'.join(treat_neg_sel),
            'control_negative_selectors': ';'.join(ctrl_neg_sel),
            'n_treatments_before_barseq_qc': len(pre_treat),
            'n_controls_before_barseq_qc': len(pre_ctrl),
            'n_treatments_after_barseq_qc': len(post_treat),
            'n_controls_after_barseq_qc': len(post_ctrl),
            'matched_treatments_before_barseq_qc': ';'.join(pre_treat),
            'matched_controls_before_barseq_qc': ';'.join(pre_ctrl),
            'matched_treatments_after_barseq_qc': ';'.join(post_treat),
            'matched_controls_after_barseq_qc': ';'.join(post_ctrl),
            'treatments_removed_by_barseq_qc': ';'.join(removed_treat),
            'controls_removed_by_barseq_qc': ';'.join(removed_ctrl),
            'treatment_removed_qc_details': '; '.join(format_sample_reason(sample_id, sample_metrics) for sample_id in removed_treat),
            'control_removed_qc_details': '; '.join(format_sample_reason(sample_id, sample_metrics) for sample_id in removed_ctrl),
            'treatment_explicit_selectors_not_found_before_barseq_qc': ';'.join(pre_treat_info['explicit_not_found']),
            'control_explicit_selectors_not_found_before_barseq_qc': ';'.join(pre_ctrl_info['explicit_not_found']),
            'treatment_negative_selectors_not_found_before_barseq_qc': ';'.join(pre_treat_info['negative_explicit_not_found']),
            'control_negative_selectors_not_found_before_barseq_qc': ';'.join(pre_ctrl_info['negative_explicit_not_found']),
            'treatments_excluded_by_negative_selection': ';'.join(pre_treat_info['excluded_by_negative_selection']),
            'controls_excluded_by_negative_selection': ';'.join(pre_ctrl_info['excluded_by_negative_selection']),
        })

    comparison_df = pd.DataFrame(comparison_rows).sort_values('comparison')
    sample_df = pd.DataFrame(sample_rows).sort_values(['comparison', 'arm', 'sample_id'])

    reason_counts_df = (
        comparison_df
        .groupby(['comparison_drop_reason', 'dropped_after_barseq_qc'], dropna=False)
        .size()
        .reset_index(name='n_comparisons')
        .sort_values(['dropped_after_barseq_qc', 'n_comparisons', 'comparison_drop_reason'], ascending=[False, False, True])
    )

    sample_contrib_rows = []
    for sample_id, count in dropped_sample_counter.most_common():
        details = dropped_sample_details[sample_id]
        sample_contrib_rows.append({
            'sample_id': sample_id,
            'n_dropped_comparisons': count,
            'qc_failure_reason': details['qc_failure_reason'],
            'total_reads': details['total_reads'],
            'median_count': details['median_count'],
        })
    sample_contrib_df = pd.DataFrame(sample_contrib_rows)

    comparison_df.to_csv(args.output_comparison_trace, sep='\t', index=False)
    sample_df.to_csv(args.output_comparison_sample_trace, sep='\t', index=False)
    reason_counts_df.to_csv(args.output_drop_reason_counts, sep='\t', index=False)

    summary_text = build_summary_text(comparison_df, reason_counts_df, sample_contrib_df)
    with open(args.output_summary_text, 'w') as handle:
        handle.write(summary_text)


if __name__ == '__main__':
    main()