#!/usr/bin/env python3

import argparse
import json

from comparison_filtering_logic import build_all_comparison_resolutions


def read_samples(sample_tsv):
    sample_map = {}
    with open(sample_tsv) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            sample_id, sample_path = line.split('\t', 1)
            sample_map[sample_id] = sample_path
    return sample_map


def main():
    parser = argparse.ArgumentParser(description='Resolve comparison selectors against a sample manifest.')
    parser.add_argument('--samples-tsv', required=True)
    parser.add_argument('--comparisons-json', required=True)
    parser.add_argument('--output-json', required=True)
    args = parser.parse_args()

    sample_map = read_samples(args.samples_tsv)
    all_ids = sorted(sample_map.keys())

    with open(args.comparisons_json) as handle:
        comparisons_payload = json.load(handle)

    resolutions = build_all_comparison_resolutions(comparisons_payload, all_ids)

    output = []
    for item in resolutions:
        output.append({
            'name': item['name'],
            'good_barcodes_file': item['good_barcodes_file'],
            'status': item['status'],
            'status_detail': item['status_detail'],
            'treat_list': [[sid, sample_map[sid]] for sid in item['treat_ids'] if sid in sample_map],
            'ctrl_list': [[sid, sample_map[sid]] for sid in item['control_ids'] if sid in sample_map],
        })

    with open(args.output_json, 'w') as handle:
        json.dump(output, handle)


if __name__ == '__main__':
    main()
