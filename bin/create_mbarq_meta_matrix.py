#!/usr/bin/env python3
import argparse
import csv

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--treat-ids", nargs="+", required=True,
                   help="Treatment sample IDs")
    p.add_argument("--ctrl-ids", nargs="+", required=True,
                   help="Control sample IDs")
    p.add_argument("--batch", default="batch1",
                   help="Batch name (default: batch1)")
    p.add_argument("-o", "--output", required=True,
                   help="Output CSV path")
    args = p.parse_args()

    rows = []

    # treatments
    for sid in args.treat_ids:
        rows.append((sid, "treated", args.batch))

    # controls
    for sid in args.ctrl_ids:
        rows.append((sid, "control", args.batch))

    with open(args.output, "w", newline="") as fh:
        w = csv.writer(fh, delimiter=",")
        w.writerow(["sampleID", "treatment", "batch"])
        w.writerows(rows)

if __name__ == "__main__":
    main()
