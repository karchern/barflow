# Known issues and decisions

Open items first, then items fixed or accepted. Dates are when the item was logged.

## Open

### `remove_all_0_barcodes` never removes anything (logged 2026-10-28, deferred)
- `bin/filter_merged_barcode_matrix.r` (~line 80) drops `#Feature` and `locus_tag` before `all(x == 0)`, but not `contig` and `position`.
  `position` is never 0, so the check is never TRUE and the step is a no-op.
- If a barcode has `NA` in `contig`/`position` and all-zero counts, `all()` returns `NA` and the row indexing produces an all-`NA` row.
- In practice: with `lowly_abundant_barcode_cutoff > 0` the mean-count filter already removes all-zero barcodes
  (their mean in the `filter_on_what` samples is 0), so the bug only matters when the cutoff is 0 and
  `remove_all_0_barcodes = true`.
- Fix: also deselect `contig` and `position`.

### mbarq failures are skipped silently in non-test mode (logged 2026-10-28, kept on purpose for now)
- `withName: mbarq { errorStrategy = 'ignore' }` when `test_run = false`. A comparison whose mbarq task fails
  (including out-of-memory) is missing from the results without an error.
- The task failure is still visible in `trace.txt` and the Nextflow log.

### The good-barcodes file is given twice (design)
- Per sample in `--sample_goodbarcodes_library_map` (used by `barseq_qc`) and per comparison in
  `--comparisons` (`good_barcodes_file`, used by `merge_barcode_matrices` and the chromosomal FC bias plot).
  The locus tags in mbarq output come from the comparisons-JSON file.
- Since 2026-10-28 the pipeline stops if the two disagree (see below). A simpler long-term design would be
  to derive the comparison's file from the sample map.

## Fixed

### Good-barcodes files from sample map and comparisons JSON could disagree (fixed 2026-10-28)
- Now checked in `barseq_qc_wf` (`check_good_barcodes_files_agree`, `modules/utils/main.nf`): for every sample of every
  comparison, the two files must be the same file or byte-identical copies; otherwise the run stops and lists
  each comparison/sample pair.
- Real case this would have caught: `results/2026_06_16_WITH_directionality` used a sample map pointing to the
  `*_WITH_directionality_in_locustag.csv` files, but `2026_04_01_comparisons.json` pointed every comparison to the
  `*_WITHOUT_directionality_in_locustag.csv` files. The merged matrices and all 810 RRA (Robust Rank Aggregation)
  result files therefore carry plain locus tags: that run contains no sense/antisense split.

### Samples missing from the sample map were dropped silently (fixed 2026-10-28)
- `reads_ch.join(sample_goodbarcodes_library_map)` drops samples without a map entry. The pipeline now logs a WARN
  listing them (`warn_samples_missing_from_sample_map`, `modules/utils/main.nf`). They are still dropped.
- Example: the `Liu_2021_QUICK` amylopectin comparison ran with 2 samples because IT043, IT032 and IT034 were not in the map.

### 32 GB memory requests for small python/R steps (fixed 2026-10-28)
- `python_basic` and `r_basic` requested 32 GB × attempt, `python_basic_quick_process` 16 GB × attempt; plot jobs sat
  PENDING ~2 h on a busy cluster. Peak RSS (resident memory) in the full Liu 2021 run (378 comparisons) was < 0.6 GB
  for every python/R step. Now 4 GB × attempt (retries: 4, 8, 12, 16 GB). mbarq (64 GB, peak 0.98 GB) is unchanged.

## Accepted / won't fix for now
- Container images `quay.io/jupyter/datascience-notebook` and `rocker/tidyverse:latest` are not pinned to a digest.
- Nextflow warns about undefined params `samplesheet` / `twofast2q_folder` (one of the two is always unset).

## Notes
- Directionality (sense/antisense) is supported by design, not by dedicated code: encode it in the `locus_tag`
  column of the good-barcodes file (e.g. `PV_ATCC8482_02977___sense`); the pipeline treats `locus_tag` as an opaque
  string (README, "TODO: You can add information to the `good_barcodes_with_locus_tag` file"). Both the sample map
  and the comparisons JSON must point to the directional files.
