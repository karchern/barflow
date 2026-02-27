barflow - Nextflow pipeline for TnSeq analysis

Structure:
- main.nf - DSL2 workflow
- modules/create_barcode_count_matrix - module that calls 2FAST2Q to extract barcodes and produce per-sample barcode counts
- modules/merge_and_mbarq - module that merges per-sample barcode counts, adds locus tag information to the merged barcode count matrix and runs mbarq for differential barcode/locus tag analysis

simple/test case:
with Nextflow installed, run `bash test_run.sh` 

input logic:
This pipeline assumes you have already characterized your input library and have a list of reliable barcodes (what we call `good barcodes` below) as well as their link to locus tags (genes or inter-genic regions).
As input, we take per-sample single-end fastq files corresponding to invididual samples/replicates that were bar-seqed.

Output:
- 2fast2q
- comparisons
- mbarq
  

Notes:
- Barcode length is parameterised via --barcode_length (default 8)
- Inputs should be single-end FASTQ files. The pipeline derives sample IDs from filenames by removing common fastq extensions.
- The merge step produces `results/barcode_count_matrix.tsv` with barcodes as rows and samples as columns.

Next steps you mentioned:
- Provide explicit barcode extraction parameters (length, mismatch allowance)
- Define comparison config format for downstream branches
- Add tests and example data
