barflow - barebones Nextflow scaffold for transposon barcode sequencing

Structure:
- main.nf - DSL2 workflow
- modules/fastp - fastp module (QC + adapter removal)
- modules/create_barcode_count_matrix - module that calls 2FAST2Q to extract barcodes and produce per-sample counts
- envs/ - conda environment YAMLs for tools
- nextflow.config - profiles for local and slurm

Usage examples:
  nextflow run main.nf -profile standard --reads './data/*.{fastq,fq,fastq.gz,fq.gz}' --outdir results
  nextflow run main.nf -profile slurm --reads '/path/to/*.fastq.gz' --outdir results

Notes:
- Barcode length is parameterised via --barcode_length (default 8)
- Inputs should be single-end FASTQ files. The pipeline derives sample IDs from filenames by removing common fastq extensions.
- The merge step produces `results/barcode_count_matrix.tsv` with barcodes as rows and samples as columns.

Next steps you mentioned:
- Provide explicit barcode extraction parameters (length, mismatch allowance)
- Define comparison config format for downstream branches
- Add tests and example data
