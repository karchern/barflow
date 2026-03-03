# barflow - Nextflow pipeline for TnSeq analysis

## Quick Start

with Nextflow installed, run `bash test_run.sh` to run pipeline with dummy data. 

## Input Logic:

This pipeline assumes you have already characterized your input library and have a list of reliable barcodes (typically those of sufficient abundance + evidence of single integration in initial library).

Specifically, as input, we need:
- Per-sample single-end fastq files corresponding to invididual samples/replicates that were bar-. Specified via parameter `--samplesheet` (see below).
- File specifying reliable barcodes. Specified via parameter `--good_barcodes_csv` (see below)
- File specifying which samples are supposed to be treated as 'control' or 'treatment' for the differential barcode analysis via parameter `--comparisons` (see below).

## Input file Specifics

### Per-sample fastq files

Each sample is assumed to be a single-end fastq file. Add the complete path of each file to a textfile and supply it the path to this file using the `--samplesheet` parameter. 
The file in each row will be sent to `2fast2q` to have the barcodes extracted. 

A dummy samplesheet looks like this:

```
/g/typas/Personal_Folders/Nic/barflow/test_data/test_reads/Sample_treatment_1.fastq.gz
/g/typas/Personal_Folders/Nic/barflow/test_data/test_reads/Sample_treatment_2.fastq.gz
/g/typas/Personal_Folders/Nic/barflow/test_data/test_reads/Sample_control_1.fastq.gz
/g/typas/Personal_Folders/Nic/barflow/test_data/test_reads/Sample_control_2.fastq.gz
```

Fastq files can be gzipped or uncompressed, and sample-IDs are dynamically extracte from the full path. The sample-IDs extracted from the paths above would be as below. Make sure your file nanmes do not contain dots, as those mess things up.

```
Sample_treatment_1
Sample_treatment_2
Sample_control_1
Sample_control_2
```

Make sure to correctly set the upstream- and downstream nucleotide sequences, either by changing them in `nextflow.config` or, alternatively, by
setting them on te command line via parameters `--upstream_seq` and `--downstream_seq`.

### Reliable barcodes

A barcoded, random transposon mutagenesis experiment generates a large number of mutants. For several reasons, not all of these mutants are typically used in downstream experiments.
This pipeline assumes that a set of good/reliable barcodes has been determined from the initial sequencing run of the library. 

Specifically, the input file (supplied via `--good_barcodes_csv`) is expected to be a 2-column csv file without column names where
- each row corresponds to a barcode
- column 1: barcode sequence
- column 2: locus tag information (so, typically, on which gene or which intergenic region the transposon with this barcode is inserted)

The information in the 2nd column will be used for (aggregation of) differential barcode abundance results down the line. 

A dummy file with barcode information would look like this

```
TGGTGTTGCGCGGCGACTGCCCACG,PV_ATCC8482_03882
CGGTCTGAGTCTACGAACGCGTACA,PV_ATCC8482_00424
AACTCCGTATAAATTGCCGTATAGA,PV_ATCC8482_01144
TCTGACGTATCAAAGGCGCTGCGGT,PV_ATCC8482_03772
ACGTGCGTGGTAAGGGTGGGACGCA,IG-between-PV_ATCC8482_00619-and-PV_ATCC8482_00620
AGTGGGTTGGATCTCTGAGGTTTAC,PV_ATCC8482_00281
AATATTCCGCTTTTGTATTGAAGAG,PV_ATCC8482_03009
CGGATTATGTCAATACAGGGGTATT,PV_ATCC8482_01688
TAGTTATTATCTACTCGGCTTATTA,PV_ATCC8482_02738
TGGGAACTCTTAGAGGCTTGCGGAC,PV_ATCC8482_00692
AGTTCGGATATTAGCATTGACTTGT,PV_ATCC8482_03676
```


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

# Quirks

- The comparisons defined in the `comparison.json` are based on the entries of the `samplesheet.txt`. What this means is that only what is in the `samplesheet.txt` can be handed over to the differential abundance computation. It is therefore crucial that you 
  - If you use the glob syntax, it will base this on the sampleIDs inferred from the `samplesheet.txt`, _not_ what has already been processed by `2fast2q`
  - If you use the explicit syntax and you specify sampleIDs that are _not_ in the `samplesheet.txt`, the pipeline will crash at the mbarq step because it is looking for a sample to load that it cannot find.
- TODO: Explain negative selection

This is due to the way Nextflow handles inputs