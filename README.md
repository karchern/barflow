# barflow - Nextflow pipeline for TnSeq analysis

## Quick Start

with Nextflow installed, run `bash test_run.sh` to run pipeline with dummy data. 

## Input Logic:

This pipeline assumes you have already characterized your input library and have a list of reliable barcodes (typically those of sufficient abundance + evidence of single integration in initial library).

Specifically, as input, we need:
- Per-sample single-end fastq files corresponding to invididual samples/replicates that were bar-. Specified via parameter `--samplesheet` (see below).
- If you do this, you also need to supply `--sample_goodbarcodes_map`: a 2-column csv with the first column containing the sample_id (matching what is in samplesheet) and the second column the path to the file containing barcodes and their genomic locus tag (see below)
  - Alternatively, you can supply a folder containing 2fast2q output files, one per sample. These can be precomputed using this pipeline (in the first step, this pipeline runs 2fast2q on all samples supplied in `--samplesheet`), or elsewhere. Files should be called `${sampleID}.2fast2q`, and path to the folder containing these files should be specified using `--2fast2q_folder`. This is mutually exclusive with `--samplesheet`, and an error will be thrown if both are supplied at the same time.
- File specifying which samples are supposed to be treated as 'control' or 'treatment' for the differential barcode analysis as well as the barcodes to be analyzed as well as their locus tag (same files as supplied to `--sample_goodbarcodes_map`) via parameter `--comparisons` (see below).

## Input file Specifics

### Per-sample fastq files

Each sample is assumed to be a single-end fastq file. Add the complete path of each file to a textfile and supply it the **absolute** path to this file using the `--samplesheet` parameter. 
The file in each row will be sent to `2fast2q` to have their barcodes extracted. 

A dummy samplesheet looks like this:

```
some/complete/path/to/test_data/test_reads/Sample_treatment_1.fastq.gz
some/complete/path/to/test_data/test_reads/Sample_treatment_2.fastq.gz
some/complete/path/to/test_data/test_reads/Sample_control_1.fastq.gz
some/complete/path/to/test_data/test_reads/Sample_control_2.fastq.gz
```

Fastq files can be gzipped or uncompressed, and sample-IDs are dynamically extracte from the full path. The sample_IDs extracted from the paths above would be as below. **Make sure your file names do not contain dots**, as those will mess up the extraction of the sample_IDs.

```
Sample_treatment_1
Sample_treatment_2
Sample_control_1
Sample_control_2
```

Make sure to correctly set the upstream- and downstream nucleotide sequences, either by changing them in `nextflow.config` or, alternatively, by
setting them on te command line via parameters `--upstream_seq` and `--downstream_seq`.

### Map between samples and their libraries

Each sample is based on a specific barcoded transposon library, and each of those has a set of 'good' barcodes that can be defined in various ways. If one wants to run 2fast2q, one needs to also supply a set of barcodes for 2fast2q to look for. To that end, supply a 2-column csv with the first column containing the sample_id (matching what is in samplesheet) and the second column the path to the file containing barcodes and their genomic locus tag (see below). See section `Reliable barcodes`.

For example:

```
Sample_treatment_1,some/complete/path/to/test_data/good_barcodes_with_locustag.csv
Sample_treatment_2,some/complete/path/to/test_data/good_barcodes_with_locustag.csv
Sample_control_1,some/complete/path/to/test_data/good_barcodes_with_locustag.csv
Sample_control_2,some/complete/path/to/test_data/good_barcodes_with_locustag.csv
```

### Alternative to Per-sample fastq files - Preprocessed 2fast2q output

It might not always be desired or efficient to rerun 2fast2q. For example, one might want to explore different parameter settings and hence want to run mbarq over a variety of parameters or comparisons. To this end, you can - as an alternative to supplying `--samplesheet`, supply a path to a folder containing precomputed 2fast2q output files using `--2fast2q_folder`. Each file in there should be the raw output of 2fast2q, and named as ${sample_id}.2fast2q. For example, using one of the sample_IDs from the dummy samplesheet directly above, a potential output file could be `Sample_treatment_1.2fast2q`.

### Comparison JSON

After extracting barcodes (and counting them), the next step of a barcoded transposon experiment is typically to compare barcode (i.e. strain) abundances between 'treatment' and 'control' conditions, in order to find differentially abundant barcodes, and hence draw conclusions about strain fitness and ultimately the potential function/role of genetic elements that were altered in those strains (by an inserting transposon).

In order to control what sample should be compared against what, design a `json` file as outlined below and supply it via `--comparisons`.

The file needs to contain a global field 'comparisons', in which each comparison is nested. Each comparison, in turn, _needs_ to contain the following fields:
- name: A string specifying the name of the comparison. For example: example, "B_uniformis__mgam_against_B_uniformis__VB"
- treatments: A list of strings, each corresponding to sample_IDs that should be used as treatments. Those need to be either contained in the samplesheet (in case `--samplesheet` is supplied) or alternatively need to be presnt as precomputed 2fast2q outputs in the path supplied by `--2fast2q_folder`
- controls: A list of strings, each corresponding to sample_IDs that should be used as controls.
- good_barcodes_file: A string specifying the path to the file containing good barcodes and their locus tag (see section `Reliable barcodes`)

For convenience, the comparison json also supports glob syntax, where asterisks are expanded based on whatever can be found in `2fast2q` output folder.

Optional fields are:
- treatments_negative_selection: A list of strings. Sample_IDs in these lists are ignored.
- controls_negative_selection: A list of strings. Sample_IDs in these lists are ignored.

A dummy comparison json looks like this:

```
{
  "comparisons": [
    {
      "name": "explicit_syntax",
      "treatments": [
        "Sample_treatment_1",
        "Sample_treatment_2"
    ],
      "controls": [
        "Sample_control_1", 
        "Sample_control_2"
    ],
     "good_barcodes_file" : 
     "some/path/to/test_data/good_barcodes_with_locustag.csv"
    },
    {
      "name": "glob_syntax",
      "treatments": [
        "Sample_treatment_*"
    ],
      "controls": [
        "Sample_control_*"
    ],
     "good_barcodes_file" : 
     "some/path/to//test_data/good_barcodes_with_locustag.csv"
    }
    ,
    {
      "name": "glob_syntax_with_negative_selection",
      "treatments": [
        "Sample_treatment_*"
    ],
      "controls": [
        "Sample_control_*"
    ],
     "treatments_negative_selection" :
     [
        "Sample_treatment_1"
     ],
     "controls_negative_selection" :
     [
        "Sample_control_2"
     ],
     "good_barcodes_file" : 
     "some/path/to/barflow/test_data/good_barcodes_with_locustag.csv"
    }
  ]
}
```

Samples can be treatment/controls of as many comparisons as you want, you don't need to be worried about them being unique.

### Reliable barcodes

A barcoded, random transposon mutagenesis experiment generates a large number of mutants. For several reasons, not all of these mutants are typically used in downstream analysis.
This pipeline assumes that a set of good/reliable barcodes has been determined from the initial sequencing run of the library, for each of the libraries used. 

Specifically, the input file (supplied via as a mandatory field `good_barcodes_file` in the comparson `json`) is expected to be a 2-column csv file without column names where each row corresponds to a barcode
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
- Do not use relative paths, only complete paths in `samplesheet` and `sample_goodbarcodes_map`.
- TODO: Explain negative selection
- TODO: You can add information to the `good_barcodes_with_locus_tag` file. For example, if you're interested in directionality biases, you can simply add this information to the `locus_tag` column. Importantly, do not add a 3rd column, but instead modify the 2nd to contain this information.

This is due to the way Nextflow handles inputs