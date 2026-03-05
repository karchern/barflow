#!/usr/bin/env Rscript
library(tidyverse)

# TODO: Implement logging statement(s) to better understand what's going on in this script? Maybe even with a global log that could be written to a file per comparison.

do_low_barcode_count_filter <- function(merged_matrices, mbarq_meta, filter_on_what, lowly_abundant_barcode_cutoff) {
    # filter_on_what can be "controls" or "treatments" or "both"
    # lowly_abundant_barcode_cutoff is a non-negative integer number
    
    if (filter_on_what == "controls") {
        # filter based on control samples only
        samples_to_consider <- mbarq_meta %>%
            filter(treatment == "control") %>%
            pull(sampleID)
    } else if (filter_on_what == "treatments") {
        # filter based on treatment samples only
        samples_to_consider <- mbarq_meta %>%
            filter(treatment == 'treated') %>%
            pull(sampleID)
    } else if (filter_on_what == "both") {
        # filter based on both control and treatment samples
        samples_to_consider <- mbarq_meta %>%
            filter(treatment %in% c('control', 'treated')) %>%
            pull(sampleID)
    } else {
        stop("filter_on_what must be either 'controls' or 'treatments' or 'both'")
    }

    stopifnot(length(samples_to_consider) > 0)

    # Now, for each barcode, we want to check if the mean of counts across the samples_to_consider is >= lowly_abundant_barcode_cutoff
    filtered_merged_matrices <- merged_matrices %>%
        # use apply
        filter(
            apply(
                select(., all_of(samples_to_consider)),
                1,
                function(x) mean(x) >= lowly_abundant_barcode_cutoff
            )
        )
    
    return(filtered_merged_matrices)
}

# parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: filter_merged_barcode_matrix.r <merged_matrices_path> <mbarq_meta_path> <filter_on_what> <lowly_abundant_barcode_cutoff> <remove_all_0_barcodes> <output_path>")
}

merged_matrices_path <- args[1]
mbarq_meta_path <- args[2]
filter_on_what <- args[3]
lowly_abundant_barcode_cutoff <- as.numeric(args[4])
remove_all_0_barcodes <- as.logical(args[5])
output_path <- args[6]

print(
    paste0(
        "filter_on_what: ", filter_on_what, "\n",
        "lowly_abundant_barcode_cutoff: ", lowly_abundant_barcode_cutoff, "\n",
        "remove_all_0_barcodes: ", remove_all_0_barcodes, "\n"
    )
)

# Some basic checks on the input parameters:
## Ensure that lowly_abundant_barcode_cutoff is a non-negative integer number
stopifnot(is.numeric(lowly_abundant_barcode_cutoff))
stopifnot(lowly_abundant_barcode_cutoff >= 0)
stopifnot(lowly_abundant_barcode_cutoff == floor(lowly_abundant_barcode_cutoff))
## Ensure that remove_all_0_barcodes is a logical (might be false or FALSE or true or TRUE)
stopifnot(is.logical(remove_all_0_barcodes))

# load merged_matrices_path
merged_matrices <- read_csv(merged_matrices_path)

# load mbarq_meta_path
mbarq_meta <- read_csv(mbarq_meta_path)

if (remove_all_0_barcodes) {
  # remove barcodes that have 0 counts across all samples
  rem_row <- merged_matrices %>%
    select(-`#Feature`) %>%
    select(-locus_tag) %>%
    apply(
        ., 
        1,
        function(x) all(x == 0)
    )
    merged_matrices <- merged_matrices[!rem_row, ]
}

merged_matrices <- do_low_barcode_count_filter(
    merged_matrices = merged_matrices,
    mbarq_meta = mbarq_meta,
    filter_on_what = filter_on_what,
    lowly_abundant_barcode_cutoff = lowly_abundant_barcode_cutoff
)

write_csv(
    merged_matrices,
    output_path
)

