#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: join_mbarq_matrices.R <treatments.tsv> <controls.tsv> <output.tsv>")
}

# Usage
# Rscript join_mbarq_matrices.R treatments.tsv controls.tsv ${comparison_name}.merged.barcode.matrices.tsv

args <- commandArgs(trailingOnly = TRUE)
treatments_file <- args[1]
controls_file   <- args[2]
good_barcodes_csv_path <- args[3]
output_file     <- args[4]

# treatments_file is a 2-column tsv without headers, first co sample-ID, second col path to matrix
# The file contains many rows, one for each sample in the treatment group
# the matrix itself is a csv and looks like this: 1-st row is header line and starts with hasthtag: 
# Feature and sample name column (mappable to sample-ID in the first col of treatments_file) and then barcode columns
# the rest of the rows are features and their counts in each barcode, e.g.:
### #Feature, Sample_treatment_2

# controls_file is the same format as treatments_file, but for control samples
# output_file is the name of the output file to write

# The function below reads treatments_file and controls_file, left_joins them (with reduce-type behaviour)
# using the sample-ID as key, and writes the output file to disk

# Read treatments and controls into data frames

# Usage:
# Rscript join_mbarq_matrices.R treatments.tsv controls.tsv ${comparison_name}.merged.barcode.matrices.tsv

read_samples_list <- function(path) {
  if (!file.exists(path)) stop(sprintf("Samples list not found: %s", path))
  df <- readr::read_tsv(path, col_names = FALSE, show_col_types = FALSE)
  if (ncol(df) != 2) stop(sprintf("Samples list %s must have exactly 2 columns: sample_id <tab> path_to_matrix", path))
  tibble(sample = as.character(df[[1]]), file = as.character(df[[2]]))
}

samples_treat <- read_samples_list(treatments_file) %>%
    mutate(
        data = map(file, ~ {
            if (!file.exists(.x)) stop(sprintf("Matrix file not found: %s", .x))
            readr::read_csv(.x, show_col_types = FALSE)
        })
    )
samples_ctrl  <- read_samples_list(controls_file) %>%
    mutate(
        data = map(file, ~ {
            if (!file.exists(.x)) stop(sprintf("Matrix file not found: %s", .x))
            readr::read_csv(.x, show_col_types = FALSE)
        })
    )
all_samples <- bind_rows(samples_treat, samples_ctrl)

message(sprintf("Found %d samples (treatments=%d controls=%d)", nrow(all_samples), nrow(samples_treat), nrow(samples_ctrl)))

all_data_reduced <- 
    reduce(
        all_samples$data,
        full_join,
        by = c("#Feature")
    )

good_barcodes_with_locustag <- read_csv(good_barcodes_csv_path, col_names = F)
colnames(good_barcodes_with_locustag) <- c("#Feature", "locus_tag")

all_data_reduced <- all_data_reduced %>%
    left_join(good_barcodes_with_locustag, by = "#Feature") %>%
    relocate(`#Feature`, locus_tag)

write_csv(all_data_reduced, output_file)