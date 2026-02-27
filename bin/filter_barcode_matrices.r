#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: filter_barcode_matrices.r <input_matrix> <output_matrix> <good_barcodes_file> <barcode_information_file>")
}

args <- commandArgs(trailingOnly = TRUE)
input_matrix <- args[1]
output_matrix   <- args[2]


