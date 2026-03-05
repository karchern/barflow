#!/usr/bin/env Rscript
library(tidyverse)

# parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: plot_volcano.r <mbarq_results_path> <output_path>")
}   

mbarq_results_path <- args[1]
output_path <- args[2]

# read mbarq results
mbarq_results <- read_tsv(mbarq_results_path)
print(mbarq_results)
