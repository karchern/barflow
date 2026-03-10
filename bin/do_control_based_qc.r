#!/usr/bin/env Rscript
library(tidyverse)

# Parse

compute_pre_mbarq_qc <- function(
    mbarq_results_path, mbarq_meta_path
) {

data <- read_csv(mbarq_results_path) %>%
  pivot_longer(
    cols = -c(`#Feature`, locus_tag),
    names_to = "sampleID",
    values_to = "barcode_count"
  ) %>%
  left_join(
    read_csv(mbarq_meta_path)
  )

#Total number of barcodes detected, as descending barplot colored by type of sample
percent_good_barcodes_detected <- data %>%
  group_by(sampleID, treatment) %>%
  summarise(
    fraction_total_barcodes_detected = mean(barcode_count >= 2)
    ) %>%
  ungroup() %>%
  arrange(fraction_total_barcodes_detected)

return(percent_good_barcodes_detected)
}





args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
    stop("Usage: do_control_based_qc.r <mbarq_results_path> <mbarq_meta_path> <comparison_name> <output_path>")
}   

mbarq_results_path <- args[1]
mbarq_meta_path <- args[2]
comparison_name <- args[3]
output_path <- args[4]

results <- compute_pre_mbarq_qc(
    mbarq_results_path,
    mbarq_meta_path
)

write_tsv(
    results,
    output_path
)







