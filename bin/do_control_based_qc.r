#!/usr/bin/env Rscript
library(tidyverse)

# Parse

compute_pre_mbarq_qc <- function(
    mbarq_results_path, mbarq_meta_path
) {

data <- read_csv(mbarq_results_path) %>%
  pivot_longer(
    cols = -any_of(c("#Feature", "locus_tag", "contig", "position")),
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
    stop("Usage: do_control_based_qc.r <mbarq_results_path> <mbarq_meta_path> <output_path> <output_path_aggregated>")
}   

mbarq_results_path <- args[1]
mbarq_meta_path <- args[2]
output_path <- args[3]
output_path_aggregated <- args[4]

results <- compute_pre_mbarq_qc(
    mbarq_results_path,
    mbarq_meta_path
)

results_aggregated <- results %>%
  summarise(
    num_treatment_samples = sum(treatment == 'treated'),
    num_control_samples = sum(treatment == 'control'),    
    mean_fraction_total_barcodes_detected_in_controls = mean(fraction_total_barcodes_detected[treatment == 'control']),
    mean_fraction_total_barcodes_detected_in_treatments = mean(fraction_total_barcodes_detected[treatment == 'treated'])
  )

if (results_aggregated$num_treatment_samples == 0) {
  exit(
    "No treatment samples found in the data. Cannot compute QC metrics. Also this should not ever happen."
  )
}

if (results_aggregated$num_control_samples == 0) {
  exit(
    "No control samples found in the data. Cannot compute QC metrics. Also this should not ever happen."
  )
}

write_tsv(
    results,
    output_path
)

write_tsv(
    results_aggregated,
    output_path_aggregated
)







