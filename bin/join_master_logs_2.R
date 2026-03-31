#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 4) {
  stop("Usage: join_master_logs.R pre_barseq_qc post_barseq_qc mbarq_qc_logs out")
}

pre_barseq_qc <- args[1]
post_barseq_qc <- args[2]
mbarq_qc_logs <- args[3]
out <- args[4]

files <- list(pre_barseq_qc, post_barseq_qc, mbarq_qc_logs)
dfs <- lapply(files, read_tsv)
# drop empty tibbles
dfs <- dfs[sapply(dfs, ncol) > 0]

if(length(dfs) == 0) {
  write_tsv(tibble(), out)
  quit(status = 0)
}

print(dfs)

master <- full_join(
    dfs[[1]],
    dfs[[2]],
    by = 'comparison', suffix = c(".pre_barseq_qc", ".post_barseq_qc")
) %>%
    full_join(
        dfs[[3]],
        by = 'comparison'
    )

print(master)

write_tsv(master, out)
