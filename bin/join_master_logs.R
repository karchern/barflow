#!/usr/bin/env Rscript

# Join up to four tabular log files by their first column (key) and write to output

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3) {
  stop("Usage: join_master_logs.R comparison pre_mbarq_log post_mbarq_log output_file")
}

comparison <- args[1]
pre_m <- args[2]
post_m <- args[3]
out <- args[4]

read_table_safe <- function(p) {
  if(!file.exists(p)) return(tibble())
  df <- read_tsv(p, col_types = cols(.default = 'c'))
  if(ncol(df) == 0) return(tibble())
  #names(df)[1] <- 'key'
  names(df)[1] <- "key"
  return(df)
}

files <- list(pre_m, post_m)
dfs <- lapply(files, read_table_safe)
# drop empty tibbles
dfs <- dfs[sapply(dfs, ncol) > 0]

if(length(dfs) == 0) {
  write_tsv(tibble(), out)
  quit(status = 0)
}

print(dfs)

master <- reduce(dfs, ~ full_join(.x, .y, by = 'key', suffix = c(".pre_mbarq", ".post_mbarq")))
colnames(master)[colnames(master) == "key"] <- "comparison"
master$comparison[1] <- comparison
stopifnot(nrow(master) == 1)
write_tsv(master, out)
