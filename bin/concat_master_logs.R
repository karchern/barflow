#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) {
    print(args)
  stop("Usage: concat_master_logs.R <file1> [<file2> ...] <output_file>")
}

out <- tail(args, 1)
files <- head(args, -1)

read_table_safe <- function(p) {
  if(!file.exists(p)) return(NULL)
  df <- read_tsv(p, col_types = cols(.default = 'c'))
  if(ncol(df) == 0) return(NULL)
  return(df)
}

dfs <- lapply(files, function(f){
  df <- read_table_safe(f)
  if(is.null(df)) return(NULL)
  name <- sub("\\.master_log\\.tsv$", "", basename(f))
  df <- df %>% mutate(comparison_name = name) %>% relocate(comparison_name)
  return(df)
})

dfs <- Filter(Negate(is.null), dfs)
if(length(dfs) == 0) {
  write_tsv(tibble(), out)
  quit(status = 0)
}

master <- bind_rows(dfs)
write_tsv(master, out)
