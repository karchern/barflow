#!/usr/bin/env Rscript
library(tidyverse)

# parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("Usage: plot_volcano.r <mbarq_results_path> <comparison_name> <output_path>")
}   

mbarq_results_path <- args[1]
comparison_name <- args[2]
output_path <- args[3]

wrap_len <- 30
comparison_name_wrapped <- gsub(paste0("(.{", wrap_len, "})"), "\\1\n", comparison_name)
comparison_name_wrapped <- sub("\n$", "", comparison_name_wrapped)

print(
    paste0(
        "mbarq_results_path: ", mbarq_results_path, "\n",
        "comparison_name: ", comparison_name, "\n",
        "output_path: ", output_path, "\n"
    )
)

# read mbarq results
all_data <- read_csv(mbarq_results_path)

all_data <- all_data %>%
  mutate(
    # Combine the two FDR columns (take the smaller of neg/pos)
    # TODO: Discuss with Carlos if this shouldn't rather be ifelse(LFC > 0, pos_selection_fdr, neg_selection_fdr)
    FDR = pmin(neg_selection_fdr, pos_selection_fdr, na.rm = TRUE),
    # Mark significance
    sig = if_else(abs(LFC) > 1.5 & FDR < 0.05, "Significant", "Not significant"),
  )

total_number_of_barcodes <- sum(all_data$number_of_barcodes)
total_number_of_genes <- nrow(all_data)

p <- ggplot(all_data %>% 
         filter(number_of_barcodes > 1),
         aes(x = LFC, y = -log10(FDR))) +
    geom_point(aes(color = sig), alpha = 0.3, size = 1.5) +
    scale_color_manual(values = c("Significant" = "red", "Not significant" = "grey70")) +
    geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black") +
    labs(
      x = "Log2 Fold Change",
      y = "-log10(FDR)",
      color = ""
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "top"
    ) +
#   scale_x_continuous(breaks = seq(-16, 16, by = 4),
#                      limits = c(-20, 20)) +
    ggtitle(
        paste0("Volcano Plot for ", comparison_name_wrapped)
    )

ggsave(
    filename = output_path,
    plot = p,
    width = 4,
    height = 4
  )
