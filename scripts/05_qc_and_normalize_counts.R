#!/usr/bin/env Rscript
# QC and RPM-normalize a raw gRNA count table before MAGeCK.
#
# Steps:
#   1. Plot the raw bulk (unsorted) count distribution per guide, marking the
#      minimum-count filtering threshold.
#   2. Filter out guides with fewer than 10 raw reads in the bulk sample.
#   3. Convert to reads-per-million (RPM) within each sample.
#   4. Plot a sample-sample correlation matrix as a replicate-concordance QC check.
#
# Usage:
#   Rscript 05_qc_and_normalize_counts.R <raw_count_table.csv> <output_rpm_table.csv> \
#     <bulk_column_name> <output_correlation_plot.png> <output_bulk_diversity_plot.png>
#
# <raw_count_table.csv> must have a "grna" column plus one column per sample
# (as produced by combining per-sample .counts files, e.g. via cbind on the
# shared "grna" column).

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(corrplot))

args <- commandArgs(trailingOnly = TRUE)
raw_counts_path <- args[1]
out_path <- args[2]
bulk_col <- args[3]
correlation_plot_path <- args[4]
diversity_plot_path <- args[5]

MIN_BULK_COUNT <- 10

raw_counts <- fread(raw_counts_path)

# 1. Library diversity in the bulk sample, before filtering, with the
#    filtering threshold marked. A +1 pseudocount avoids -Inf on the log10
#    axis for guides with zero raw reads.
diversity_plot <- ggplot(raw_counts, aes(x = .data[[bulk_col]] + 1)) +
  geom_histogram(bins = 60, fill = "#6699CC", color = NA) +
  geom_vline(xintercept = MIN_BULK_COUNT, linetype = "dashed", color = "#EC5F67") +
  scale_x_log10() +
  theme_classic() +
  labs(
    x = paste0(bulk_col, " raw read count + 1 (log10)"),
    y = "Number of guides",
    title = "Bulk sample library representation"
  )
ggsave(diversity_plot_path, diversity_plot, width = 6, height = 4.5, dpi = 200)

# 2. Filter out guides with low bulk representation
raw_counts <- raw_counts[raw_counts[[bulk_col]] >= MIN_BULK_COUNT, ]

# 3. RPM-normalize each sample column
sample_cols <- setdiff(colnames(raw_counts), c("grna", "id"))
total_counts <- colSums(raw_counts[, ..sample_cols], na.rm = TRUE)

rpm <- as.data.table(raw_counts)
for (col in sample_cols) {
  rpm[[col]] <- (raw_counts[[col]] / total_counts[[col]]) * 1e6
}

fwrite(rpm, out_path)

# 4. Replicate-concordance QC: sample-sample correlation heatmap
png(correlation_plot_path, width = 1400, height = 1400, res = 200)
cor(rpm[, ..sample_cols], method = "pearson", use = "complete.obs") %>%
  corrplot(order = "hclust", tl.col = "black", addCoef.col = "black", number.cex = 0.5)
dev.off()
