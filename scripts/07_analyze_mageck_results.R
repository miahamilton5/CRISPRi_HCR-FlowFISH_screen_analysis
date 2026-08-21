#!/usr/bin/env Rscript
# Call significant elements from a MAGeCK RRA gene_summary.txt and plot a volcano.
#
# For each element, takes the more significant of the positive/negative
# selection p-values and log2 fold-change, applies Bonferroni p-value
# correction, and classifies elements as significantly enriched (repression
# increases the sorted signal), significantly depleted (repression decreases
# it), or not significant.
#
# Usage:
#   Rscript 07_analyze_mageck_results.R <gene_summary.txt> <output_table.csv>

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))

args <- commandArgs(trailingOnly = TRUE)
gene_summary_path <- args[1]
out_path <- args[2]

results <- fread(gene_summary_path) %>%
  select(id, `pos|p-value`, `neg|p-value`, `pos|lfc`, `pos|fdr`, `neg|fdr`, `neg|rank`, `pos|rank`) %>%
  mutate(
    pval = pmin(`pos|p-value`, `neg|p-value`),
    fdr = pmin(`pos|fdr`, `neg|fdr`),
    # negate MAGeCK's pos|lfc (enrichment in the bottom/treatment bin) so that
    # sign matches biological direction: negative = repression decreases
    # expression, positive = repression increases expression
    l2fc = -`pos|lfc`,
    pval_adj = p.adjust(pval, method = "bonferroni"),
    significant = case_when(
      pval_adj <= 0.05 & l2fc > 0 ~ "significant_positive_selection",
      pval_adj <= 0.05 & l2fc < 0 ~ "significant_negative_selection",
      TRUE ~ "not_significant"
    ),
    label = case_when(
      significant != "not_significant" & (`neg|rank` < 3 | `pos|rank` < 3) ~ id
    )
  )

fwrite(results, out_path)

cat("Significant elements:\n")
print(results %>% count(significant))

ggplot(results, aes(x = l2fc, y = -log10(pval_adj), color = significant)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(
    aes(label = label),
    size = 3,
    max.overlaps = Inf,
    box.padding = 1,
    point.padding = 0.5,
    min.segment.length = 0,
    segment.size = 0.3,
    force = 5,
    seed = 1,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c(
    "not_significant" = "lightgray",
    "significant_negative_selection" = "#EC5F67",
    "significant_positive_selection" = "#6699CC"
  )) +
  scale_x_continuous(limits = c(NA, 2)) +
  theme_classic() +
  labs(x = "log2(fold-change)", y = "-log10(adjusted p-value)")
