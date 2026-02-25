#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/summarize_deseq2.R",
      "<samples.tsv> <contrasts.tsv> <tx2gene.tsv> <salmon_dir> <deseq2_results.tsv> <out_dir>"
    )
  )
}

samples_file <- args[[1]]
contrasts_file <- args[[2]]
tx2gene_file <- args[[3]]
salmon_dir <- args[[4]]
results_file <- args[[5]]
out_dir <- args[[6]]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- read.delim(samples_file, stringsAsFactors = FALSE)
contrasts <- read.delim(contrasts_file, stringsAsFactors = FALSE)
tx2gene <- read.delim(tx2gene_file, stringsAsFactors = FALSE)
res <- read.delim(results_file, stringsAsFactors = FALSE)

# ---- Basic summary statistics ----
res_valid <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
sig_005 <- res_valid[res_valid$padj < 0.05, ]
sig_010 <- res_valid[res_valid$padj < 0.10, ]
sig_005_lfc1 <- sig_005[abs(sig_005$log2FoldChange) >= 1, ]

summary_df <- data.frame(
  metric = c(
    "total_rows",
    "tested_rows_non_na",
    "padj_lt_0_10",
    "padj_lt_0_05",
    "padj_lt_0_05_abs_log2fc_ge_1",
    "up_padj_lt_0_05",
    "down_padj_lt_0_05"
  ),
  value = c(
    nrow(res),
    nrow(res_valid),
    nrow(sig_010),
    nrow(sig_005),
    nrow(sig_005_lfc1),
    sum(sig_005$log2FoldChange > 0),
    sum(sig_005$log2FoldChange < 0)
  )
)

write.table(
  summary_df,
  file = file.path(out_dir, "summary_stats.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ---- Top hits table ----
up_hits <- sig_005[sig_005$log2FoldChange > 0, ]
down_hits <- sig_005[sig_005$log2FoldChange < 0, ]

up_hits <- up_hits[
  order(up_hits$padj),
  c("gene_id", "log2FoldChange", "padj", "baseMean", "pvalue"),
  drop = FALSE
]
down_hits <- down_hits[
  order(down_hits$padj),
  c("gene_id", "log2FoldChange", "padj", "baseMean", "pvalue"),
  drop = FALSE
]

top_up <- head(up_hits, 20)
top_down <- head(down_hits, 20)
top_up$direction <- rep("up_in_group1", nrow(top_up))
top_down$direction <- rep("down_in_group1", nrow(top_down))
top_hits <- rbind(top_up, top_down)

write.table(
  top_hits,
  file = file.path(out_dir, "top_hits.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ---- Volcano plot ----
plot_df <- res_valid
plot_df$neglog10_padj <- -log10(plot_df$padj)
plot_df$significant <- ifelse(plot_df$padj < 0.05 & abs(plot_df$log2FoldChange) >= 1, "yes", "no")

p_volcano <- ggplot(plot_df, aes(x = log2FoldChange, y = neglog10_padj, color = significant)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c("no" = "grey70", "yes" = "#d62728")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(
    title = "Volcano plot",
    subtitle = paste0(contrasts$group1[1], " vs ", contrasts$group2[1]),
    x = "log2 fold change",
    y = "-log10 adjusted p-value",
    color = "padj<0.05\n|log2FC|>=1"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(out_dir, "volcano.png"),
  plot = p_volcano,
  width = 7,
  height = 5,
  dpi = 150
)

# ---- PCA plot (rebuild dds from Salmon outputs) ----
files <- file.path(salmon_dir, samples$sample, "quant.sf")
names(files) <- samples$sample

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)

samples$condition <- as.factor(samples$condition)
dds <- DESeqDataSetFromTximport(txi = txi, colData = samples, design = ~ condition)
vsd <- vst(dds, blind = TRUE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, name = rownames(pca_data))) +
  geom_point(size = 3) +
  labs(
    title = "PCA of samples (VST)",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(out_dir, "pca.png"),
  plot = p_pca,
  width = 7,
  height = 5,
  dpi = 150
)

