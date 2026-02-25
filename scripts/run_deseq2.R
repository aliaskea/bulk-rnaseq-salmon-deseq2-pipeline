suppressPackageStartupMessages({
  library(optparse)
  library(tximport)
  library(DESeq2)
})

option_list <- list(
  make_option("--samples", type = "character"),
  make_option("--contrasts", type = "character"),
  make_option("--tx2gene", type = "character"),
  make_option("--salmon_dir", type = "character"),
  make_option("--out", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

samples <- read.delim(opt$samples, stringsAsFactors = FALSE)
contrasts <- read.delim(opt$contrasts, stringsAsFactors = FALSE)
tx2gene <- read.delim(opt$tx2gene, stringsAsFactors = FALSE)

files <- file.path(opt$salmon_dir, samples$sample, "quant.sf")
names(files) <- samples$sample

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)

dds <- DESeqDataSetFromTximport(
  txi = txi,
  colData = samples,
  design = ~ condition
)
dds <- DESeq(dds)

contrast_row <- contrasts[1, ]
res <- results(dds, contrast = c("condition", contrast_row$group1, contrast_row$group2))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
write.table(res_df, file = opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
