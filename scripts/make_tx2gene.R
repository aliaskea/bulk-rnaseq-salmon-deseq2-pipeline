args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript scripts/make_tx2gene.R <annotation.gtf> <out.tsv>")
}

gtf_file <- args[[1]]
out_file <- args[[2]]

suppressPackageStartupMessages({
  library(rtracklayer)
})

gtf <- import(gtf_file)
gtf_df <- as.data.frame(gtf)

if (!all(c("type", "transcript_id", "gene_id") %in% colnames(gtf_df))) {
  stop("GTF is missing required fields: type, transcript_id, gene_id")
}

tx2gene <- unique(
  gtf_df[gtf_df$type == "transcript", c("transcript_id", "gene_id")]
)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write.table(
  tx2gene,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
