# Output: results/Fig2D_kegg.csv; the dot-plot is rendered by Fig2_plot.py.

options(timeout = 300)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
})

.args <- commandArgs(trailingOnly = FALSE)
.fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
fig_dir <- normalizePath(file.path(dirname(.fa[1]), ".."))
proj_dir <- normalizePath(file.path(fig_dir, ".."))

input_file <- file.path(proj_dir, "figure1", "results", "all_diff_ac4c_sites.csv")
res_dir <- file.path(fig_dir, "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

peak_df <- read.csv(input_file, stringsAsFactors = FALSE)
genes_all <- unique(na.omit(peak_df$GeneName))
cat("Total differential genes:", length(genes_all), "\n")

eg <- bitr(genes_all, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
cat("Mapped Entrez IDs:", nrow(eg), "\n")

kegg_res <- enrichKEGG(
  gene          = eg$ENTREZID,
  organism      = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

if (!is.null(kegg_res) && nrow(kegg_res@result) > 0) {
  res <- kegg_res@result[kegg_res@result$p.adjust < 0.05, ]
  write.csv(res, file.path(res_dir, "Fig2D_kegg.csv"), row.names = FALSE)
} else {
  cat("No significant KEGG terms\n")
}

cat("Figure 2 KEGG enrichment complete\n")
