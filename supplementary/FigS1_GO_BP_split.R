# Supplementary Figure S1 GO BP enrichment, up- and down-acetylated genes.

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
})

.args <- commandArgs(trailingOnly = FALSE)
.fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
here <- normalizePath(dirname(.fa[1]))
proj_dir <- normalizePath(file.path(here, ".."))

input_file <- file.path(proj_dir, "figure1", "results", "all_diff_ac4c_sites.csv")
res_dir <- file.path(proj_dir, "supplementary", "results")

dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

peak_df <- read.csv(input_file, stringsAsFactors = FALSE)

genes_up   <- unique(na.omit(peak_df$GeneName[peak_df$direction == "up"]))
genes_down <- unique(na.omit(peak_df$GeneName[peak_df$direction == "down"]))

cat("Up-regulated genes:", length(genes_up), "\n")
cat("Down-regulated genes:", length(genes_down), "\n")

run_go <- function(genes, csv_name) {
  if (length(genes) < 5) {
    cat(csv_name, ": gene count < 5, skipped\n")
    return(invisible(NULL))
  }

  eg <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  cat("Mapped Entrez IDs:", nrow(eg), "\n")

  go <- enrichGO(
    gene          = eg$ENTREZID,
    OrgDb         = org.Mm.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )

  if (is.null(go) || nrow(go@result) == 0) {
    cat(csv_name, ": no significant GO terms\n")
    return(invisible(NULL))
  }

  res <- go@result[go@result$p.adjust < 0.05, ]
  write.csv(res, file.path(res_dir, csv_name), row.names = FALSE)

   invisible(go)
}

cat("\n--- Supplementary Fig. S1A (up-regulated) ---\n")
run_go(genes_up, "SourceData_FigS1A_GO_BP_up.csv")

cat("\n--- Supplementary Fig. S1B (down-regulated) ---\n")
run_go(genes_down, "SourceData_FigS1B_GO_BP_down.csv")

cat("\nFigure S1 GO BP enrichment complete\n")
