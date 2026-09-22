# Contrast 1: in-house mRNA-seq, IR (n=3) vs Sham (n=3)
# Contrast 2: GSE322673 cold-exposure mRNA-seq, Cold (n=3) vs RT (n=3)
# Full result tables, thresholded DEG tables, PCA coordinates and summaries are
# written to figure3/results; panels are rendered by Fig3_plot.py.

suppressPackageStartupMessages({
  library(readxl)
  library(DESeq2)
  library(org.Mm.eg.db)
})

.args <- commandArgs(trailingOnly = FALSE)
.fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
fig_dir <- normalizePath(file.path(dirname(.fa[1]), ".."))
res_dir <- file.path(fig_dir, "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# Contrast 1: in-house mRNA-seq IR vs Sham
raw_ir <- file.path(fig_dir, "raw_data", "Supplementary_Table_S2_differential_expression.xlsx")
df <- as.data.frame(read_excel(raw_ir, sheet = "All_genes_expression", skip = 1,
                               .name_repair = "minimal"))
samples_ir <- c("IR278.Input", "IR279.Input", "IR283.Input",
                "Sham286.Input", "Sham295.Input", "Sham296.Input")
cts <- as.matrix(df[, 2:7])
colnames(cts) <- samples_ir
rownames(cts) <- df$gene_id
storage.mode(cts) <- "integer"
cts <- cts[rowSums(cts) >= 10, ]

coldata_ir <- data.frame(row.names = samples_ir,
                         group = factor(c("IR", "IR", "IR", "Sham", "Sham", "Sham"),
                                        levels = c("Sham", "IR")))
stopifnot(all(rownames(coldata_ir) == colnames(cts)))

dds <- DESeqDataSetFromMatrix(cts, coldata_ir, ~ group)
dds <- DESeq(dds)
res <- as.data.frame(results(dds, contrast = c("group", "IR", "Sham")))
res$Ensembl <- rownames(res)
res$Symbol <- df$gene_name[match(res$Ensembl, df$gene_id)]
res <- res[, c("Ensembl", "Symbol", "baseMean", "log2FoldChange", "lfcSE",
               "stat", "pvalue", "padj")]
res <- res[order(res$pvalue), ]
write.csv(res, file.path(res_dir, "ir_mrna_deseq2_full.csv"), row.names = FALSE)

deg <- res[!is.na(res$padj) & abs(res$log2FoldChange) > 1 & res$padj < 0.05, ]
write.csv(deg, file.path(res_dir, "Fig3B_ir_deg.csv"), row.names = FALSE)
n_up <- sum(deg$log2FoldChange > 1)
n_dn <- sum(deg$log2FoldChange < -1)
sink(file.path(res_dir, "Fig3B_ir_summary.txt"), split = TRUE)
cat("F3AB mRNA-seq DESeq2: IR vs Sham (3v3)\n")
cat("Genes after low-count filter:", nrow(res), "\n")
cat("DEG (|log2FC|>1 & P<0.05):", nrow(deg),
    " Up-regulated:", n_up, " Down-regulated:", n_dn, "\n")
sink()

vsd <- vst(dds, blind = FALSE)
pca <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
pv <- round(100 * attr(pca, "percentVar"))
pca$PC1_percentVar <- pv[1]
pca$PC2_percentVar <- pv[2]
write.csv(pca, file.path(res_dir, "Fig3A_pca.csv"), row.names = FALSE)

# Contrast 2: GSE322673 Cold vs RT
raw_cold <- file.path(fig_dir, "raw_data", "GSE322673", "GSE322673_raw_counts.csv")
cts2 <- read.csv(raw_cold, check.names = FALSE)
samples_cold <- c("Cold_1", "Cold_2", "Cold_3", "RT_1", "RT_2", "RT_3")
cts2$geneID <- sub("\\..*$", "", cts2$geneID)
cts2 <- cts2[!duplicated(cts2$geneID), ]
rownames(cts2) <- cts2$geneID
cts2 <- as.matrix(cts2[, samples_cold])
storage.mode(cts2) <- "integer"
cts2 <- cts2[rowSums(cts2) >= 10, ]

coldata_cold <- data.frame(row.names = samples_cold,
                           condition = factor(c("Cold", "Cold", "Cold", "RT", "RT", "RT"),
                                              levels = c("RT", "Cold")))
stopifnot(all(rownames(coldata_cold) == colnames(cts2)))

dds2 <- DESeqDataSetFromMatrix(cts2, coldata_cold, ~ condition)
dds2 <- DESeq(dds2)
res2 <- as.data.frame(results(dds2, contrast = c("condition", "Cold", "RT")))
res2$Ensembl <- rownames(res2)
sym <- tryCatch({
  suppressMessages(mapIds(org.Mm.eg.db, res2$Ensembl, "SYMBOL", "ENSEMBL", multiVals = "first"))
}, error = function(e) NULL)
res2$Symbol <- if (!is.null(sym)) unname(sym[res2$Ensembl]) else NA
res2 <- res2[, c("Ensembl", "Symbol", "baseMean", "log2FoldChange", "lfcSE",
                 "stat", "pvalue", "padj")]
res2 <- res2[order(res2$pvalue), ]
write.csv(res2, file.path(res_dir, "cold_mrna_deseq2_full.csv"), row.names = FALSE)

deg2 <- res2[!is.na(res2$padj) & abs(res2$log2FoldChange) > 1 & res2$padj < 0.05, ]
write.csv(deg2, file.path(res_dir, "Fig3C_cold_deg.csv"), row.names = FALSE)
n_up2 <- sum(deg2$log2FoldChange > 1)
n_dn2 <- sum(deg2$log2FoldChange < -1)
sink(file.path(res_dir, "Fig3C_cold_summary.txt"), split = TRUE)
cat("GSE322673 DESeq2: Cold vs RT (3v3)\n")
cat("Genes after low-count filter:", nrow(res2), "\n")
cat("DEG (|log2FC|>1 & P<0.05):", nrow(deg2),
    " Up-regulated:", n_up2, " Down-regulated:", n_dn2, "\n")
sink()

cat("Figure 3 DESeq2 data complete\n")
