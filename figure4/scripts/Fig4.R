# Figure 4 human evaluation of the six direction-concordant core genes.
#   A: core genes ranked by ac4C fold change (reads figure3/results/Fig3D_core_genes.csv)
#   B: single-gene ROC in GSE66360 (discovery)
#   C: single-gene ROC in GSE48060 (external validation)
#   D/E: S100A9/S100A8/CXCR2 boxplots in GSE66360 and GSE48060

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(GEOquery)
  library(pROC)
  library(hgu133plus2.db)
  library(reshape2)
})

.args <- commandArgs(trailingOnly = FALSE)
.fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
fig_dir <- normalizePath(file.path(dirname(.fa[1]), ".."))
proj_dir <- normalizePath(file.path(fig_dir, ".."))
raw_dir <- file.path(fig_dir, "raw_data")
res_dir <- file.path(fig_dir, "results")
fig_dir_out <- file.path(fig_dir, "report")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir_out, showWarnings = FALSE, recursive = TRUE)

genes_h <- c("S100A9", "CXCR2", "S100A8", "F10", "SNHG11", "ARG1")
box_genes <- c("S100A9", "S100A8", "CXCR2")

# Panel A: direction-concordant core genes ranked by ac4C fold change
core <- read.csv(file.path(proj_dir, "figure3", "results", "Fig3D_core_genes.csv"),
                 stringsAsFactors = FALSE)
core$Gene <- factor(core$Gene, levels = rev(core$Gene))

p_a <- ggplot(core, aes(x = ac4C_fold, y = Gene)) +
  geom_bar(stat = "identity", width = 0.65, fill = "#3A7CA6") +
  geom_text(aes(label = sprintf("%.1f", ac4C_fold)), hjust = -0.15, size = 3.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "ac4C Fold change (IR vs Sham)",
    y = NULL,
    title = "Direction-consistent core genes\n(ac4C up + mRNA up + cold up, adjP<0.05, n=6)"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 13),
    plot.title = element_text(size = 15, hjust = 0.5, face = "bold")
  )
ggsave(file.path(fig_dir_out, "Fig4A_core_genes.png"), p_a, width = 6, height = 4, dpi = 300)
write.csv(core, file.path(res_dir, "Fig4A_core_genes.csv"), row.names = FALSE)

# Shared GEO series-matrix loader and single-gene ROC (panels B and C)
load_mat <- function(gse) {
  g <- getGEO(filename = file.path(raw_dir, paste0(gse, "_series_matrix.txt")),
              GSEMatrix = TRUE, getGPL = FALSE)
  e <- exprs(g)
  if (max(e, na.rm = TRUE) > 100) e <- log2(e)
  syms <- unname(mapIds(hgu133plus2.db, rownames(e), "SYMBOL", "PROBEID", multiVals = "first"))
  ok <- !is.na(syms)
  e <- e[ok, ]
  syms <- syms[ok]
  o <- order(rowMeans(e), decreasing = TRUE)
  e <- e[o, ]
  syms <- syms[o]
  dup <- duplicated(syms)
  e <- e[!dup, ]
  rownames(e) <- syms[!dup]
  list(e = e, pd = pData(g))
}

run_roc <- function(gse, label, csv_name, png_name) {
  dat <- load_mat(gse)
  e <- dat$e
  pd <- dat$pd
  ch <- apply(pd[, grep("disease", colnames(pd), ignore.case = TRUE), drop = FALSE], 1,
              function(x) paste(x, collapse = "|"))
  grp <- ifelse(grepl("control", ch, ignore.case = TRUE), "Control",
                ifelse(grepl("infarction|patient", ch, ignore.case = TRUE), "AMI", NA))
  keep <- !is.na(grp)
  e <- e[, keep]
  grp <- factor(grp[keep], levels = c("Control", "AMI"))
  cat(gse, "group:", paste(names(table(grp)), table(grp), collapse = " vs "), "\n")

  res <- data.frame(Gene = genes_h, n_probe = NA, AUC = NA, CI_low = NA, CI_high = NA,
                    P_wilcox = NA, log2FC_AMI_vs_Ctrl = NA)
  roc_list <- list()
  for (i in seq_along(genes_h)) {
    gh <- genes_h[i]
    if (!(gh %in% rownames(e))) next
    x <- e[gh, ]
    r <- roc(grp, x, quiet = TRUE, levels = c("Control", "AMI"), direction = "<")
    ci <- ci.auc(r)
    res$n_probe[i] <- 1
    res$AUC[i] <- as.numeric(auc(r))
    res$CI_low[i] <- ci[1]
    res$CI_high[i] <- ci[3]
    res$P_wilcox[i] <- wilcox.test(x ~ grp)$p.value
    res$log2FC_AMI_vs_Ctrl[i] <- mean(x[grp == "AMI"]) - mean(x[grp == "Control"])
    roc_list[[gh]] <- r
  }
  res <- res[!is.na(res$AUC), ]
  res$Padj <- p.adjust(res$P_wilcox, "BH")
  write.csv(res, file.path(res_dir, csv_name), row.names = FALSE)

  png(file.path(fig_dir_out, png_name), width = 1600, height = 1400, res = 150, pointsize = 23)
  par(cex.main = 1.1)
  plot(roc_list[[1]], col = 1, main = paste("ROC -", gse, label), legacy.axes = TRUE)
  cols <- rainbow(length(roc_list))
  for (i in seq_along(roc_list)) plot(roc_list[[i]], col = cols[i], add = TRUE)
  legend("bottomright",
         legend = paste0(names(roc_list), " AUC=",
                         sapply(roc_list, function(r) round(as.numeric(auc(r)), 3))),
         col = cols, lty = 1, cex = 0.87, bty = "n")
  dev.off()
}

run_roc("GSE66360", "(CEC, discovery)", "Fig4B_roc_gse66360.csv", "Fig4B_ROC_curve_GSE66360.png")
run_roc("GSE48060", "(blood, validation)", "Fig4C_roc_gse48060.csv", "Fig4C_ROC_curve_GSE48060.png")

# Panels D/E: expression boxplots of S100A9/S100A8/CXCR2
COL_CTRL <- "#4575B4"
COL_AMI <- "#D73027"

plot_box <- function(expr_df, dataset_label, group_label, out_file) {
  expr_df$group <- factor(expr_df$group, levels = c("Control", "AMI"))
  expr_df$gene <- factor(expr_df$gene, levels = box_genes)
  ann <- do.call(rbind, lapply(levels(expr_df$gene), function(g) {
    sub <- expr_df[expr_df$gene == g, ]
    x_ctrl <- sub$value[sub$group == "Control"]
    x_ami <- sub$value[sub$group == "AMI"]
    p <- wilcox.test(x_ami, x_ctrl)$p.value
    stars <- ifelse(p < 0.001, "***",
                    ifelse(p < 0.01, "**",
                           ifelse(p < 0.05, "*", "ns")))
    y_max <- max(sub$value, na.rm = TRUE)
    data.frame(gene = g, p = p, stars = stars, y_pos = y_max * 1.05,
               stringsAsFactors = FALSE)
  }))
  ann$gene <- factor(ann$gene, levels = levels(expr_df$gene))

  p <- ggplot(expr_df, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 1.5, outlier.color = "black",
                 outlier.fill = "white", color = "black", width = 0.6) +
    scale_fill_manual(values = c("Control" = COL_CTRL, "AMI" = COL_AMI)) +
    facet_wrap(~ gene, nrow = 1, scales = "free_y") +
    geom_text(data = ann, aes(x = 1.5, y = y_pos, label = stars),
              inherit.aes = FALSE, size = 7, fontface = "bold", vjust = 0) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(title = paste0(dataset_label, " | ", group_label), x = NULL, y = "log2 expression") +
    theme_classic(base_size = 17) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 19),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold.italic", size = 17),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "none",
      panel.spacing = unit(1.2, "lines")
    )
  ggsave(out_file, p, width = 8, height = 4, dpi = 300)
  cat("Saved:", out_file, "\n")
  ann
}

load_and_extract <- function(gse) {
  g <- getGEO(filename = file.path(raw_dir, paste0(gse, "_series_matrix.txt")),
              GSEMatrix = TRUE, getGPL = FALSE)
  e <- exprs(g)
  if (max(e, na.rm = TRUE) > 100) e <- log2(e)
  syms <- unname(mapIds(hgu133plus2.db, rownames(e), "SYMBOL", "PROBEID", multiVals = "first"))
  ok <- !is.na(syms)
  e <- e[ok, ]
  syms <- syms[ok]
  o <- order(rowMeans(e), decreasing = TRUE)
  e <- e[o, ]
  syms <- syms[o]
  dup <- duplicated(syms)
  e <- e[!dup, ]
  rownames(e) <- syms[!dup]
  pd <- pData(g)
  ch <- apply(pd[, grep("disease", colnames(pd), ignore.case = TRUE), drop = FALSE], 1,
              function(x) paste(x, collapse = "|"))
  grp <- ifelse(grepl("control", ch, ignore.case = TRUE), "Control",
                ifelse(grepl("infarction|patient", ch, ignore.case = TRUE), "AMI", NA))
  keep <- !is.na(grp)
  e <- e[, keep, drop = FALSE]
  grp <- grp[keep]
  list(e = e, grp = grp, gse = gse)
}

d66 <- load_and_extract("GSE66360")
ev <- lapply(box_genes, function(g) {
  if (!(g %in% rownames(d66$e))) return(NULL)
  data.frame(sample = colnames(d66$e), group = d66$grp, gene = g,
             value = d66$e[g, ], stringsAsFactors = FALSE)
})
ev <- do.call(rbind, ev[!sapply(ev, is.null)])
ann66 <- plot_box(ev, "GSE66360 (CEC, n=99)", "AMI vs Control",
                  file.path(fig_dir_out, "Fig4D_boxplot_GSE66360.png"))
write.csv(ann66, file.path(res_dir, "Fig4D_boxstats_gse66360.csv"), row.names = FALSE)

d48 <- load_and_extract("GSE48060")
ev48 <- lapply(box_genes, function(g) {
  if (!(g %in% rownames(d48$e))) return(NULL)
  data.frame(sample = colnames(d48$e), group = d48$grp, gene = g,
             value = d48$e[g, ], stringsAsFactors = FALSE)
})
ev48 <- do.call(rbind, ev48[!sapply(ev48, is.null)])
ann48 <- plot_box(ev48, "GSE48060 (whole blood, n=52)", "post-first-MI vs Control",
                  file.path(fig_dir_out, "Fig4E_boxplot_GSE48060.png"))
write.csv(ann48, file.path(res_dir, "Fig4E_boxstats_gse48060.csv"), row.names = FALSE)

cat("Figure 4 complete\n")
