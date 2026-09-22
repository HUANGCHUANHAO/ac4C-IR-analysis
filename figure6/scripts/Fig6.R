# Figure 6 S100A9 + CXCR2 combined diagnostic model for AMI.
# Training: GSE66360 (49 AMI / 50 control); external validation: GSE48060 (31 AMI / 21 control).
# Single-gene ROC references: S100A8, S100A9, CXCR2.
# Panels A-D and all result tables are written to figure6/report and figure6/results;
# the one-/two-/three-gene sensitivity comparison is exported as Table1_model_variants.csv.

suppressPackageStartupMessages({
  library(GEOquery)
  library(pROC)
  library(hgu133plus2.db)
  library(rms)
})

set.seed(20260905)

.args <- commandArgs(trailingOnly = FALSE)
.fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
proj_dir <- normalizePath(file.path(dirname(.fa[1]), "..", ".."))
raw_dir <- file.path(proj_dir, "figure4", "raw_data")
out_dir <- file.path(proj_dir, "figure6", "results")
fig_dir <- file.path(proj_dir, "figure6", "report")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

marker_genes <- c("S100A8", "S100A9", "CXCR2")

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
  pd <- pData(g)
  ch <- apply(pd[, grep("disease", colnames(pd), ignore.case = TRUE), drop = FALSE],
              1, function(x) paste(x, collapse = "|"))
  grp <- ifelse(grepl("control", ch, ignore.case = TRUE), 0,
                ifelse(grepl("infarction|patient", ch, ignore.case = TRUE), 1, NA))
  keep <- !is.na(grp)
  d <- as.data.frame(t(e[marker_genes, keep, drop = FALSE]))
  d$y <- grp[keep]
  rownames(d) <- colnames(e)[keep]
  d
}

tr <- load_mat("GSE66360")
va <- load_mat("GSE48060")
cat("train", sum(tr$y == 1), "/", sum(tr$y == 0),
    "; valid", sum(va$y == 1), "/", sum(va$y == 0), "\n")

# Combined two-gene logistic model and external logistic recalibration
fit <- glm(y ~ S100A9 + CXCR2, family = binomial(), data = tr)
tr$p <- predict(fit, type = "response")
va$p_raw <- predict(fit, newdata = va, type = "response")

recal <- glm(va$y ~ I(qlogis(va$p_raw)), family = binomial())
va$p <- predict(recal, type = "response")
cat("recalibration intercept, slope:", coef(recal), "\n")

coef_tab <- as.data.frame(coef(summary(fit)))
coef_tab$Term <- rownames(coef_tab)
coef_tab$OR <- exp(coef_tab[, 1])
write.csv(coef_tab, file.path(out_dir, "model_coefficients.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
print(coef_tab)

mk <- function(d, v) roc(d$y, d[[v]], quiet = TRUE, levels = c(0, 1), direction = "<")
Rtr <- list(S100A8 = mk(tr, "S100A8"), S100A9 = mk(tr, "S100A9"),
            CXCR2 = mk(tr, "CXCR2"), Model = mk(tr, "p"))
Rva <- list(S100A8 = mk(va, "S100A8"), S100A9 = mk(va, "S100A9"),
            CXCR2 = mk(va, "CXCR2"), Model = mk(va, "p_raw"))

auc_tab <- do.call(rbind, lapply(names(Rtr), function(nm) {
  ct <- ci.auc(Rtr[[nm]])
  cv <- ci.auc(Rva[[nm]])
  data.frame(Feature = nm,
             Train_AUC = as.numeric(auc(Rtr[[nm]])), Train_CIlo = ct[1], Train_CIhi = ct[3],
             Valid_AUC = as.numeric(auc(Rva[[nm]])), Valid_CIlo = cv[1], Valid_CIhi = cv[3])
}))
write.csv(auc_tab, file.path(out_dir, "Fig6A_auc.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
print(auc_tab, digits = 3)
cat("DeLong train Model vs S100A9:", roc.test(Rtr$Model, Rtr$S100A9)$p.value, "\n")
cat("DeLong valid Model vs S100A9:", roc.test(Rva$Model, Rva$S100A9)$p.value, "\n")

coords_tr <- coords(Rtr$Model, "best", best.method = "closest.topleft",
                    ret = c("threshold", "sensitivity", "specificity"))
write.csv(coords_tr, file.path(out_dir, "Fig6A_cutpoint.csv"), fileEncoding = "UTF-8")

save_png <- function(file, w, h, fn, ps = 17.4) {
  png(file.path(fig_dir, file), width = w * 300, height = h * 300, res = 300, pointsize = ps)
  fn()
  dev.off()
}
plot_roc <- function(R, title) {
  cols <- c("grey60", "grey30", "steelblue3", "red2")
  ltys <- c(2, 2, 2, 1)
  plot(R[[1]], legacy.axes = TRUE, col = cols[1], lty = ltys[1], lwd = 2, main = title)
  for (i in 2:4) plot(R[[i]], col = cols[i], lty = ltys[i], lwd = 2, add = TRUE)
  abline(0, 1, col = "grey85", lty = 3)
  legend("bottomright", bty = "n", lwd = 2, col = cols, lty = ltys, cex = 0.95,
         legend = paste0(names(R), " AUC=",
                         sapply(R, function(x) sprintf("%.3f", as.numeric(auc(x))))))
}
save_png("Fig6A_ROC_train.png", 6.2, 6, function() plot_roc(Rtr, "Training cohort GSE66360 (n=99)"))
save_png("Fig6A_ROC_valid.png", 6.2, 6, function() plot_roc(Rva, "External validation GSE48060 (n=52)"))

# Panel B nomogram
dd <- datadist(tr)
options(datadist = "dd")
lrm_fit <- lrm(y ~ S100A9 + CXCR2, data = tr, x = TRUE, y = TRUE)
sink(file.path(out_dir, "lrm_model.txt"))
print(lrm_fit)
sink()
nom <- nomogram(lrm_fit, fun = plogis, funlabel = "Probability of AMI", lp = FALSE,
                fun.at = c(0.1, 0.3, 0.5, 0.7, 0.9))
save_png("Fig6B_nomogram.png", 8.5, 4.6,
         function() plot(nom, xfrac = 0.32, cex.axis = 0.95, cex.var = 1.0), ps = 13)

# Panel C calibration
cal <- calibrate(lrm_fit, method = "boot", B = 1000, data = tr)
vp_raw <- val.prob(va$p_raw, va$y, pl = FALSE)
vp_rc <- val.prob(va$p, va$y, pl = FALSE)
sink(file.path(out_dir, "Fig6C_calibration.txt"))
cat("== external raw frozen ==\n")
print(vp_raw)
cat("\n== external logistic-recalibrated ==\n")
print(vp_rc)
cat("\n== train bootstrap calibrate ==\n")
print(cal)
sink()
save_png("Fig6C_calibration.png", 10.5, 5.25, function() {
  par(mfrow = c(1, 2), mar = c(5, 4.8, 3.2, 1.2))
  plot(cal, legend = FALSE, subtitles = FALSE,
       main = "Train GSE66360 (bootstrap B=1000)",
       xlab = "Predicted probability", ylab = "Observed probability", cex.main = 1.0)
  leg_lab <- c("Apparent", "Bias-corrected", "Ideal", "C.L.")
  leg_cex <- 0.85
  leg_x <- 0.97 - max(strwidth(leg_lab, cex = leg_cex, units = "user")) -
    strwidth("MMM", cex = leg_cex, units = "user")
  legend(leg_x, 0.33, leg_lab, lty = c(3, 1, 2, 1),
         col = c(rep("black", 3), "gray80"), bty = "n", cex = leg_cex)
  val.prob(va$p, va$y, statloc = list(x = 0.0, y = 0.98),
           legendloc = list(x = 0.50, y = 0.27),
           cex = 0.58, xlab = "Predicted Probability", ylab = "Actual Probability")
  title(main = "External GSE48060 (recalibrated)", cex.main = 1.0)
}, ps = 14)

# Panel D decision curve analysis
dca <- function(p, y, grid = seq(0.01, 0.99, 0.01)) {
  n <- length(y)
  prev <- mean(y)
  out <- t(sapply(grid, function(pt) {
    tp <- sum(p >= pt & y == 1)
    fp <- sum(p >= pt & y == 0)
    c(pt, tp / n - fp / n * (pt / (1 - pt)), prev - (1 - prev) * (pt / (1 - pt)))
  }))
  data.frame(threshold = out[, 1], model = out[, 2],
             treat_all = out[, 3], treat_none = 0)
}
dc_tr <- dca(tr$p, tr$y)
dc_va <- dca(va$p, va$y)
write.csv(dc_tr, file.path(out_dir, "Fig6D_dca_train.csv"), row.names = FALSE)
write.csv(dc_va, file.path(out_dir, "Fig6D_dca_valid.csv"), row.names = FALSE)
saveRDS(list(tr = tr, va = va, cal = cal, vp_rc = vp_rc, dc_tr = dc_tr,
             dc_va = dc_va, nom = nom, Rtr = Rtr, Rva = Rva),
        file.path(out_dir, "fig6_plot_objects.rds"))
plot_dca <- function(dc, title) {
  plot(dc$threshold, dc$model, type = "l", lwd = 2.2, col = "red2", xlim = c(0, 1),
       ylim = c(-0.05, max(c(dc$model, dc$treat_all), na.rm = TRUE) * 1.05),
       xlab = "Threshold probability", ylab = "Net benefit", main = title)
  lines(dc$threshold, dc$treat_all, lwd = 1.8, lty = 2, col = "steelblue3")
  abline(h = 0, lwd = 1.6, lty = 3, col = "grey40")
  legend("topright", bty = "n", cex = 0.9, lwd = c(2.2, 1.8, 1.6), lty = c(1, 2, 3),
         col = c("red2", "steelblue3", "grey40"),
         legend = c("2-gene model", "Treat all", "Treat none"))
}
save_png("Fig6D_DCA_train.png", 5.4, 5, function() plot_dca(dc_tr, "DCA - train GSE66360"))
save_png("Fig6D_DCA_valid.png", 5.4, 5, function() plot_dca(dc_va, "DCA - external GSE48060"))
save_png("Fig6D_DCA_combined.png", 11, 5.2, function() {
  par(mfrow = c(1, 2), mar = c(5, 5, 3.4, 1.2))
  plot_dca(dc_tr, "DCA - train GSE66360")
  plot_dca(dc_va, "DCA - external GSE48060")
})

write.csv(data.frame(sample = rownames(va), y = va$y,
                     p_raw_frozen = va$p_raw, p_recalibrated = va$p),
          file.path(out_dir, "external_probabilities.csv"), row.names = FALSE)

# Sensitivity: one-/two-/three-gene model variants (Table 1)
cat("cor S100A8/A9 train:", cor(tr$S100A8, tr$S100A9),
    " valid:", cor(va$S100A8, va$S100A9), "\n")
zmean <- sapply(marker_genes, function(g) mean(tr[[g]]))
zsd <- sapply(marker_genes, function(g) sd(tr[[g]]))
add_composite <- function(d) {
  d$SA <- ((d$S100A8 - zmean[1]) / zsd[1] + (d$S100A9 - zmean[2]) / zsd[2]) / 2
  d$CXz <- (d$CXCR2 - zmean[3]) / zsd[3]
  d
}
tr <- add_composite(tr)
va <- add_composite(va)

variants <- list(
  "S100A9 only"        = y ~ S100A9,
  "A8+A9"              = y ~ S100A8 + S100A9,
  "A9+CXCR2 (final)"   = y ~ S100A9 + CXCR2,
  "A8+A9+CXCR2"        = y ~ S100A8 + S100A9 + CXCR2,
  "SA-composite+CXCR2" = y ~ SA + CXz
)

res <- data.frame()
for (nm in names(variants)) {
  f <- glm(variants[[nm]], family = binomial(), data = tr)
  pt <- predict(f, type = "response")
  pv <- predict(f, newdata = va, type = "response")
  at <- as.numeric(auc(roc(tr$y, pt, quiet = TRUE, direction = "<")))
  av <- as.numeric(auc(roc(va$y, pv, quiet = TRUE, direction = "<")))
  rc <- glm(va$y ~ I(qlogis(pv)), family = binomial())
  pr <- predict(rc, type = "response")
  vr <- val.prob(pr, va$y, pl = FALSE)
  ct <- coef(summary(f))
  res <- rbind(res, data.frame(
    Variant = nm, Train_AUC = round(at, 3), Valid_AUC = round(av, 3),
    ValidRecal_Brier = round(vr[["Brier"]], 3),
    ValidRecal_Intercept = round(vr[["Intercept"]], 3),
    ValidRecal_Slope = round(vr[["Slope"]], 3),
    Coefs = paste(rownames(ct), round(ct[, 1], 3), collapse = "; ")))
}
print(res, right = FALSE)
write.csv(res, file.path(out_dir, "Table1_model_variants.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("Figure 6 complete\n")
