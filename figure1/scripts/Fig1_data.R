# Figure 1 data preparation: differential ac4C site tables.
# Exports results/all_diff_ac4c_sites.csv and results/Fig1C_chrom_distribution.csv.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
})

.args <- commandArgs(trailingOnly = FALSE)
.fa <- sub("^--file=", "", .args[grep("^--file=", .args)])
fig_dir <- normalizePath(file.path(dirname(.fa[1]), ".."))

raw_file <- file.path(fig_dir, "raw_data", "Supplementary_Table_S1_differential_ac4C_sites.xlsx")
res_dir <- file.path(fig_dir, "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

up_sites <- read_excel(raw_file, sheet = "up.IR_vs_Sham", skip = 28)
down_sites <- read_excel(raw_file, sheet = "down.IR_vs_Sham", skip = 28)
up_sites$direction <- "up"
down_sites$direction <- "down"

diff_sites <- bind_rows(up_sites, down_sites) %>%
  filter(!is.na(GeneName), GeneName != "")

write.csv(diff_sites, file.path(res_dir, "all_diff_ac4c_sites.csv"), row.names = FALSE)

chrom_counts <- diff_sites %>%
  count(chrom, direction) %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0)

write.csv(chrom_counts, file.path(res_dir, "Fig1C_chrom_distribution.csv"), row.names = FALSE)

cat("Figure 1 data done:", nrow(diff_sites), "differential sites\n")
