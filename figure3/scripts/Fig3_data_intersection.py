# Inputs: figure1/results/all_diff_ac4c_sites.csv,
#         figure3/results/ir_mrna_deseq2_full.csv, cold_mrna_deseq2_full.csv
# Outputs: results/Fig3D_triple_genes.csv, Fig3D_core_genes.csv, Fig3D_intersection_summary.txt
from pathlib import Path

import pandas as pd
from scipy.stats import hypergeom

FIG3 = Path(__file__).resolve().parents[1]
PROJ = FIG3.parent
res_dir = FIG3 / "results"

ac4c_df = pd.read_csv(PROJ / "figure1" / "results" / "all_diff_ac4c_sites.csv")
ac4c_genes = set(ac4c_df["GeneName"].dropna().astype(str).str.strip().unique())

ir_df = pd.read_csv(res_dir / "ir_mrna_deseq2_full.csv")
ir_deg = ir_df[(ir_df["padj"].notna()) & (abs(ir_df["log2FoldChange"]) > 1) & (ir_df["padj"] < 0.05)]
ir_genes = set(ir_deg["Symbol"].dropna().astype(str).str.strip().unique())

cold_df = pd.read_csv(res_dir / "cold_mrna_deseq2_full.csv")
cold_deg = cold_df[(cold_df["padj"].notna()) & (abs(cold_df["log2FoldChange"]) > 1) & (cold_df["padj"] < 0.05)]
cold_genes = set(cold_deg["Symbol"].dropna().astype(str).str.strip().unique())

# Universe: genes detected in both transcriptomes
ir_all = set(ir_df["Symbol"].dropna().astype(str).str.strip().unique())
cold_all = set(cold_df["Symbol"].dropna().astype(str).str.strip().unique())
universe_set = ir_all & cold_all
M = len(universe_set)

K_ac4c = len(ac4c_genes & universe_set)
K_ir = len(ir_genes & universe_set)
K_cold = len(cold_genes & universe_set)

ac_ir = len((ac4c_genes & ir_genes) & universe_set)
ac_cold = len((ac4c_genes & cold_genes) & universe_set)
ir_cold = len((ir_genes & cold_genes) & universe_set)
triple = len((ac4c_genes & ir_genes & cold_genes) & universe_set)


def hyper_p(M, K, n, k):
    return hypergeom.sf(k - 1, M, K, n)


p_ac4c_ir = hyper_p(M, K_ac4c, K_ir, ac_ir)
p_ac4c_cold = hyper_p(M, K_ac4c, K_cold, ac_cold)
p_ir_cold = hyper_p(M, K_ir, K_cold, ir_cold)
p_triple = hyper_p(M, ac_ir, K_cold, triple)

triple_genes = sorted((ac4c_genes & ir_genes & cold_genes) & universe_set)
pd.DataFrame({"Symbol": triple_genes}).to_csv(res_dir / "Fig3D_triple_genes.csv", index=False)

# Core genes: concordant up-regulation of ac4C, IR and cold (input to figure4)
max_fold = (ac4c_df.assign(Gene=ac4c_df["GeneName"].dropna().astype(str).str.strip())
            .groupby("Gene")["Foldchange"].max())
ir_up = ir_deg[ir_deg["log2FoldChange"] > 0].copy()
ir_up["Symbol"] = ir_up["Symbol"].astype(str).str.strip()
ir_up = ir_up.drop_duplicates("Symbol").set_index("Symbol")
cold_up = cold_deg[cold_deg["log2FoldChange"] > 0].copy()
cold_up["Symbol"] = cold_up["Symbol"].astype(str).str.strip()
cold_up = cold_up.drop_duplicates("Symbol").set_index("Symbol")
core_rows = []
for g in triple_genes:
    if g in ir_up.index and g in cold_up.index:
        core_rows.append({
            "Gene": g,
            "ac4C_fold": max_fold.get(g),
            "IR_lfc": ir_up.at[g, "log2FoldChange"],
            "IR_padj": ir_up.at[g, "padj"],
            "Cold_lfc": cold_up.at[g, "log2FoldChange"],
            "Cold_padj": cold_up.at[g, "padj"]})
core_out = pd.DataFrame(core_rows).sort_values("ac4C_fold", ascending=False)
core_out.to_csv(res_dir / "Fig3D_core_genes.csv", index=False)

with open(res_dir / "Fig3D_intersection_summary.txt", "w", encoding="utf-8") as f:
    f.write(f"ac4C set size: {len(ac4c_genes)}\n")
    f.write(f"IR set size: {len(ir_genes)}\n")
    f.write(f"Cold set size: {len(cold_genes)}\n")
    f.write(f"Universe (both transcriptomes detected): {M}\n")
    f.write(f"ac4C genes in universe: {K_ac4c}\n")
    f.write(f"IR genes in universe: {K_ir}\n")
    f.write(f"Cold genes in universe: {K_cold}\n")
    f.write(f"ac4C intersect IR (in universe): {ac_ir} (P={p_ac4c_ir:.4g})\n")
    f.write(f"ac4C intersect Cold (in universe): {ac_cold} (P={p_ac4c_cold:.4g})\n")
    f.write(f"IR intersect Cold (in universe): {ir_cold} (P={p_ir_cold:.4g})\n")
    f.write(f"Triple intersection (in universe): {triple} (P={p_triple:.4g})\n")

print(f"Universe: {M}")
print(f"ac4C-IR P: {p_ac4c_ir:.4g}")
print(f"IR-cold P: {p_ir_cold:.4g}")
print(f"ac4C-cold P: {p_ac4c_cold:.4g}")
print(f"Triple P: {p_triple:.4g}")
print("Figure 3 intersection data complete")
