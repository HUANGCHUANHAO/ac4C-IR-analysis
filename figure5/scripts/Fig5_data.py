# Loads the Sham (0h) and IR6h matrices once, runs QC/normalization/embedding/
# leiden clustering and marker-score cell-type annotation, then writes
#   results/GSE247139_processed.h5ad, results/Fig5B_abundance.csv,
#   results/Fig5C_expr_by_celltype.csv and results/run_summary.txt
# Panels A-D are rendered from these outputs by Fig5_plot.py.

import os
import gzip
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

FIG = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RAW_DIR = os.path.join(FIG, "raw_data", "GSE247139")
RES_DIR = os.path.join(FIG, "results")
os.makedirs(RES_DIR, exist_ok=True)

GROUPS = {"0h": "Sham", "6h": "IR6h"}

CORE_GENES = ["S100a9", "Cxcr2", "S100a8", "F10", "Snhg11", "Mmp8", "Arg1"]
INFLAM_GENES = ["Tnf", "Il6", "Ccl2", "Tlr4", "Tgfb1", "Mmp9", "Fn1", "Stat3"]
PLOT_GENES = CORE_GENES + INFLAM_GENES

CELLTYPE_MARKERS = {
    "Cardiomyocyte": ["Tnnt2", "Myh7", "Myh6", "Tnni3"],
    "Fibroblast": ["Col1a1", "Col3a1", "Pdgfra", "Dcn"],
    "Endothelial": ["Pecam1", "Cdh5", "Kdr"],
    "SMC_Pericyte": ["Acta2", "Myh11", "Rgs5", "Pdgfrb"],
    "MonoMacro": ["Adgre1", "Cd68", "Lgals3", "Mertk", "Lyz2", "C1qa"],
    "Neutrophil": ["S100a8", "S100a9", "Ly6g", "Csf3r"],
    "B_cell": ["Cd79a", "Ms4a1", "Cd19"],
    "T_cell": ["Cd3d", "Cd3e", "Cd8a", "Cd4"],
    "NK_cell": ["Nkg7", "Klrb1c", "Gzma"],
    "DC": ["Itgax", "Flt3", "Clec9a"],
}

CT_DISP = {
    "Cardiomyocyte": "Cardiomyocyte",
    "Fibroblast": "Fibroblast",
    "Endothelial": "Endothelial",
    "SMC_Pericyte": "SMC/Pericyte",
    "MonoMacro": "Mono/Macrophage",
    "Neutrophil": "Neutrophil",
    "B_cell": "B cell",
    "T_cell": "T cell",
    "NK_cell": "NK cell",
    "DC": "Dendritic cell",
    "Other": "Other",
}

CT_ORDER = [
    "Neutrophil", "Mono/Macrophage", "B cell", "T cell", "NK cell", "Dendritic cell",
    "Cardiomyocyte", "Fibroblast", "Endothelial", "SMC/Pericyte", "Other"
]

print("Loading GSE247139 matrices...")
adatas = []
for h, cond in sorted(GROUPS.items()):
    mtx_path = os.path.join(RAW_DIR, f"GSE247139_{h}_matrix.mtx.gz")
    bc_path = os.path.join(RAW_DIR, f"GSE247139_{h}_barcodes.tsv.gz")
    feat_path = os.path.join(RAW_DIR, f"GSE247139_{h}_features.tsv.gz")
    with gzip.open(mtx_path, "rb") as fh:
        mat = mmread(fh).tocsr().T
    with gzip.open(bc_path, "rt") as fh:
        bc = [l.strip().split("\t")[0] for l in fh]
    with gzip.open(feat_path, "rt") as fh:
        feat = [l.strip().split("\t") for l in fh]
    genes = [c[1] if len(c) > 1 else c[0] for c in feat]
    a = sc.AnnData(X=mat, obs=pd.DataFrame(index=[f"{h}_{b}" for b in bc]),
                   var=pd.DataFrame(index=genes))
    a.obs["condition"] = cond
    a.obs["bc_suffix"] = [b.rsplit("-", 1)[-1] for b in bc]
    a.var_names_make_unique()
    adatas.append(a)
    print(f"  {cond}: {a.n_obs} cells x {a.n_vars} genes")

adata = sc.concat(adatas, join="outer", index_unique=None)
adata.var_names_make_unique()
print(f"Combined: {adata.n_obs} cells x {adata.n_vars} genes")

adata.var["mt"] = adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
n_cells_before = adata.n_obs
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[(adata.obs.n_genes_by_counts < 5000) & (adata.obs.pct_counts_mt < 20)].copy()
print(f"QC: {n_cells_before} -> {adata.n_obs} cells")

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata.copy()

marker_genes = set().union(*CELLTYPE_MARKERS.values()) | set(PLOT_GENES)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="condition")
keep = adata.var.highly_variable | adata.var_names.isin(marker_genes)
adata = adata[:, keep].copy()
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_neighbors=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.6)

ad_score = adata.raw[:, adata.var_names].to_adata()
for ct, gs in CELLTYPE_MARKERS.items():
    use = [g for g in gs if g in ad_score.var_names]
    if use:
        sc.tl.score_genes(ad_score, use, score_name=f"score_{ct}")
score_cols = [c for c in ad_score.obs.columns if c.startswith("score_")]
adata.obs[score_cols] = ad_score.obs[score_cols].values

cl = adata.obs.groupby("leiden")[score_cols].mean()
SCORE_MIN = 0.05
assign = {}
for clu, row in cl.iterrows():
    best = row.idxmax()
    assign[clu] = best.replace("score_", "") if row[best] >= SCORE_MIN else "Other"
adata.obs["celltype"] = adata.obs["leiden"].map(assign)
adata.obs["celltype_disp"] = adata.obs["celltype"].map(CT_DISP)

cell_order = adata.obs["celltype_disp"].value_counts().index.tolist()
final_order = [c for c in CT_ORDER if c in cell_order] + \
              [c for c in cell_order if c not in CT_ORDER]
adata.obs["celltype_disp"] = pd.Categorical(
    adata.obs["celltype_disp"], categories=final_order, ordered=True)

adata.write(os.path.join(RES_DIR, "GSE247139_processed.h5ad"))

abund = pd.crosstab(adata.obs["celltype_disp"], adata.obs["condition"])
abund_pct = pd.crosstab(adata.obs["celltype_disp"], adata.obs["condition"],
                        normalize="columns") * 100
abund_out = pd.DataFrame({
    "celltype": abund.index,
    "Sham_n": abund["Sham"].values,
    "IR6h_n": abund["IR6h"].values,
    "Sham_pct": abund_pct["Sham"].values.round(2),
    "IR6h_pct": abund_pct["IR6h"].values.round(2),
})
abund_out.to_csv(os.path.join(RES_DIR, "Fig5B_abundance.csv"), index=False)
print("Saved: Fig5B_abundance.csv")

raw_ad = adata.raw.to_adata()
expr_rows = []
for ct in cell_order:
    sub = raw_ad[raw_ad.obs["celltype_disp"] == ct]
    for g in CORE_GENES:
        if g not in sub.var_names:
            continue
        vec = sub[:, g].X
        if not isinstance(vec, np.ndarray):
            vec = vec.toarray()
        vals = np.asarray(vec).ravel()
        expr_rows.append({
            "celltype": ct,
            "gene": g,
            "mean_expr": float(vals.mean()),
            "pct_expr": (vals > 0).mean() * 100,
            "n_cells": sub.n_obs})
pd.DataFrame(expr_rows).to_csv(
    os.path.join(RES_DIR, "Fig5C_expr_by_celltype.csv"), index=False)
print("Saved: Fig5C_expr_by_celltype.csv")

core_present = [g for g in CORE_GENES if g in adata.var_names]
core_missing = [g for g in CORE_GENES if g not in adata.var_names]
with open(os.path.join(RES_DIR, "run_summary.txt"), "w", encoding="utf-8") as f:
    f.write("Fig5 run summary\n\n")
    f.write(f"Total cells after QC: {adata.n_obs}\n")
    f.write(f"Genes after QC/HVG: {adata.n_vars}\n")
    f.write(f"Conditions: {adata.obs['condition'].value_counts().to_dict()}\n")
    f.write("Cell types:\n")
    f.write(abund_out.to_string(index=False) + "\n\n")
    f.write(f"Core genes present: {core_present}\n")
    f.write(f"Core genes missing: {core_missing}\n")
print("Saved: run_summary.txt")
print("Figure 5 data complete")
