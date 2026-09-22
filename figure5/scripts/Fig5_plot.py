# Inputs : figure5/results/GSE247139_processed.h5ad (built by Fig5_data.py),
#          figure5/results/Fig5B_abundance.csv
# Outputs: report/Fig5A_UMAP.png, Fig5B_celltype_abundance.png,
#          Fig5C_core_genes_dotplot.png, Fig5D_myeloid_violin.png
import os
import sys
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "common")))
from plot_style import W_THIRD

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

FIG = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RES = os.path.join(FIG, "results")
OUT = os.path.join(FIG, "report")
os.makedirs(OUT, exist_ok=True)
H5AD = os.path.join(RES, "GSE247139_processed.h5ad")

CORE6 = ["S100a9", "Cxcr2", "S100a8", "F10", "Snhg11", "Arg1"]
MYE4 = ["S100a9", "S100a8", "Cxcr2", "Arg1"]
SHAM_C, IR_C = "#1B9E77", "#D73027"


def _font(t):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "figure.dpi": 300, "savefig.dpi": 300,
        "axes.titleweight": "bold",
        "axes.titlesize": t["title"], "axes.labelsize": t["label"],
        "xtick.labelsize": t["tick"], "ytick.labelsize": t["tick"],
        "legend.fontsize": t["legend"],
    })


def panel_umap():
    _font(dict(title=19, label=16, tick=15, legend=14))
    adata = sc.read_h5ad(H5AD)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sc.pl.umap(adata, color="condition", ax=axes[0], show=False,
               title="Condition", frameon=False,
               palette={"Sham": SHAM_C, "IR6h": IR_C})
    sc.pl.umap(adata, color="celltype_disp", ax=axes[1], show=False,
               title="Cell type", frameon=False, legend_loc="right margin")
    fig.suptitle("GSE247139: single-cell landscape of heart IR",
                 fontsize=20, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(OUT, "Fig5A_UMAP.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved: Fig5A_UMAP.png")


def panel_abundance():
    _font(dict(title=18, label=16, tick=16, legend=16))
    plot_abund = pd.read_csv(os.path.join(RES, "Fig5B_abundance.csv"))
    plot_abund = plot_abund[plot_abund["celltype"] != "Other"]
    plot_abund = plot_abund.sort_values("IR6h_pct", ascending=True)
    y = np.arange(len(plot_abund))
    width = 0.35
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(y - width/2, plot_abund["Sham_pct"], width, label="Sham", color=SHAM_C)
    ax.barh(y + width/2, plot_abund["IR6h_pct"], width, label="IR6h", color=IR_C)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_abund["celltype"])
    ax.set_xlabel("Percentage of cells (%)")
    ax.set_title("Cell-type composition in Sham vs IR6h")
    ax.legend(frameon=False)
    sns.despine(ax=ax)
    plt.tight_layout()
    fig.savefig(os.path.join(OUT, "Fig5B_celltype_abundance.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved: Fig5B_celltype_abundance.png")


def panel_dotplot():
    _font(dict(title=17, label=16, tick=13, legend=13))
    adata = sc.read_h5ad(H5AD)
    gcol = "celltype_disp" if "celltype_disp" in adata.obs.columns else "celltype"
    ad = adata.raw.to_adata()
    genes_c = [g for g in CORE6 if g in ad.var_names]
    sc.pl.dotplot(ad, genes_c, groupby=gcol, use_raw=False, color_map="Reds",
                  title="Core gene expression across cell types",
                  save=False, show=False)
    fig = plt.gcf()
    fig.set_size_inches(W_THIRD, W_THIRD / 1.42)
    fig.savefig(os.path.join(OUT, "Fig5C_core_genes_dotplot.png"),
                dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("saved Fig5C")


def panel_violin():
    _font(dict(title=17, label=17, tick=14, legend=14))
    adata = sc.read_h5ad(H5AD)
    mye = adata[adata.obs["celltype"].isin(["Neutrophil", "MonoMacro"])].copy()
    genes_d = [g for g in MYE4 if g in mye.raw.var_names]
    fig, axes = plt.subplots(2, 2, figsize=(W_THIRD, W_THIRD / 1.32))
    axes = axes.flatten()
    for i, g in enumerate(genes_d):
        ax = axes[i]
        sc.pl.violin(mye, g, groupby="condition", ax=ax, show=False,
                     stripplot=False, use_raw=True,
                     palette={"Sham": SHAM_C, "IR6h": IR_C})
        ax.set_title(g, fontsize=17, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("")
        vir = mye[mye.obs["condition"] == "IR6h"].raw[:, g].X.toarray().flatten()
        vsh = mye[mye.obs["condition"] == "Sham"].raw[:, g].X.toarray().flatten()
        ymax = max(vir.max(), vsh.max())
        ax.set_ylim(0, ymax * 1.2)
        ax.tick_params(labelsize=14)
    for i in range(len(genes_d), len(axes)):
        axes[i].axis("off")
    fig.suptitle("Core genes in neutrophils, monocytes/macrophages",
                 fontsize=17, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(OUT, "Fig5D_myeloid_violin.png"),
                dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("saved Fig5D descriptive violin (no statistical test)")


if __name__ == "__main__":
    panel_umap()
    panel_abundance()
    panel_dotplot()
    panel_violin()
