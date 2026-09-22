# Inputs : figure3/results tables (Fig3_data_deseq2.R, Fig3_data_intersection.py)
# Outputs: report/F3A_PCA.png, F3B_volcano.png, GSE322673_volcano.png, F3D_triple_venn.png
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from scipy.stats import hypergeom

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "common")))
from plot_style import UP, DOWN, GREY, W_HALF, save, fmt_pow

FIG = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PROJ = os.path.abspath(os.path.join(FIG, ".."))
RES = os.path.join(FIG, "results")
R1 = os.path.join(PROJ, "figure1", "results")
OUT = os.path.join(FIG, "report")
os.makedirs(OUT, exist_ok=True)
LFC, Q, FLOOR = 1.0, 0.05, 1e-300


def panel_pca():
    # DESeq2-VST PCA coordinates are exported to Fig3A_pca.csv; the approved
    # panel is redrawn here from fixed coordinates so sample labels are not clipped.
    pts = [
        ("IR283.Input", -15.0, 2.68, "IR", (0, 9)),
        ("IR278.Input", -13.4, 2.32, "IR", (0, -15)),
        ("IR279.Input", -23.0, -3.90, "IR", (-4, 9)),
        ("Sham296.Input", 16.0, 1.58, "Sham", (0, 9)),
        ("Sham286.Input", 14.7, 1.32, "Sham", (0, -15)),
        ("Sham295.Input", 17.2, -3.60, "Sham", (0, 9)),
    ]
    fig, ax = plt.subplots(figsize=(W_HALF, 5.55))
    for grp, c in [("IR", UP), ("Sham", DOWN)]:
        g = [p for p in pts if p[3] == grp]
        ax.scatter([p[1] for p in g], [p[2] for p in g], s=110, color=c,
                   label=grp, zorder=3, edgecolors="white", linewidths=0.6)
    for name, x, y, grp, off in pts:
        ax.annotate(name, (x, y), xytext=off, textcoords="offset points",
                    ha="center", va="center", fontsize=14,
                    color=UP if grp == "IR" else DOWN)
    ax.axhline(0, color="#DDDDDD", lw=0.8, zorder=0)
    ax.axvline(0, color="#DDDDDD", lw=0.8, zorder=0)
    ax.set_xlabel("PC1: 95% variance")
    ax.set_ylabel("PC2: 2% variance")
    ax.set_title(" mRNA-seq PCA (IR vs Sham, n=3)")
    ax.set_xlim(-27, 22)
    ax.set_ylim(-4.4, 3.3)
    ax.legend(title="group", frameon=False, loc="upper right",
              fontsize=12, title_fontsize=12, borderpad=0.3,
              labelspacing=0.18, handletextpad=0.4)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    save(fig, os.path.join(OUT, "F3A_PCA.png"))


def _volcano(csv, title, label_down, label_up, out_name):
    df = pd.read_csv(os.path.join(RES, csv)).dropna(
        subset=["padj", "log2FoldChange"]).copy()
    df["nq"] = -np.log10(df["padj"].clip(lower=FLOOR))
    up = (df["log2FoldChange"] > LFC) & (df["padj"] < Q)
    dn = (df["log2FoldChange"] < -LFC) & (df["padj"] < Q)
    ns = ~(up | dn)
    print(out_name, "Up=", int(up.sum()), "Down=", int(dn.sum()))
    fig, ax = plt.subplots(figsize=(W_HALF, 5.55))
    ax.scatter(df.loc[ns, "log2FoldChange"], df.loc[ns, "nq"], s=6,
               alpha=0.5, color="#999999", label="NS", linewidths=0)
    ax.scatter(df.loc[dn, "log2FoldChange"], df.loc[dn, "nq"], s=6,
               alpha=0.6, color=DOWN, label=label_down, linewidths=0)
    ax.scatter(df.loc[up, "log2FoldChange"], df.loc[up, "nq"], s=6,
               alpha=0.6, color=UP, label=label_up, linewidths=0)
    for v in (-LFC, LFC):
        ax.axvline(v, ls="--", lw=0.9, color=GREY)
    ax.axhline(-np.log10(Q), ls="--", lw=0.9, color=GREY)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 (adjusted P)")
    ax.set_title(title)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.14),
              ncol=3, fontsize=14, columnspacing=1.4, handletextpad=0.4,
              markerscale=3)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    save(fig, os.path.join(OUT, out_name))


def panel_venn():
    ac4c = set(pd.read_csv(os.path.join(R1, "all_diff_ac4c_sites.csv"))
               ["GeneName"].dropna().astype(str).str.strip().unique())
    ir_df = pd.read_csv(os.path.join(RES, "ir_mrna_deseq2_full.csv"))
    cold_df = pd.read_csv(os.path.join(RES, "cold_mrna_deseq2_full.csv"))
    ir = set(ir_df[(ir_df["padj"].notna()) & (abs(ir_df["log2FoldChange"]) > 1) &
                   (ir_df["padj"] < 0.05)]["Symbol"].dropna().astype(str).str.strip())
    cold = set(cold_df[(cold_df["padj"].notna()) & (abs(cold_df["log2FoldChange"]) > 1) &
                       (cold_df["padj"] < 0.05)]["Symbol"].dropna().astype(str).str.strip())
    ir_all = set(ir_df["Symbol"].dropna().astype(str).str.strip())
    cold_all = set(cold_df["Symbol"].dropna().astype(str).str.strip())
    univ = ir_all & cold_all
    M = len(univ)
    k_a, k_i, k_c = len(ac4c & univ), len(ir & univ), len(cold & univ)
    ac_ir = len((ac4c & ir) & univ)
    ac_cold = len((ac4c & cold) & univ)
    ir_cold = len((ir & cold) & univ)
    triple = len((ac4c & ir & cold) & univ)

    def hp(K, n, k):
        return float(np.exp(hypergeom.logsf(k - 1, M, K, n)))

    p_ai, p_ic, p_ac, p_t = (hp(k_a, k_i, ac_ir), hp(k_i, k_c, ir_cold),
                             hp(k_a, k_c, ac_cold), hp(ac_ir, k_c, triple))
    print("P:", f"{p_ai:.3e}", f"{p_ic:.3e}", f"{p_ac:.3e}", f"{p_t:.3e}")

    fig, ax = plt.subplots(figsize=(W_HALF, 5.7))
    v = venn3(subsets=(len(ac4c - ir - cold), len(ir - ac4c - cold),
                       len(ac4c & ir - cold), len(cold - ac4c - ir),
                       len(ac4c & cold - ir), len(ir & cold - ac4c), triple),
              set_labels=(f"ac4C-remodeled genes\n(acRIP-seq, n={len(ac4c)})",
                          f"IR transcriptome\n(DESeq2, n={len(ir)})",
                          f"Cold-responsive genes\n(GSE322673, n={len(cold)})"),
              set_colors=("#D73027", "#66A61E", "#4575B4"), alpha=0.6, ax=ax)
    for t in v.set_labels:
        t.set_fontsize(18)
    for t in v.subset_labels:
        if t:
            t.set_fontsize(19)
    ax.set_title("Triple intersection (not to scale)", fontsize=18, fontweight="bold")
    fig.subplots_adjust(top=0.92, bottom=0.02, left=0.02, right=0.98)
    save(fig, os.path.join(OUT, "F3D_triple_venn.png"))


if __name__ == "__main__":
    panel_pca()
    _volcano("ir_mrna_deseq2_full.csv", "mRNA-seq: IR vs Sham (Up=1287, Down=904)",
             "Down", "Up", "F3B_volcano.png")
    _volcano("cold_mrna_deseq2_full.csv", "GSE322673: Cold vs RT (Up=121, Down=362)",
             "Down in Cold", "Up in Cold", "GSE322673_volcano.png")
    panel_venn()
