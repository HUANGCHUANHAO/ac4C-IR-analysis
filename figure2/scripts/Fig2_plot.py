# Inputs : figure2/results tables built by Fig2_data.py and Fig2_data_kegg.R
# Outputs: report/Fig2A_genomic_distribution.png, Fig2B_5mer_enrichment.png,
#          Fig2C_CCW_density.png, Fig2D_KEGG.png
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "common")))
from plot_style import UP, DOWN, W_HALF, save, fmt_pow, dotplot_enrichment

FIG = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RES = os.path.join(FIG, "results")
OUT = os.path.join(FIG, "report")
os.makedirs(OUT, exist_ok=True)


def panel_region():
    reg = pd.read_csv(os.path.join(RES, "Fig2A_region_proportions.csv"))
    reg = reg.rename(columns={reg.columns[0]: "region"})
    regions = list(reg["region"])
    x = np.arange(len(regions))
    w = 0.38
    fig, ax = plt.subplots(figsize=(W_HALF, 5.0))
    b1 = ax.bar(x - w/2, reg["up"], w, color=UP, label="Up-acetylated peaks (IR)")
    b2 = ax.bar(x + w/2, reg["down"], w, color=DOWN, label="Down-acetylated peaks (IR)")
    for b in (b1, b2):
        for r_ in b:
            ax.text(r_.get_x() + r_.get_width()/2, r_.get_height() + 0.8,
                    f"{r_.get_height():.1f}", ha="center", va="bottom", fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(regions, rotation=30, ha="right")
    ax.set_ylabel("Proportion (%)")
    ax.set_xlabel("Genomic region")
    ax.set_ylim(0, 92)
    ax.set_title("Genomic distribution of differential ac4C peaks")
    ax.legend(frameon=False, loc="upper right")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    save(fig, os.path.join(OUT, "Fig2A_genomic_distribution.png"))


def panel_5mer():
    mer = pd.read_csv(os.path.join(RES, "Fig2B_5mer_enrich.csv"))
    mer = mer[mer["padj"] < 0.05]
    top = mer.sort_values("log2OR", ascending=False).head(10)
    bot = mer.sort_values("log2OR").head(10)
    sel = pd.concat([top, bot]).sort_values("log2OR")
    colors = [UP if v > 0 else DOWN for v in sel["log2OR"]]
    fig, ax = plt.subplots(figsize=(W_HALF, 5.0))
    ax.barh(sel["mer"], sel["log2OR"], color=colors, height=0.7)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("log2 odds ratio")
    ax.set_title("Top enriched 5-mers: up- vs down-acetylated peaks (all peaks)")
    ax.tick_params(axis="y", labelsize=14)
    ax.legend(handles=[mpatches.Patch(color=UP, label="Up-enriched"),
                       mpatches.Patch(color=DOWN, label="Down-enriched")],
              loc="lower right", framealpha=0.9)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    save(fig, os.path.join(OUT, "Fig2B_5mer_enrichment.png"))


def panel_ccw():
    plt.rcParams.update({
        "axes.titlesize": 17, "axes.labelsize": 15,
        "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 13})
    ccw = pd.read_csv(os.path.join(RES, "Fig2C_ccw_density.csv"))
    p_up, p_dn = 1.260e-82, 7.016e-74
    pos = [0.7, 2.3]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.9), sharey=True)
    for ax, direction, p, ttl in [
        (axes[0], "up", p_up, "Up-acetylated peaks"),
        (axes[1], "down", p_dn, "Down-acetylated peaks")]:
        sub = ccw[ccw["direction"] == direction]
        data = [sub["ccw_density_real"], sub["ccw_density_shuffled"]]
        bp = ax.boxplot(data, positions=pos, widths=0.5, patch_artist=True,
                        showfliers=True, flierprops=dict(marker="o", ms=3, alpha=0.6),
                        medianprops=dict(color="black"))
        bp["boxes"][0].set_facecolor(UP if direction == "up" else DOWN)
        bp["boxes"][1].set_facecolor("#BFBFBF")
        ax.set_xticks(pos)
        ax.set_xticklabels(["Real peaks", "Shuffled background"])
        ax.set_title(f"{ttl}\nWilcoxon P={fmt_pow(p)}", pad=8)
        ax.set_ylim(0, 20)
        ax.set_xlim(0, 3)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("CCW density (per 100 bp)")
    fig.subplots_adjust(left=0.075, right=0.99, bottom=0.14, top=0.84, wspace=0.12)
    save(fig, os.path.join(OUT, "Fig2C_CCW_density.png"))


def panel_kegg():
    dotplot_enrichment(RES, "Fig2D_kegg.csv",
                       os.path.join(OUT, "Fig2D_KEGG.png"),
                       "KEGG enrichment \u2013 differential ac4C genes", 20, 6.0)


if __name__ == "__main__":
    panel_region()
    panel_5mer()
    panel_ccw()
    panel_kegg()
