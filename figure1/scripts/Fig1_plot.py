# Inputs : results/all_diff_ac4c_sites.csv, results/Fig1C_chrom_distribution.csv (built by Fig1_data.R)
# Outputs: report/Fig1A_volcano.png, Fig1B_peak_length_distribution.png, Fig1C_chromosomal_distribution.png
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "common")))
from plot_style import UP, DOWN, GREY, W_HALF, W_FULL, save, fmt_pow

FIG = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RES = os.path.join(FIG, "results")
OUT = os.path.join(FIG, "report")
os.makedirs(OUT, exist_ok=True)


def panel_volcano(df):
    d = df.copy()
    d["log2FC"] = np.log2(d["Foldchange"]) * np.where(d["direction"] == "up", 1, -1)
    d["neglogP"] = -np.log10(d["P_value"])
    fig, ax = plt.subplots(figsize=(W_HALF, 5.55))
    for direction, c, lab in [("down", DOWN, "Down"), ("up", UP, "Up")]:
        s = d[d["direction"] == direction]
        ax.scatter(s["log2FC"], s["neglogP"], s=6, alpha=0.55, color=c,
                   label=lab, linewidths=0)
    ax.axvline(-1, ls="--", lw=0.9, color=GREY)
    ax.axvline(1, ls="--", lw=0.9, color=GREY)
    ax.axhline(-np.log10(1e-5), ls="--", lw=0.9, color=GREY)
    top8 = d.sort_values("Foldchange", ascending=False).head(8).copy()
    top8 = top8.sort_values("neglogP").reset_index(drop=True)
    lab_y, last = [], -1e9
    for _, r in top8.iterrows():
        yy = max(r["neglogP"], last + 0.9)
        lab_y.append(yy)
        last = yy
    top8["lab_y"] = lab_y
    xmax = top8["log2FC"].max()
    for _, r in top8.iterrows():
        ax.annotate(r["GeneName"], (r["log2FC"], r["neglogP"]),
                    xytext=(xmax + 0.55, r["lab_y"]), fontsize=14,
                    ha="left", va="center",
                    arrowprops=dict(arrowstyle="-", lw=0.4, color=GREY))
    ax.set_xlim(-9.5, xmax + 2.4)
    ax.set_xlabel("log2 fold change (IR vs Sham)")
    ax.set_ylabel("-log10 P value")
    ax.set_title(f"Differential ac4C sites (diffreps, FC\u22652, P<{fmt_pow(1e-5)})")
    ax.legend(frameon=False, loc="upper left")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    save(fig, os.path.join(OUT, "Fig1A_volcano.png"))


def panel_peak_length(df):
    fig, ax = plt.subplots(figsize=(W_HALF, 4.7))
    bins = np.linspace(0, df["Peak_length"].quantile(0.995), 60)
    for direction, c in [("up", UP), ("down", DOWN)]:
        v = df.loc[df["direction"] == direction, "Peak_length"]
        ax.hist(v, bins=bins, density=True, alpha=0.55, color=c,
                label=f"{direction} (median {int(round(v.median()))} bp)")
    ax.set_xlabel("Peak length (bp)")
    ax.set_ylabel("Density")
    ax.set_title("Peak length distribution")
    ax.legend(frameon=False)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    save(fig, os.path.join(OUT, "Fig1B_peak_length_distribution.png"))


def panel_chromosome():
    cc = pd.read_csv(os.path.join(RES, "Fig1C_chrom_distribution.csv"))

    def chr_key(name):
        s = name.replace("chr", "")
        if s == "X":
            return 23
        if s == "Y":
            return 24
        return int(s)

    cc["k"] = cc["chrom"].map(chr_key)
    cc = cc.sort_values("k")
    fig, ax = plt.subplots(figsize=(W_FULL, 3.95))
    x = np.arange(len(cc))
    ax.bar(x, cc["up"], color=UP, width=0.62, label="up")
    ax.bar(x, cc["down"], bottom=cc["up"], color=DOWN, width=0.62, label="down")
    ax.set_xticks(x)
    ax.set_xticklabels(cc["chrom"], rotation=45, ha="right", fontsize=15)
    ax.set_ylabel("Differential ac4C sites")
    ax.set_title("Chromosomal distribution of differential ac4C sites")
    ax.legend(frameon=False, ncol=2, loc="upper right")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    save(fig, os.path.join(OUT, "Fig1C_chromosomal_distribution.png"))


if __name__ == "__main__":
    sites = pd.read_csv(os.path.join(RES, "all_diff_ac4c_sites.csv"))
    panel_volcano(sites)
    panel_peak_length(sites)
    panel_chromosome()
