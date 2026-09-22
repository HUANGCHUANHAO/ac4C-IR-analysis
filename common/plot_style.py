"""Shared plotting style imported by the FigXX_plot.py panel scripts.
All single panels are redrawn from the result CSVs each data script writes.
Arial; unified font sizes; 300 dpi; scientific P values are rendered as inline
'coefficient x 10^exponent' (caret, normal-size exponent), matching the
manuscript text (e.g. 1.26x10^-82), never 1e-82 or mathtext superscripts.
Standard export widths at 300 dpi for a 4200 px composite inner width:
  full row = 4200 px (14.0 in), half column = 2080 px (~6.93 in),
  one-third column = 1372 px (~4.57 in).
"""
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

UP = "#D73027"
DOWN = "#4575B4"
GREY = "#666666"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "axes.unicode_minus": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.titlesize": 19,
    "axes.titleweight": "bold",
    "axes.labelsize": 17,
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "legend.fontsize": 15,
    "axes.linewidth": 1.0,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

DPI = 300
W_HALF = 2080 / DPI     # half-column panel width, inches
W_THIRD = 1372 / DPI    # one-third column
W_FULL = 4200 / DPI     # full-width row


def fmt_pow(x, digits=2):
    """Inline scientific notation: 1.26e-82 -> '1.26 x 10^-82'."""
    try:
        x = float(x)
    except Exception:
        return str(x)
    if x == 0:
        return "0"
    if x >= 1e-3:
        return f"{x:.3f}".rstrip("0").rstrip(".")
    e = math.floor(math.log10(abs(x)))
    m = x / 10.0 ** e
    s = f"{m:.{digits}f}".rstrip("0").rstrip(".")
    return f"{s}×10^{e}"


def log10_pow_breaks(vmin, vmax, every=5):
    """Round-decade breaks + inline labels covering [vmin, vmax]."""
    lo = math.floor(math.log10(vmin))
    hi = math.ceil(math.log10(vmax))
    vals, labs = [], []
    for e in range(lo, hi + 1):
        if (e - lo) % every == 0 or e == hi:
            v = 10.0 ** e
            vals.append(v)
            labs.append(f"1×10^{e}")
    return vals, labs


def save(fig, path):
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("saved", path)


# --- helper for the enrichment bubble plots (Fig2D, FigS1A, FigS1B) ---
def dotplot_enrichment(RES, csv, out, title, topn=15, h=5.4):
    import os, math
    import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib.colors as mcolors
    r = pd.read_csv(os.path.join(RES, csv))
    r = r.sort_values("p.adjust").head(topn).copy()
    r["gr"] = r["GeneRatio"].map(
        lambda s: float(s.split("/")[0]) / float(s.split("/")[1]))
    r["sig"] = -np.log10(r["p.adjust"])
    r["Description"] = r["Description"].str.replace(
        r" - Mus musculus \(house mouse\)", "", regex=True)
    r = r.sort_values("p.adjust", ascending=False)
    y = np.arange(len(r))
    vlo = int(np.floor(r["sig"].min())); vhi = int(np.ceil(r["sig"].max()))
    norm = mcolors.Normalize(vmin=vlo, vmax=vhi)
    smin, smax = r["Count"].min(), r["Count"].max()
    def pt_sz(c):
        return 90 + (np.clip(c, smin, smax) - smin) / max(smax - smin, 1) * 320
    fig, ax = plt.subplots(figsize=(W_HALF, h))
    sc = ax.scatter(r["gr"], y, s=pt_sz(r["Count"]), c=r["sig"], cmap="RdBu_r",
                    norm=norm, edgecolors="black", linewidths=0.45, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(r["Description"], fontsize=15)
    ax.set_xlabel("GeneRatio"); ax.set_title(title)
    ax.grid(axis="x", alpha=0.3, zorder=0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    cb = fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.03)
    cb.set_label("p.adjust", fontsize=15)
    step = max(2, math.ceil((vhi - vlo) / 6))
    tv = list(range(vlo, vhi + 1, step))
    if vhi not in tv:
        tv.append(vhi)
    cb.set_ticks(tv); cb.set_ticklabels([f"1\u00d710^-{t}" for t in tv])
    qs = sorted(set([int(round(smin)), int(round((smin + smax) / 2)), int(round(smax))]))
    handles = [plt.scatter([], [], s=pt_sz(q), c="#E8E8E8", edgecolors="black",
                           linewidths=0.45) for q in qs]
    ax.legend(handles, [str(q) for q in qs], title="Count", framealpha=0.9,
              loc="lower right", labelspacing=1.3, handletextpad=1.1, fontsize=14)
    save(fig, out)
