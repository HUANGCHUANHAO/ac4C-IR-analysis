# FigS1A GO BP enrichment dot-plot, up-acetylated genes.
# Input : supplementary/results/SourceData_FigS1A_GO_BP_up.csv (built by FigS1_GO_BP_split.R)
# Output: supplementary/FigS1A_GO_BP_up.png
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "common"))
from plot_style import dotplot_enrichment

RES = os.path.join(ROOT, "supplementary/results")
OUT = HERE
os.makedirs(OUT, exist_ok=True)
dotplot_enrichment(RES, "SourceData_FigS1A_GO_BP_up.csv",
                   os.path.join(OUT, "FigS1A_GO_BP_up.png"),
                   "GO BP enrichment – up-acetylated genes", 15, 5.4)
