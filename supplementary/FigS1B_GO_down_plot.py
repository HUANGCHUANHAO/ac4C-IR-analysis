# FigS1B GO BP enrichment dot-plot, down-acetylated genes.
# Input : supplementary/results/SourceData_FigS1B_GO_BP_down.csv (built by FigS1_GO_BP_split.R)
# Output: supplementary/FigS1B_GO_BP_down.png
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "common"))
from plot_style import dotplot_enrichment

RES = os.path.join(ROOT, "supplementary/results")
OUT = HERE
os.makedirs(OUT, exist_ok=True)
dotplot_enrichment(RES, "SourceData_FigS1B_GO_BP_down.csv",
                   os.path.join(OUT, "FigS1B_GO_BP_down.png"),
                   "GO BP enrichment – down-acetylated genes", 15, 5.4)
