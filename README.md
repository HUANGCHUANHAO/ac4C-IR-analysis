# ac4C–IR

Analysis code for a study on ac4C RNA modification and cold exposure in ischemic
myocardial injury. The pipeline starts from acRIP-seq and mRNA-seq of mouse
myocardium after ischemia–reperfusion, intersects the differential ac4C sites and
differential genes with a public cold-exposure transcriptome, and narrows the
result to an S100a8/a9–Cxcr2 neutrophil signature, which is then checked in
single-cell data and in two independent human cohorts.

The repository holds the analysis code and the one annotation file no script can
run without, the mm10 GTF under `figure2/raw_data/annotation/`. Input data,
intermediate tables and figures are
not included. Every `raw_data/` folder carries a `DOWNLOAD.txt` naming the files
its scripts expect and where to get them, and each script writes its output into
a `results/` or `report/` folder that it creates on the first run.

## Scope

This is the code behind the accompanying manuscript, deposited so that the
analysis can be read and rerun. It is not a maintained software package, so use
the package versions listed under Environment below; later releases may not
reproduce the published tables exactly.

## Repository layout

| Path | Contents |
|---|---|
| `common/` | `plot_style.py`, the colours, sizes and helpers that every plotting script imports |
| `figure1/` ... `figure6/` | One folder per figure, each with `scripts/` and usually `raw_data/`. Those `raw_data/` folders hold nothing but a `DOWNLOAD.txt`, apart from the mm10 GTF bundled in `figure2/raw_data/annotation/`. `figure6/` has no `raw_data/` of its own; `Fig6.R` reads the two series matrices from `figure4/raw_data/` |
| `supplementary/` | Scripts for Figure S1 |
| `requirements.txt` | Python packages used by the Python scripts |
| `R_REQUIREMENTS.txt` | R packages used by the R scripts |
| `LICENSE` | MIT |

The plotting scripts in `figure1/`, `figure2/`, `figure3/`, `figure5/` and
`supplementary/` work out the path to `common/` from their own file, add it to
`sys.path` and then import `plot_style`, so that folder has to stay at the
repository root. Every other path
a script needs is worked out from the script's own location as well (the R
scripts read theirs from `--file=`). Nothing depends on the working directory, so
the scripts can be run from anywhere and the repository can be cloned and renamed
freely, as long as files are not moved between folders inside it.

## Data sources

| Dataset | Material | Groups | Used for | Figure |
|---|---|---|---|---|
| In-house acRIP-seq | Mouse myocardium | Sham vs IR 45 min + 24 h, n = 3 each | Differential ac4C sites | 1, 2 |
| In-house mRNA-seq | Mouse myocardium | Sham vs IR 45 min + 24 h, n = 3 each | Differential expression | 3 |
| GSE322673 | Mouse myocardium | Cold-MI 3 / RT-MI 3, 16 °C vs 22 °C for 7 d | Cold-exposure contrast (exploratory) | 3 |
| GSE247139 | Mouse heart, single cell | Sham 11196 / IR6h 7881 cells | Cell-type localisation | 5 |
| GSE66360 | Human CD146+ circulating endothelial cells | AMI 49 / control 50 | Single-gene ROC, model training | 4, 6 |
| GSE48060 | Human whole blood | 31 patients after a first-time MI / 21 controls | Single-gene ROC, external validation | 4, 6 |
| UCSC mm10 (genome FASTA) and mm10 GTF | Mouse genome | – | Peak annotation, 5-mer and motif analysis | 2 |

The in-house acRIP-seq and mRNA-seq data are not deposited in a public repository.
Their processed result tables accompany the manuscript as Supplementary Tables S1
(differential ac4C sites) and S2 (differential expression), and the scripts here
read those two supplementary tables directly. The raw libraries and the
vendor-supplied count and site tables are available from the authors on request.

Nothing in the table above is stored here apart from the bundled mm10 GTF. Put the
files where the scripts expect them, using each folder's `DOWNLOAD.txt` as the
checklist.

- `figure1/raw_data/` gets `Supplementary_Table_S1_differential_ac4C_sites.xlsx`,
  the in-house ac4C site table (Supplementary Table S1)
- `figure2/raw_data/MM10/` gets the mm10 genome FASTA and its `samtools faidx`
  index. The annotation GTF in `annotation/` ships with the repository.
- `figure3/raw_data/` gets `Supplementary_Table_S2_differential_expression.xlsx`
  (Supplementary Table S2, the in-house mRNA-seq counts) and
  `GSE322673/GSE322673_raw_counts.csv`, which is public
- `figure4/raw_data/` gets `GSE66360_series_matrix.txt` and
  `GSE48060_series_matrix.txt`, both public. `figure4/scripts/Fig4.R` and
  `figure6/scripts/Fig6.R` both read them.
- `figure5/raw_data/GSE247139/` gets the 10x matrices for the 0 h and 6 h samples,
  which are public

## Environment

We used Python 3.14 and R 4.4.2.

```
python -m pip install -r requirements.txt
```

For R, install what `R_REQUIREMENTS.txt` lists, from CRAN and Bioconductor. The R
scripts are meant to be run with `Rscript`. They read their own path from the
command line to work out where the repository root is, which does not happen if
they are `source()`-ed in an interactive session.

## Running the scripts

Run these from the repository root. Each script creates its figure's `results/`
and `report/` folders on the first run. Some figures read tables produced by an
earlier one, so the order below is not arbitrary.

1. Figure 1: differential ac4C sites
   - `Rscript figure1/scripts/Fig1_data.R` → `figure1/results/all_diff_ac4c_sites.csv`, `figure1/results/Fig1C_chrom_distribution.csv`
   - `python figure1/scripts/Fig1_plot.py` → `figure1/report/Fig1A_volcano.png`, `Fig1B_peak_length_distribution.png`, `Fig1C_chromosomal_distribution.png`

2. Figure 2: peak annotation, 5-mer, CCW density, KEGG
   - `python figure2/scripts/Fig2_data.py` → `figure2/results/Fig2A_region_counts.csv`, `Fig2A_region_proportions.csv`, `Fig2A_region_prop_long.csv`, `Fig2A_region_summary.txt`, `Fig2B_5mer_enrich.csv`, `Fig2C_ccw_density.csv`, `Fig2C_ccw_summary.txt`, `peak_sequences.csv`
   - `Rscript figure2/scripts/Fig2_data_kegg.R` → `figure2/results/Fig2D_kegg.csv`
   - `python figure2/scripts/Fig2_plot.py` → `figure2/report/Fig2A_genomic_distribution.png`, `Fig2B_5mer_enrichment.png`, `Fig2C_CCW_density.png`, `Fig2D_KEGG.png`

3. Figure 3: mRNA-seq differential expression and the three-way intersection
   - `Rscript figure3/scripts/Fig3_data_deseq2.R` → `figure3/results/Fig3A_pca.csv`, `Fig3B_ir_deg.csv`, `Fig3B_ir_summary.txt`, `Fig3C_cold_deg.csv`, `Fig3C_cold_summary.txt`, `ir_mrna_deseq2_full.csv`, `cold_mrna_deseq2_full.csv`
   - `python figure3/scripts/Fig3_data_intersection.py` → `figure3/results/Fig3D_triple_genes.csv`, `Fig3D_core_genes.csv`, `Fig3D_intersection_summary.txt`
   - `python figure3/scripts/Fig3_plot.py` → `figure3/report/F3A_PCA.png`, `F3B_volcano.png`, `F3D_triple_venn.png`, `GSE322673_volcano.png`

4. Figure 4: single-gene ROC in the human cohorts
   - `Rscript figure4/scripts/Fig4.R` → `figure4/results/Fig4A_core_genes.csv`, `Fig4B_roc_gse66360.csv`, `Fig4C_roc_gse48060.csv`, `Fig4D_boxstats_gse66360.csv`, `Fig4E_boxstats_gse48060.csv` and the matching panels in `figure4/report/`

5. Figure 6: combined model, nomogram, calibration, DCA
   - `Rscript figure6/scripts/Fig6.R` → `figure6/results/Fig6A_auc.csv`, `Fig6A_cutpoint.csv`, `Fig6C_calibration.txt`, `Fig6D_dca_train.csv`, `Fig6D_dca_valid.csv`, `Table1_model_variants.csv`, `external_probabilities.csv`, `lrm_model.txt`, `model_coefficients.csv` and the panels in `figure6/report/`

6. Figure 5: single-cell localisation (independent of figures 1–4; needs only `figure5/raw_data/GSE247139/`)
   - `python figure5/scripts/Fig5_data.py` → `figure5/results/GSE247139_processed.h5ad`, `figure5/results/Fig5B_abundance.csv`, `figure5/results/Fig5C_expr_by_celltype.csv`, `figure5/results/run_summary.txt`
   - `python figure5/scripts/Fig5_plot.py` → `figure5/report/Fig5A_UMAP.png`, `Fig5B_celltype_abundance.png`, `Fig5C_core_genes_dotplot.png`, `Fig5D_myeloid_violin.png`

7. Figure S1: GO BP enrichment of up- and down-acetylated genes (needs the Figure 1 output)
   - `Rscript supplementary/FigS1_GO_BP_split.R` → `supplementary/results/SourceData_FigS1A_GO_BP_up.csv`, `supplementary/results/SourceData_FigS1B_GO_BP_down.csv`
   - `python supplementary/FigS1A_GO_up_plot.py` → `supplementary/FigS1A_GO_BP_up.png`
   - `python supplementary/FigS1B_GO_down_plot.py` → `supplementary/FigS1B_GO_BP_down.png`

`figure5/scripts/Fig5_data.py` is the slow step and writes a ~600 MB
`figure5/results/GSE247139_processed.h5ad` that the plotting script re-reads.
Everything else finishes in a few minutes per script.

## Notes

- The in-house bulk libraries were collected 24 h after reperfusion, whereas
  GSE247139 was collected 6 h after IR. The manuscript says so as well. The two
  datasets are not compared quantitatively.
- GSE48060 cases are patients sampled after a first-time myocardial infarction.
- `figure6/scripts/Fig6.R` sets `set.seed(20260905)` before the bootstrap and
  cross-validation steps, so the model figures are reproducible.
- The R/Python split is deliberate. Differential expression, enrichment, ROC and
  the logistic model are in R, while the peak-level and single-cell work is in
  Python.
- `figure2/scripts/Fig2_data_kegg.R` queries the KEGG REST API at
  `rest.kegg.jp`, so it needs a network connection. The pathway set it gets back
  is a live database: as KEGG is updated, the gene universe behind the background
  changes, and the p and adjusted p values in `Fig2D_kegg.csv` can shift in their
  fifth significant digit. The pathways reported, their order and the gene counts
  stay the same. Figures regenerated from a later KEGG release will therefore
  differ from the published `Fig2D_kegg.csv` at that level.

## License

MIT. See `LICENSE`.
