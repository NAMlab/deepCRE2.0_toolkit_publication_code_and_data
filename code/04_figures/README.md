# Main and Root/Cross-Tissue Publication Figures (`code/04_figures/`)

Python scripts (Python 3.11 + `requirements_analysis.txt`) generating Figure 1, Figure 2 (leaf), Supplementary Figure 1 (root), Supplementary Figure 2 (cross-tissue evaluation), and the Figure 2 / Supplementary Figures 1–2 Source Data Excel workbook.

> **Offline reproduction:** Everything except Figure 1 (`build_updated_fig1.py` and `build_figure1_setup()`) runs offline directly from [`results/evaluation/`](../../results/README.md). Figure 1 reads UI screenshot/crop PNGs under `EXTERNAL/figure_inputs/` (`fig1_ab_crop.png`, `deepcre_random100_*.png`, `tutorial_canvas_10.png`) that are not distributed in this repository. (Figure 3 / Supplementary Figure 3 benchmark against PhytoExpr and PlantCAD2 is maintained outside this repository.)

---

## Scripts

| Script | Purpose | Manuscript item | CLI Usage |
|---|---|---|---|
| `build_updated_fig1.py` | Renders the annotated Figure 1a web-tool workflow + cleaned Figure 1b schematic (`No. of genes`, `log10(mean TPM+1)`, `gDR`) from external `EXTERNAL/figure_inputs/` screenshots and crops. | **Figure 1** | `python3 code/04_figures/build_updated_fig1.py` *(requires `EXTERNAL/figure_inputs/*.png`)* |
| `make_pub_figures_split.7.py` | Shared data loader, UPGMA/Mantel computation, cladogram/dot/heatmap drawing primitives, and base figure splits (`h-q` + `SuppFig_MSR`). Imported at runtime via `importlib` by all scripts below. | Shared plotting & stats engine | `FIG_OUT_DIR=/path/to/out python3 code/04_figures/make_pub_figures_split.7.py` |
| `make_fig2_comprehensive.0.py` | Layout engine (`draw()`) arranging SSR, MSR-all, MSR-eudicots, and MSR-grasses panels (`a–q`) with shared color scale. | Imported by `.2.py` | `FIG_OUT_DIR=/path/to/out python3 code/04_figures/make_fig2_comprehensive.0.py [leaf\|root]` |
| `make_fig2_comprehensive.2.py` | Final **Figure 2** builder (`a–s`): adds panel `f` on MSR-all cladogram, prints chromosome fold counts `n` beside panel `a`, and standardizes `Dcar` heatmap labels. | **Figure 2** (leaf `.pdf`/`.png`, `_values.tsv`, `_legend.txt`) | `FIG_OUT_DIR=/path/to/out python3 code/04_figures/make_fig2_comprehensive.2.py [leaf\|root]` |
| `make_suppfig1_root.2.py` | Publishes the root counterpart built by `make_fig2_comprehensive.2.py` as `Supplementary_Figure_1_root_as.*` with its root-specific legend. | **Supplementary Figure 1** | `FIG_OUT_DIR=/path/to/out python3 code/04_figures/make_suppfig1_root.2.py` |
| `make_suppfig2_tissue.7.py` | Plots leaf vs. root SSR and MSR performance, 4-category class composition (`stable`, `swapped`, `high/low to intermediate`), probability contrast on swapped genes, and auROC on swapped genes (`a–o`), printing exact fold counts `n`. | **Supplementary Figure 2** (`SuppFig_tissue_v7.*`) | `FIG_OUT_DIR=/path/to/out python3 code/04_figures/make_suppfig2_tissue.7.py` |
| `build_source_data_fig2.5.py` | Exports Source Data Excel sheets (`a–s` for Figure 2 / Supplementary Figure 1 across leaf & root, plus Supplementary Figure 2 panels from `<suppfig2_values.tsv>`), reusing UPGMA/Mantel functions from `make_pub_figures_split.7.py`. | **Source Data** (Figure 2, Supp. Figs. 1–2) | `python3 code/04_figures/build_source_data_fig2.5.py results/evaluation code/04_figures <out/SuppFig_tissue_v7_values.tsv> <out/Source_Data.xlsx>` |

---

## How to Change or Extend

* **Redirect output without editing code:** Set `FIG_OUT_DIR=/your/output/dir` (defaults to `EXTERNAL/figure_output/figures_working_20261015/figures`).
* **Add a species or reorder phylogeny:** Update `ORDER`, `LABELS`, `TOPO_ALL`, `TOPO_DIC`, `TOPO_GRA`, and `CLADES` in `make_pub_figures_split.7.py`; downstream `.2.py`, `make_suppfig1_root.2.py`, `make_suppfig2_tissue.7.py`, and `build_source_data_fig2.5.py` inherit those definitions automatically.
