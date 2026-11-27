#!/usr/bin/env python3
"""
Purpose:
    Build Supplementary Figure 1 (the root counterpart of Figure 2, panels `a–s` with exact
    chromosome fold counts `n`) by invoking `make_fig2_comprehensive.2.py` for `root` and
    exporting `Supplementary_Figure_1_root_as.{pdf,png,_values.tsv,_legend.txt}`.

Manuscript (Peleke et al., 2026):
    Produces Supplementary Figure 1 offline directly from `results/evaluation/`.

Usage & I/O (run from repo root):
    FIG_OUT_DIR=/path/to/out python3 code/04_figures/make_suppfig1_root.2.py

How to change or extend:
    Redirect output via `FIG_OUT_DIR`; panel geometry and species lists are inherited from
    `make_fig2_comprehensive.2.py` and `make_pub_figures_split.7.py`.
"""
import importlib.util
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("fig2comp2", HERE / "make_fig2_comprehensive.2.py")
C = importlib.util.module_from_spec(spec)
spec.loader.exec_module(C)
F = C.F

TISSUE = "root"
SRC = F.OUT / f"Figure2_{TISSUE}_with_clade_MSR_as"
DST = F.OUT / f"Supplementary_Figure_1_{TISSUE}_as"


def write_legend(path, t, nmap):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(
            "Supplementary Figure 1 | Flowering plant deepCRE model performance for single-species (SSR) and "
            "multi-species (MSR) training in root tissue.\n"
            "The root counterpart of Figure 2. Panels, species order, markers and colour scale are identical to "
            "Figure 2, which shows the same analyses for leaf tissue; the panel descriptions there apply here "
            "unchanged, and only the values differ.\n"
            "a, The thirteen flowering plant species and their phylogenetic relationships, as in Figure 2a. The "
            "cladogram gives the species order of all panels; n is the number of chromosome folds (SSR models) per "
            "species shown in b and c.\n"
            "b, c, auROC (b) and accuracy (c) of the SSR deepCRE models for root tissue. Grey dots, one chromosome "
            "fold each on its hold-out chromosome; black dot with printed value, the fold provided on the web "
            "toolkit (highest validation auROC); open circle and horizontal line, mean and standard deviation "
            f"across folds. Number of chromosome folds per species: n = {C.n_sentence(nmap)}.\n"
            "d, Cross-species evaluation of the SSR models: mean predicted probability of high expression for genes "
            "labelled as lowly expressed (left) and highly expressed (right), for every combination of model species "
            "(rows) and target species (columns), as in Figure 2d.\n"
            f"e, UPGMA tree of Euclidean distances between the 52-value cross-prediction profiles of the thirteen "
            f"species (rows and columns of d). Mantel test against the topological distances of the cladogram in a: "
            f"r = {t['mantel_r']:.3f}, p = {t['mantel_p']:.4f} (Pearson, 9,999 permutations, 78 species pairs). "
            f"Leaf order: {', '.join(t['order'])}.\n"
            "f-i, MSR models trained on all 13 species (rows: held-out species): cladogram (f); performance on the "
            "hold-out species using auROC (g) and accuracy (h) (black dot with printed value), in comparison to the "
            "SSR web-tool model of that species (open diamond) and to the cross-species predictions of the SSR models "
            "of the other 12 species (light grey dots). As in Figure 2, the MSR model is scored on all low- and "
            "high-expression genes of the held-out species, the SSR web-tool model on the paralog-filtered genes of "
            "its own held-out chromosome. (i) Cross-species evaluation of the all-species MSR models, shown as in d: "
            "rows are the held-out species of each model, columns the target species. The diagonal is the held-out "
            "species, i.e. the only column of that row the model has not seen in training.\n"
            "j-n, MSR models trained on the 9 eudicot species only (rows: held-out eudicot): cladogram of the "
            "eudicots (j); hold-out species auROC (k) and accuracy (l) as in g and h; cross-species evaluation as in "
            "d for the eudicot-only models within their own clade (m) and on the grasses (n).\n"
            "o-s, MSR models trained on the 4 grass species only (rows: held-out grass), shown as in j-n: cladogram "
            "of the grasses (o); hold-out species auROC (p) and accuracy (q); cross-species evaluation within the "
            "grasses (r) and on the eudicots (s).\n"
            "Source data are provided as Source data files.\n"
        )


def main():
    t_info, nmap = C.main((TISSUE,))[TISSUE]
    for ext in (".pdf", ".png", "_values.tsv"):
        src = Path(str(SRC) + ext)
        assert src.exists(), src
        shutil.copyfile(src, str(DST) + ext)
        print(f"  {DST.name}{ext}")
    write_legend(str(DST) + "_legend.txt", t_info, nmap)
    print(f"  {DST.name}_legend.txt")


if __name__ == "__main__":
    main()
