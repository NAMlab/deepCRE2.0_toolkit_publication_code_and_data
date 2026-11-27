#!/usr/bin/env python3
"""
Purpose:
    Layout engine (`draw()`) assembling the multi-block Figure 2 architecture (SSR, MSR-all,
    MSR-eudicots, MSR-grasses) using data loaders, UPGMA/Mantel computation, and plotting
    primitives imported from `make_pub_figures_split.7.py`.

Manuscript (Peleke et al., 2026):
    Imported at runtime via `importlib` by `make_fig2_comprehensive.2.py` to draw Figure 2 and
    Supplementary Figure 1 without duplicating statistical or plotting logic.

Usage & I/O (run from repo root):
    FIG_OUT_DIR=/path/to/out python3 code/04_figures/make_fig2_comprehensive.0.py [leaf|root]
    - Reads `results/evaluation/final/` and `results/evaluation/msr/`
    - Writes `Figure2_<tissue>_with_clade_MSR.{pdf,png,_values.tsv,_legend.txt}` into `F.OUT`

How to change or extend:
    Use `make_fig2_comprehensive.2.py` for the final manuscript layout (`a–s` + fold counts `n`).
    Override output directory via environment variable `FIG_OUT_DIR`.
"""
import importlib.util
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("figs7", HERE / "make_pub_figures_split.7.py")
F = importlib.util.module_from_spec(spec)
spec.loader.exec_module(F)

MM, CELL, FIG_W = F.MM, F.CELL, F.FIG_W
ORDER, DARK, MID, LIGHT = F.ORDER, F.DARK, F.MID, F.LIGHT
FS, FS_S, FS_L = F.FS, F.FS_S, F.FS_L

X_TREE, X_AUROC, X_ACC = (0.0, 24.5), (26.0, 18.0), (51.0, 18.0)
X_HEAT0 = 76.0
X_PT = 148.0                      # performance-similarity tree; wider than v7 because the colour bar moved out
HEAD_A, FOOT_A = 11.0, 27.0       # foot of the SSR block holds brackets, colour bar and the b/c marker key
HEAD, FOOT_FLAT, FOOT_HEAT = 10.0, 4.0, 16.5   # 16.5: the "eudicots" bracket label ends 15.5 mm below the heatmap
GAP = 2.0
BRACKET_Y, CB_Y = 15.5, 19.5      # mm below the heatmap bottom

L_SSR = [
    (dict(marker="o", ms=3.0, mfc=DARK, mec=DARK), "SSR web-tool model (best fold; value)"),
    (dict(marker="o", ms=1.7, mfc=MID, mec=MID), "SSR chromosome folds"),
    (dict(meansd=True), "mean ± s.d. of folds"),
]
L_MSR = [
    (dict(marker="o", ms=3.2, mfc=DARK, mec=DARK), "MSR model on its held-out species (value)"),
    (dict(marker="D", ms=2.6, mfc="white", mec=DARK, mew=0.6), "SSR web-tool model of this species"),
    (dict(marker="o", ms=1.8, mfc=LIGHT, mec=MID, mew=0.2), "SSR models of the other species (cross-species)"),
]


def main(tissues=("leaf", "root")):
    F.setup_rc()
    pf_all = F.read_tsv_dicts(F.FINAL / "S2v8_per_fold.tsv")
    pf = [r for r in pf_all
          if r["status"] in ("converged", "converged_rescue") and r["code"] in ORDER
          and r["auROC_v8"] not in ("", "nan", "NA")]
    D, G = ORDER[:9], ORDER[9:]

    for tissue in tissues:
        full_rows = F.read_tsv_dicts(F.MSR_DIR / f"msr_eval_{tissue}_full.tsv")
        summary, oth_rows = F.compute_msr_summary(full_rows, pf_all, tissue)
        e0, _ = F.read_matrix_tsv(F.FINAL / f"crossspecies_mean_{tissue}_class0.tsv")
        e1, _ = F.read_matrix_tsv(F.FINAL / f"crossspecies_mean_{tissue}_class1.tsv")
        M = {}
        for s in ("dicots", "monocots"):
            for k in (0, 1):
                M[(s, k)], _ = F.read_matrix_tsv(F.MSR_DIR / f"suppfig_msr_mean_{tissue}_{s}_class{k}.tsv")

        vals = []
        blocks = [
            dict(name="Single-species models (SSR)",
                 topo=F.TOPO_ALL, rows=ORDER, head=HEAD_A, foot=FOOT_A,
                 letters=["a", "b", "c", "d"], extra_letter="e", heat_letters="first_only",
                 heats=[(e0, ORDER, "True lowly expressed genes", "Cross-species evaluation (row: model, column: target)"),
                        (e1, ORDER, "True highly expressed genes", None)],
                 auroc=lambda ax: F.ssr_dots(ax, pf, tissue, ORDER, "auROC_v8", "auROC", vals, "b"),
                 acc=lambda ax: F.ssr_dots(ax, pf, tissue, ORDER, "accuracy_v8", "Accuracy", vals, "c"),
                 extra=lambda ax: F.draw_perf_tree(ax, e0, e1),
                 legend_under_b=L_SSR, colourbar=True),
            dict(name="Multi-species models (MSR) — all 13 species (row: held-out species)",
                 topo=F.TOPO_ALL, rows=ORDER, head=HEAD, foot=FOOT_FLAT,
                 letters=["", "f", "g"], heats=None,
                 auroc=lambda ax: F.msr_dots(ax, summary, oth_rows, ORDER, "all", None, "auROC", "auROC", vals, "f"),
                 acc=lambda ax: F.msr_dots(ax, summary, oth_rows, ORDER, "all", None, "accuracy", "Accuracy", vals, "g"),
                 legend_right=(L_MSR, X_HEAT0)),
            dict(name="Multi-species models (MSR) — eudicots only (row: held-out eudicot)",
                 topo=F.TOPO_DIC, rows=D, head=HEAD, foot=FOOT_HEAT,
                 letters=["h", "i", "j", "k", "l"], heat_letters="pair_first",
                 heats=[(M[("dicots", 0)], D, "True low", "on eudicots (own clade)"),
                        (M[("dicots", 1)], D, "True high", None),
                        (M[("dicots", 0)], G, "True low", "on grasses (other clade; cross-species)"),
                        (M[("dicots", 1)], G, "True high", None)],
                 auroc=lambda ax: F.msr_dots(ax, summary, oth_rows, D, "dicots", None, "auROC", "auROC", vals, "i"),
                 acc=lambda ax: F.msr_dots(ax, summary, oth_rows, D, "dicots", None, "accuracy", "Accuracy", vals, "j")),
            dict(name="Multi-species models (MSR) — grasses only (row: held-out grass)",
                 topo=F.TOPO_GRA, rows=G, head=HEAD, foot=FOOT_HEAT,
                 letters=["m", "n", "o", "p", "q"], heat_letters="pair_first",
                 heats=[(M[("monocots", 0)], G, "True low", "on grasses (own clade)"),
                        (M[("monocots", 1)], G, "True high", None),
                        (M[("monocots", 0)], D, "True low", "on eudicots (other clade; cross-species)"),
                        (M[("monocots", 1)], D, "True high", None)],
                 auroc=lambda ax: F.msr_dots(ax, summary, oth_rows, G, "monocots", None, "auROC", "auROC", vals, "n"),
                 acc=lambda ax: F.msr_dots(ax, summary, oth_rows, G, "monocots", None, "accuracy", "Accuracy", vals, "o")),
        ]
        t_info = draw(blocks, F.OUT / f"Figure2_{tissue}_with_clade_MSR")
        F.write_tsv(F.OUT / f"Figure2_{tissue}_with_clade_MSR_values.tsv", vals)
        write_legend(F.OUT / f"Figure2_{tissue}_with_clade_MSR_legend.txt", tissue, t_info)


def draw(blocks, out_prefix):
    FW_mm = FIG_W
    FH_mm = sum(b["head"] + len(b["rows"]) * CELL + b["foot"] for b in blocks) + GAP * (len(blocks) - 1)
    FW, FH = FW_mm * MM, FH_mm * MM
    fig = plt.figure(figsize=(FW, FH))
    ax_at = lambda x0, w, y0, h: fig.add_axes([x0 / FW_mm, y0 / FH_mm, w / FW_mm, h / FH_mm])
    top, extra, checks = FH_mm, {}, []

    for b in blocks:
        bh = len(b["rows"]) * CELL
        top -= b["head"]
        y0 = top - bh

        ax_tree = ax_at(*X_TREE, y0, bh)
        ax_auroc = ax_at(*X_AUROC, y0, bh)
        ax_acc = ax_at(*X_ACC, y0, bh)
        F.draw_tree(ax_tree, b["topo"], b["rows"])
        b["auroc"](ax_auroc)
        b["acc"](ax_acc)
        axes, xs = [ax_tree, ax_auroc, ax_acc], [X_TREE[0], X_AUROC[0], X_ACC[0]]

        heats = b.get("heats")
        x = X_HEAT0
        if heats:
            mode = b.get("heat_letters", "first_only")
            for j, (m_dict, cols, title, group) in enumerate(heats):
                hw = len(cols) * CELL
                ax_h = ax_at(x, hw, y0, bh)
                F.draw_heat(ax_h, m_dict, b["rows"], cols, title)
                F.draw_brackets(ax_at(x, hw, y0 - BRACKET_Y, 8.0), cols)
                if group:
                    fig.text(x / FW_mm, (top + 7.6) / FH_mm, group, fontsize=FS, va="bottom")
                axes.append(ax_h)
                if mode == "each" or (mode == "first_only" and j == 0) or (mode == "pair_first" and j % 2 == 0):
                    xs.append(x)
                x += max(hw, 10.0 if len(heats) > 2 else 30.0) + (2.5 if j % 2 == 0 else 6.0)

        for xx, lab in zip(xs, b["letters"]):
            if lab:
                fig.text(xx / FW_mm, (top + 3.2) / FH_mm, lab, fontsize=FS_L, fontweight="bold", va="bottom")
        fig.text(X_TREE[0] / FW_mm, (top + 7.6) / FH_mm, b["name"], fontsize=FS + 0.5, fontweight="bold", va="bottom")

        if b.get("extra"):
            axp = ax_at(X_PT, FW_mm - X_PT, y0, bh)
            extra = b["extra"](axp)
            fig.text(X_PT / FW_mm, (top + 3.2) / FH_mm, b["extra_letter"], fontsize=FS_L, fontweight="bold", va="bottom")

        if b.get("legend_under_b"):
            F.draw_legend(fig, b["legend_under_b"], X_AUROC[0], y0 - 7.5, FW, FH)
        if b.get("legend_right"):
            items, lx = b["legend_right"]
            F.draw_legend(fig, items, lx, y0 + bh - 1.5, FW, FH)

        if b.get("colourbar"):
            cb_w = 36.0
            cb_x = X_HEAT0 + ((x - 2.5) - X_HEAT0 - cb_w) / 2.0
            axcb = ax_at(cb_x, cb_w, y0 - CB_Y, 2.0)
            cb = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, 1), cmap=F.CMAP),
                              cax=axcb, ticks=[0, 0.5, 1], orientation="horizontal")
            cb.ax.tick_params(labelsize=FS_S, length=1.5, width=0.5, pad=1)
            cb.outline.set_linewidth(0.4)
            cb.set_label("Mean predicted probability of high expression", fontsize=FS_S, labelpad=2)

        checks.append((axes, len(b["rows"])))
        top = y0 - b["foot"] - GAP

    for axs, n in checks:
        rr = np.arange(n)
        ys = [a.transData.transform(np.c_[np.zeros_like(rr), rr])[:, 1] for a in axs]
        assert max(float(np.abs(y - ys[0]).max()) for y in ys[1:]) < 0.5, "rows not aligned"
    # after the last block one trailing GAP has been subtracted that is not part of the canvas
    assert top + GAP >= -0.6, f"blocks overflow the canvas by {-(top + GAP):.2f} mm"

    fig.savefig(str(out_prefix) + ".pdf")
    fig.savefig(str(out_prefix) + ".png", dpi=600)
    plt.close(fig)
    print(f"Saved {out_prefix} ({FW_mm:.0f} x {FH_mm:.0f} mm)")
    return extra


def write_legend(path, tissue, t):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(
            f"Figure 2 | Flowering plant deepCRE model performance for single-species (SSR) and multi-species (MSR) training ({tissue}).\n"
            f"a, A broad phylogenetic spectrum of flowering plants and crop species was selected, encompassing rosids (Arabidopsis thaliana (thale cress; Atha), Camelina sativa (false flax; Csat), Brassica oleracea (cabbage; Bole), Cucumis sativus (cucumber; Cucu), Glycine max (soy; Gmax)), asterids (Daucus carota (carrot; Dcar), Solanum lycopersicum (tomato; Slyc)), Caryophyllales (Chenopodium quinoa (quinoa; Cqui), Beta vulgaris (sugar beet; Bvul)), and grasses/cereals (Oryza sativa (rice; Osat), Hordeum vulgare (barley; Hvul), Sorghum bicolor (sorghum; Sbic), Zea mays (maize; Zmay)).\n"
            f"b, c, Area under the receiver operating characteristic curve (auROC, b) and accuracy (c) of single-species reference training (SSR) deepCRE models for {tissue} tissue show excellent predictive performance. Following the deepCRE SSR training strategy, cross-validations were performed on individual chromosomes. Grey dots show the performance of each chromosome fold on its hold-out chromosome validation set; the best performing fold per species (selected by maximum validation auROC) is the model provided on the web toolkit and is shown as a black dot with its value printed at the right. The mean and standard deviation across all folds are displayed as an open circle and a horizontal line, respectively (n = number of chromosome folds per species, see Supplementary Table 1).\n"
            f"d, Cross-species evaluation of the SSR models: mean predicted probability (PP) of high expression for genes labelled as lowly expressed (left) and for genes labelled as highly expressed (right), for every combination of model species (rows) and target species (columns). A mean PP near 0 indicates a gene is predicted to be lowly expressed, a PP near 1 highly expressed. The black diagonal line marks models predicting their own native species (including training data). Dashed rectangles emphasise cross-prediction within phylogenetically closely related groups. In the left heatmap, blue cells indicate that more than 50% of lowly expressed genes were predicted correctly, while orange cells show the opposite.\n"
            f"e, UPGMA tree of Euclidean distances between the 52-value cross-prediction profiles of the 13 species (rows and columns of d). Mantel test against the topological distances of the cladogram in a: r = {t['mantel_r']:.3f}, p = {t['mantel_p']:.4f} (Pearson, 9,999 permutations, 78 species pairs) shows that relatedness correlates with cross-species prediction performance. Leaf order: {', '.join(t['order'])}.\n"
            f"f, g, Multi-species reference training (MSR) performance on the hold-out species using auROC (f) and accuracy (g), for MSR models trained on all 13 species (black dot with value at right), in comparison to the SSR web-tool model of that species (open diamond) and to the cross-species predictions of the SSR models of the other 12 species (light grey dots).\n"
            f"h-l, MSR models trained on the 9 eudicot species only (rows: held-out eudicot): cladogram of the eudicots (h); hold-out species auROC (i) and accuracy (j) as in f and g; cross-species evaluation as in d for the eudicot-only models within their own clade (k) and on the grasses (l).\n"
            f"m-q, MSR models trained on the 4 grass species only (rows: held-out grass), shown as in h-l: cladogram of the grasses (m); hold-out species auROC (n) and accuracy (o); cross-species evaluation within the grasses (p) and on the eudicots (q).\n"
        )


if __name__ == "__main__":
    main((sys.argv[1],) if len(sys.argv) > 1 else ("leaf", "root"))
