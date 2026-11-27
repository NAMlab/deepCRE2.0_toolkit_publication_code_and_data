#!/usr/bin/env python3
"""
Purpose:
    Generate the final publication Figure 2 (`Figure2_<tissue>_with_clade_MSR_as`): adds panel
    letter `f` on the all-species MSR cladogram, displays exact chromosome fold counts `n` next
    to panel `a` and in the caption, and standardizes `Dcar` heatmap axis labels.

Manuscript (Peleke et al., 2026):
    Produces Figure 2 (leaf: `a–s` panels, `.pdf`, `.png`, `_values.tsv`, `_legend.txt`). Also
    imported by `make_suppfig1_root.2.py` and `make_suppfig2_tissue.7.py`.

Usage & I/O (run from repo root; offline runnable from `results/evaluation/`):
    FIG_OUT_DIR=/path/to/out python3 code/04_figures/make_fig2_comprehensive.2.py [leaf|root]
    - Defaults to both `leaf` and `root`; writes into `$FIG_OUT_DIR` (or `EXTERNAL/figure_output/...`)

How to change or extend:
    Set `FIG_OUT_DIR` to redirect outputs. Edit `ORDER`, `LABELS`, `TOPO_ALL`, or `CLADES` in
    `make_pub_figures_split.7.py` to add species or adjust phylogeny layout.
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("fig2comp0", HERE / "make_fig2_comprehensive.0.py")
C = importlib.util.module_from_spec(spec)
spec.loader.exec_module(C)
F = C.F

# (3) carrot tick label: v7 printed "Dcau" under the heatmaps while every legend uses "Dcar"
F.COLLAB = {c: c for c in F.ORDER}

ORDER = C.ORDER
HEAD_A, FOOT_A = C.HEAD_A, C.FOOT_A
HEAD, FOOT_HEAT = C.HEAD, C.FOOT_HEAT
L_SSR, L_MSR = C.L_SSR, C.L_MSR

TREE_W = 0.50     # x of the tips in axis units (v7: 1.0); frees room for the n column
N_X = 2.25        # right edge of the n column = right edge of the tree axis (24.5 mm)


def fold_counts(pf, tissue):
    n = {}
    for r in pf:
        if r["tissue"] == tissue:
            n[r["code"]] = n.get(r["code"], 0) + 1
    assert set(n) == set(ORDER), (tissue, sorted(n))
    return n


def make_tree_with_n(nmap):
    """Cladogram as F.draw_tree, with the tips at TREE_W and a right-aligned column 'n' (folds per species)."""
    def draw_tree_n(ax, topo, rows):
        y = {c: i for i, c in enumerate(rows)}
        depth = lambda node: 0 if isinstance(node, str) else 1 + max(depth(c) for c in node)
        D = depth(topo)

        def draw(node):
            if isinstance(node, str):
                return TREE_W, y[node]
            kids = [draw(c) for c in node]
            x = TREE_W - depth(node) * 0.9 * TREE_W / max(D, 1)
            ys = [k[1] for k in kids]
            ax.plot([x, x], [min(ys), max(ys)], color=F.DARK, lw=0.7)
            for kx, ky in kids:
                ax.plot([x, kx], [ky, ky], color=F.DARK, lw=0.7, solid_capstyle="butt")
            return x, (ys[1] if len(kids) == 3 else (min(ys) + max(ys)) / 2)

        x, yy = draw(topo)
        ax.plot([0.0, x], [yy, yy], color=F.DARK, lw=0.7)
        labels, nums = [], []
        for c in rows:
            labels.append(ax.text(TREE_W + 0.05, y[c], F.LABELS[c], va="center", style="italic", fontsize=F.FS))
            nums.append(ax.text(N_X, y[c], str(nmap[c]), va="center", ha="right", fontsize=F.FS_S))
        ax.text(N_X, -0.95, "n", va="center", ha="right", style="italic", fontsize=F.FS_S)
        ax.set_ylim(len(rows) - 0.5, -0.5)
        ax.set_xlim(0, 2.25)
        ax.axis("off")
        # species names must not run into the n column
        r = ax.figure.canvas.get_renderer()
        gap = min(nb.get_window_extent(r).x0 - lb.get_window_extent(r).x1 for lb, nb in zip(labels, nums))
        assert gap > 8, f"species label touches n column (gap {gap:.1f} px)"
    return draw_tree_n


def main(tissues=("leaf", "root")):
    F.setup_rc()
    pf_all = F.read_tsv_dicts(F.FINAL / "S2v8_per_fold.tsv")
    pf = [r for r in pf_all
          if r["status"] in ("converged", "converged_rescue") and r["code"] in ORDER
          and r["auROC_v8"] not in ("", "nan", "NA")]
    D, G = ORDER[:9], ORDER[9:]
    out = {}

    for tissue in tissues:
        full_rows = F.read_tsv_dicts(F.MSR_DIR / f"msr_eval_{tissue}_full.tsv")
        summary, oth_rows = F.compute_msr_summary(full_rows, pf_all, tissue)
        e0, _ = F.read_matrix_tsv(F.FINAL / f"crossspecies_mean_{tissue}_class0.tsv")
        e1, _ = F.read_matrix_tsv(F.FINAL / f"crossspecies_mean_{tissue}_class1.tsv")
        M = {}
        for s in ("all", "dicots", "monocots"):
            for k in (0, 1):
                M[(s, k)], _ = F.read_matrix_tsv(F.MSR_DIR / f"suppfig_msr_mean_{tissue}_{s}_class{k}.tsv")
        nmap = fold_counts(pf, tissue)

        vals = []
        blocks = [
            dict(name="Single-species models (SSR)",
                 topo=F.TOPO_ALL, rows=ORDER, head=HEAD_A, foot=FOOT_A, tree=make_tree_with_n(nmap),
                 letters=["a", "b", "c", "d"], extra_letter="e", heat_letters="first_only",
                 heats=[(e0, ORDER, "True lowly expressed genes", "Cross-species evaluation (row: model, column: target)"),
                        (e1, ORDER, "True highly expressed genes", None)],
                 auroc=lambda ax: F.ssr_dots(ax, pf, tissue, ORDER, "auROC_v8", "auROC", vals, "b"),
                 acc=lambda ax: F.ssr_dots(ax, pf, tissue, ORDER, "accuracy_v8", "Accuracy", vals, "c"),
                 extra=lambda ax: F.draw_perf_tree(ax, e0, e1),
                 legend_under_b=L_SSR, colourbar=True),
            dict(name="Multi-species models (MSR) — all 13 species",
                 topo=F.TOPO_ALL, rows=ORDER, head=HEAD, foot=FOOT_HEAT,
                 letters=["f", "g", "h", "i"], heat_letters="first_only",
                 heats=[(M[("all", 0)], ORDER, "True lowly expressed genes",
                         "Cross-species evaluation (row: held-out species, column: target)"),
                        (M[("all", 1)], ORDER, "True highly expressed genes", None)],
                 auroc=lambda ax: F.msr_dots(ax, summary, oth_rows, ORDER, "all", None, "auROC", "auROC", vals, "g"),
                 acc=lambda ax: F.msr_dots(ax, summary, oth_rows, ORDER, "all", None, "accuracy", "Accuracy", vals, "h"),
                 legend_under_b=L_MSR),
            dict(name="Multi-species models (MSR) — eudicots only (row: held-out eudicot)",
                 topo=F.TOPO_DIC, rows=D, head=HEAD, foot=FOOT_HEAT,
                 letters=["j", "k", "l", "m", "n"], heat_letters="pair_first",
                 heats=[(M[("dicots", 0)], D, "True low", "on eudicots (own clade)"),
                        (M[("dicots", 1)], D, "True high", None),
                        (M[("dicots", 0)], G, "True low", "on grasses (other clade; cross-species)"),
                        (M[("dicots", 1)], G, "True high", None)],
                 auroc=lambda ax: F.msr_dots(ax, summary, oth_rows, D, "dicots", None, "auROC", "auROC", vals, "k"),
                 acc=lambda ax: F.msr_dots(ax, summary, oth_rows, D, "dicots", None, "accuracy", "Accuracy", vals, "l")),
            dict(name="Multi-species models (MSR) — grasses only (row: held-out grass)",
                 topo=F.TOPO_GRA, rows=G, head=HEAD, foot=FOOT_HEAT,
                 letters=["o", "p", "q", "r", "s"], heat_letters="pair_first",
                 heats=[(M[("monocots", 0)], G, "True low", "on grasses (own clade)"),
                        (M[("monocots", 1)], G, "True high", None),
                        (M[("monocots", 0)], D, "True low", "on eudicots (other clade; cross-species)"),
                        (M[("monocots", 1)], D, "True high", None)],
                 auroc=lambda ax: F.msr_dots(ax, summary, oth_rows, G, "monocots", None, "auROC", "auROC", vals, "p"),
                 acc=lambda ax: F.msr_dots(ax, summary, oth_rows, G, "monocots", None, "accuracy", "Accuracy", vals, "q")),
        ]
        prefix = F.OUT / f"Figure2_{tissue}_with_clade_MSR_as"
        t_info = draw(blocks, prefix)
        F.write_tsv(Path(str(prefix) + "_values.tsv"), vals)
        write_legend(Path(str(prefix) + "_legend.txt"), tissue, t_info, nmap)
        out[tissue] = (t_info, nmap)
    return out


def draw(blocks, out_prefix):
    """v0's layout engine; a block may bring its own cladogram function ('tree')."""
    orig = F.draw_tree
    trees = iter([b.get("tree", orig) for b in blocks])
    F.draw_tree = lambda ax, topo, rows: next(trees)(ax, topo, rows)
    try:
        return C.draw(blocks, out_prefix)
    finally:
        F.draw_tree = orig


def n_sentence(nmap):
    """'n = 5 (A. thaliana), 7 (Cu. sativus, H. vulgare), ...' grouped by fold count, species in figure order."""
    by = {}
    for c in ORDER:
        by.setdefault(nmap[c], []).append(F.LABELS[c])
    return "; ".join(f"{k} ({', '.join(v)})" for k, v in sorted(by.items()))


def write_legend(path, tissue, t, nmap):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(
            f"Figure 2 | Flowering plant deepCRE model performance for single-species (SSR) and multi-species (MSR) training ({tissue}).\n"
            f"a, A broad phylogenetic spectrum of flowering plants and crop species was selected, encompassing rosids (Arabidopsis thaliana (thale cress; Atha), Camelina sativa (false flax; Csat), Brassica oleracea (cabbage; Bole), Cucumis sativus (cucumber; Cucu), Glycine max (soy; Gmax)), asterids (Daucus carota (carrot; Dcar), Solanum lycopersicum (tomato; Slyc)), Caryophyllales (Chenopodium quinoa (quinoa; Cqui), Beta vulgaris (sugar beet; Bvul)), and grasses/cereals (Oryza sativa (rice; Osat), Hordeum vulgare (barley; Hvul), Sorghum bicolor (sorghum; Sbic), Zea mays (maize; Zmay)). The cladogram gives the species order of all panels; n is the number of chromosome folds (SSR models) per species shown in b and c.\n"
            f"b, c, Area under the receiver operating characteristic curve (auROC, b) and accuracy (c) of single-species reference training (SSR) deepCRE models for {tissue} tissue. Following the deepCRE SSR training strategy, cross-validations were performed on individual chromosomes. Grey dots show the performance of each chromosome fold on its hold-out chromosome validation set; the best performing fold per species (selected by maximum validation auROC) is the model provided on the web toolkit and is shown as a black dot with its value printed at the right. The mean and standard deviation across all folds are displayed as an open circle and a horizontal line, respectively. Number of chromosome folds per species: n = {n_sentence(nmap)}.\n"
            f"d, Cross-species evaluation of the SSR models: mean predicted probability (PP) of high expression for genes labelled as lowly expressed (left) and for genes labelled as highly expressed (right), for every combination of model species (rows) and target species (columns). A mean PP near 0 indicates a gene is predicted to be lowly expressed, a PP near 1 highly expressed. The black diagonal line marks models predicting their own native species (including training data). Dashed rectangles emphasise cross-prediction within phylogenetically closely related groups. In the left heatmap, blue cells indicate that more than 50% of lowly expressed genes were predicted correctly, while orange cells show the opposite.\n"
            f"e, UPGMA tree of Euclidean distances between the 52-value cross-prediction profiles of the 13 species (rows and columns of d). Mantel test against the topological distances of the cladogram in a: r = {t['mantel_r']:.3f}, p = {t['mantel_p']:.4f} (Pearson, 9,999 permutations, 78 species pairs) shows that relatedness correlates with cross-species prediction performance. Leaf order: {', '.join(t['order'])}.\n"
            f"f-i, MSR models trained on all 13 species (rows: held-out species): cladogram (f); performance on the hold-out species using auROC (g) and accuracy (h) (black dot with value at right), in comparison to the SSR web-tool model of that species (open diamond) and to the cross-species predictions of the SSR models of the other 12 species (light grey dots). The MSR model is scored on all low- and high-expression genes of the held-out species, the SSR web-tool model on the paralog-filtered genes of its own held-out chromosome. (i) Cross-species evaluation of the all-species MSR models, shown as in d: mean predicted probability of high expression for true lowly (left) and truly highly expressed genes (right), for every combination of held-out species of the model (rows) and target species (columns). The diagonal is the held-out species of each model, i.e. the only column of that row that the model has not seen in training.\n"
            f"j-n, MSR models trained on the 9 eudicot species only (rows: held-out eudicot): cladogram of the eudicots (j); hold-out species auROC (k) and accuracy (l) as in g and h; cross-species evaluation as in d for the eudicot-only models within their own clade (m) and on the grasses (n).\n"
            f"o-s, MSR models trained on the 4 grass species only (rows: held-out grass), shown as in j-n: cladogram of the grasses (o); hold-out species auROC (p) and accuracy (q); cross-species evaluation within the grasses (r) and on the eudicots (s).\n"
            f"Source data are provided as Source data files.\n"
        )


if __name__ == "__main__":
    main((sys.argv[1],) if len(sys.argv) > 1 else ("leaf", "root"))
