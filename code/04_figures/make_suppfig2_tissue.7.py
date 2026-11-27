#!/usr/bin/env python3
"""
Purpose:
    Plot Supplementary Figure 2 (`SuppFig_tissue_v7.{pdf,png,_values.tsv,_legend.txt}`): leaf vs.
    root SSR and MSR cross-tissue evaluations (`a–i` with fold counts `n`), 4-category gene
    label composition (`j`), probability contrasts on swapped genes (`k,m`), and classification
    auROC on swapped genes (`l,n,o`).

Manuscript (Peleke et al., 2026):
    Produces Supplementary Figure 2 offline directly from `results/evaluation/final/` and
    `results/evaluation/tissue_cross/`.

Usage & I/O (run from repo root):
    FIG_OUT_DIR=/path/to/out python3 code/04_figures/make_suppfig2_tissue.7.py

How to change or extend:
    Redirect output via `FIG_OUT_DIR`. Modify `TC` path or panel widths at the top of the module,
    or update `ORDER` in `make_pub_figures_split.7.py` to change the species set.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("figs7", HERE / "make_pub_figures_split.7.py")
F = importlib.util.module_from_spec(spec)
spec.loader.exec_module(F)
spec2 = importlib.util.spec_from_file_location("fig2comp2", HERE / "make_fig2_comprehensive.2.py")
C2 = importlib.util.module_from_spec(spec2)
spec2.loader.exec_module(C2)
C2.TREE_W = 0.30     # the SSR cladograms here are 22 mm wide (Figure 2: 24.5 mm), so the branches are compressed further


def fold_counts(tissue):
    pf = [r for r in F.read_tsv_dicts(F.FINAL / "S2v8_per_fold.tsv")
          if r["status"] in ("converged", "converged_rescue") and r["code"] in ORDER
          and r["auROC_v8"] not in ("", "nan", "NA")]
    n = C2.fold_counts(pf, tissue)
    assert sum(n.values()) in (115, 116), (tissue, sum(n.values()))
    return n

MM, CELL, FIG_W = F.MM, F.CELL, F.FIG_W
ORDER, DARK, MID, LIGHT = F.ORDER, F.DARK, F.MID, F.LIGHT
FS, FS_S, FS_L = F.FS, F.FS_S, F.FS_L
TC = Path("results/evaluation/tissue_cross")
OTHER = {"leaf": "root", "root": "leaf"}

# panels are 14 mm wide with a 14 mm gutter: the "loss" column of one panel ends ~4 mm before
# the axis of the next, which it did not at 15 mm / 10 mm
L_TREE, L_A, L_B = (0.0, 22.0), (25.0, 14.0), (53.0, 14.0)
R_TREE, R_A, R_B = (95.0, 22.0), (120.0, 14.0), (148.0, 14.0)
WIDE = 42.0
L_C, L_D = (25.0, 21.0), (52.0, 21.0)      # bottom block: contrast + performance, left half
R_C, R_D = (120.0, 21.0), (147.0, 21.0)    # and right half
HEAD = 9.0
FOOT_SSR, FOOT_MSR, FOOT_COMP, FOOT_CONTRAST = 15.0, 9.0, 15.0, 25.0   # last block: four-item key
GAP = 3.0

KEY_SSR = [
    (dict(marker="o", ms=3.2, mfc=DARK, mec=DARK), "own species, own tissue (value, loss)"),
    (dict(marker="D", ms=2.6, mfc="white", mec=DARK, mew=0.6), "own species, other tissue"),
    (dict(marker="s", ms=2.8, mfc=MID, mec=MID), "other species, own tissue (median)"),
    (dict(marker="o", ms=1.8, mfc=LIGHT, mec=MID, mew=0.2), "other species, own tissue (n = 12 models)"),
]
KEY_MSR = [
    (dict(marker="o", ms=3.2, mfc=DARK, mec=DARK), "held-out species, own tissue (value, loss)"),
    (dict(marker="D", ms=2.6, mfc="white", mec=DARK, mew=0.6), "held-out species, other tissue"),
]
KEY_COMP = [   # white -> black -> dark grey -> light grey, so neighbouring segments stay distinct
    (dict(facecolor="white", edgecolor=MID, linewidth=0.5), "stable class in both tissues"),
    (dict(facecolor=DARK, edgecolor=DARK, linewidth=0.3), "swap between high and low"),
    (dict(facecolor=MID, edgecolor=MID, linewidth=0.3), "high in one tissue, intermediate in the other"),
    (dict(facecolor=LIGHT, edgecolor=MID, linewidth=0.3), "low in one tissue, intermediate in the other"),
]
# circles always mean a group of genes, squares always mean a model, so the filled/open
# convention never has two meanings in the same block
KEY_CONTRAST = [
    (dict(marker="o", ms=3.0, mfc=DARK, mec=DARK), "genes high in leaf, low in root"),
    (dict(marker="o", ms=3.0, mfc="white", mec=DARK, mew=0.6), "genes low in leaf, high in root"),
    (dict(marker="s", ms=2.8, mfc=DARK, mec=DARK), "leaf model"),
    (dict(marker="s", ms=2.8, mfc="white", mec=DARK, mew=0.6), "root model"),
]


def load():
    ssr = pd.read_csv(TC / "tissue_cross_full.tsv", sep="\t")
    msr = pd.read_csv(TC / "tissue_cross_msr_full.tsv", sep="\t")
    msr = msr[msr.msr_set == "all"]
    sw = pd.read_csv(TC / "tissue_switch_direction.tsv", sep="\t")
    swm = pd.read_csv(TC / "tissue_cross_msr_switch_direction.tsv", sep="\t")
    swm = swm[swm.msr_set == "all"].rename(columns={"held_out_species": "species"})
    comp = pd.read_csv(TC / "tissue_label_composition.tsv", sep="\t").set_index("species").reindex(ORDER)
    perf = pd.read_csv(TC / "tissue_switch_performance.tsv", sep="\t")
    perm = pd.read_csv(TC / "tissue_switch_permutation.tsv", sep="\t")
    g = lambda d, k: d[d.gene_group == k].set_index("species").reindex(ORDER)
    return (ssr, msr, {"HL": g(sw, "higher in leaf"), "LH": g(sw, "higher in root")},
            {"HL": g(swm, "high in leaf only"), "LH": g(swm, "high in root only")}, comp, perf, perm)


def value_loss(ax, i, own, other):
    """Own-tissue value and, beside it, what is lost when the other tissue's classes are used."""
    ax.text(1.02, i, f"{own:.2f}", fontsize=4.9, va="center", ha="left",
            transform=ax.get_yaxis_transform(), clip_on=False)
    ax.text(1.40, i, f"{own - other:+.2f}", fontsize=4.9, va="center", ha="left", color=MID,
            transform=ax.get_yaxis_transform(), clip_on=False)


def header(fig, ax, FW_mm, FH_mm, x_mm, w_mm, y_mm):
    fig.text((x_mm + w_mm * 1.02) / FW_mm, y_mm / FH_mm, "own", fontsize=4.6, color=DARK, va="bottom")
    fig.text((x_mm + w_mm * 1.40) / FW_mm, y_mm / FH_mm, "loss", fontsize=4.6, color=MID, va="bottom")


def ssr_panel(ax, d, tissue, metric, xlabel, vals, panel):
    rng = np.random.default_rng(7 if metric == "auROC" else 17)
    for i, code in enumerate(ORDER):
        s = d[(d.model_species == code) & (d.model_tissue == tissue)]
        own = float(s[(s.target_species == code) & (s.label_tissue == tissue)][metric].iloc[0])
        cross = float(s[(s.target_species == code) & (s.label_tissue == OTHER[tissue])][metric].iloc[0])
        others = s[(s.target_species != code) & (s.label_tissue == tissue)][metric].to_numpy()
        ax.scatter(others, i + rng.uniform(-0.22, 0.22, len(others)),
                   s=1.8, color=LIGHT, edgecolor=MID, linewidths=0.2, zorder=2)
        ax.scatter([np.median(others)], [i], s=8, marker="s", color=MID, zorder=3, linewidths=0)
        ax.scatter([cross], [i], s=8, marker="D", facecolor="white", edgecolor=DARK, linewidths=0.6, zorder=4)
        ax.scatter([own], [i], s=11, color=DARK, zorder=5, linewidths=0)
        value_loss(ax, i, own, cross)
        vals.append(dict(panel=panel, model_type="SSR", model_tissue=tissue, species=code, metric=xlabel,
                         own_tissue=round(own, 4), other_tissue=round(cross, 4),
                         loss_to_other_tissue=round(own - cross, 4),
                         other_species_own_tissue_median=round(float(np.median(others)), 4),
                         n_other_species=len(others)))
    ax.set_ylim(len(ORDER) - 0.5, -0.5)
    F.metric_axis(ax, xlabel=xlabel)


def msr_panel(ax, msr, tissue, metric, xlabel, vals, panel):
    for i, code in enumerate(ORDER):
        m = msr[(msr.held_out_species == code) & (msr.model_tissue == tissue)]
        own = float(m[m.label_tissue == tissue][metric].iloc[0])
        cross = float(m[m.label_tissue == OTHER[tissue]][metric].iloc[0])
        ax.plot([min(own, cross), max(own, cross)], [i, i], color=LIGHT, lw=0.5,
                solid_capstyle="butt", zorder=2)
        ax.scatter([cross], [i], s=8, marker="D", facecolor="white", edgecolor=DARK, linewidths=0.6, zorder=4)
        ax.scatter([own], [i], s=11, color=DARK, zorder=5, linewidths=0)
        value_loss(ax, i, own, cross)
        vals.append(dict(panel=panel, model_type="MSR (all 13 species)", model_tissue=tissue, species=code,
                         metric=xlabel, own_tissue=round(own, 4), other_tissue=round(cross, 4),
                         loss_to_other_tissue=round(own - cross, 4)))
    ax.set_ylim(len(ORDER) - 0.5, -0.5)
    F.metric_axis(ax, xlabel=xlabel)


def composition(ax, comp, vals):
    seg = [("stable", lambda r: r.n_stable, dict(facecolor="white", edgecolor=MID, linewidth=0.5)),
           ("swapped", lambda r: r.n_swapped, dict(color=DARK, linewidth=0)),
           ("high to intermediate", lambda r: r.n_high_to_intermediate, dict(color=MID, linewidth=0)),
           ("low to intermediate", lambda r: r.n_low_to_intermediate,
            dict(facecolor=LIGHT, edgecolor=MID, linewidth=0.3))]
    for i, code in enumerate(ORDER):
        r = comp.loc[code]
        tot = r.n_classed_in_at_least_one
        left = 0.0
        for _, f, kw in seg:
            w = 100.0 * f(r) / tot
            ax.barh(i, w, left=left, height=0.62, zorder=3, **kw)
            left += w
        ax.text(101.5, i, f"{int(tot):,}", fontsize=4.9, va="center", ha="left")
        vals.append(dict(panel="j", species=code, n_classed_in_at_least_one_tissue=int(tot),
                         n_stable=int(r.n_stable), n_swapped=int(r.n_swapped),
                         n_high_to_intermediate=int(r.n_high_to_intermediate),
                         n_low_to_intermediate=int(r.n_low_to_intermediate),
                         pct_stable=round(100.0 * r.n_stable / tot, 2),
                         pct_swapped=round(100.0 * r.n_swapped / tot, 2),
                         pct_high_to_intermediate=round(100.0 * r.n_high_to_intermediate / tot, 2),
                         pct_low_to_intermediate=round(100.0 * r.n_low_to_intermediate / tot, 2)))
    ax.set_ylim(len(ORDER) - 0.5, -0.5)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 50, 100])
    ax.set_xticklabels(["0", "50", "100"], fontsize=FS_S)
    ax.set_yticks([])
    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(axis="x", length=1.5, width=0.5, pad=1)
    ax.set_xlabel("Genes classed in at least one tissue (%)", fontsize=FS, labelpad=1)


def flipped_bars(ax, sw, vals):
    L, R = sw["HL"], sw["LH"]
    lim = 1.30 * max(L.n_genes.max(), R.n_genes.max())
    off = 0.04 * lim
    for i, code in enumerate(ORDER):
        nl, nr = int(L.loc[code, "n_genes"]), int(R.loc[code, "n_genes"])
        ax.barh(i, -nl, height=0.55, color=DARK, linewidth=0, zorder=3)
        ax.barh(i, nr, height=0.55, facecolor="white", edgecolor=DARK, linewidth=0.6, zorder=3)
        ax.text(-nl - off, i, f"{nl:,}", fontsize=4.9, va="center", ha="right")
        ax.text(nr + off, i, f"{nr:,}", fontsize=4.9, va="center", ha="left")
        vals.append(dict(panel="k", species=code, n_high_leaf_low_root=nl, n_low_leaf_high_root=nr))
    ax.axvline(0, color=DARK, lw=0.5, zorder=4)
    ax.set_ylim(len(ORDER) - 0.5, -0.5)
    ax.set_xlim(-lim, lim)
    t = [-600, -300, 0, 300, 600]
    ax.set_xticks(t)
    ax.set_xticklabels([f"{abs(x)}" for x in t], fontsize=FS_S)
    ax.set_yticks([])
    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(axis="x", length=1.5, width=0.5, pad=1)
    ax.set_xlabel("Genes that swap expression class (n)", fontsize=FS, labelpad=1)
    ax.text(0.0, -0.215, "high in leaf, low in root", fontsize=FS_S, ha="left", va="top", color=MID,
            transform=ax.transAxes)
    ax.text(1.0, -0.215, "low in leaf, high in root", fontsize=FS_S, ha="right", va="top", color=MID,
            transform=ax.transAxes)


def contrast(ax, sw, vals, panel, label, perm=None):
    L, R = sw["HL"], sw["LH"]
    for i, code in enumerate(ORDER):
        a = float(L.loc[code, "contrast_leaf_minus_root"])
        b = float(R.loc[code, "contrast_leaf_minus_root"])
        ax.plot([min(a, b), max(a, b)], [i, i], color=MID, lw=0.5, solid_capstyle="butt", zorder=2)
        ax.scatter([a], [i], s=10, color=DARK, zorder=4, linewidths=0)
        ax.scatter([b], [i], s=9, facecolor="white", edgecolor=DARK, linewidths=0.6, zorder=4)
        vals.append(dict(panel=panel, model_type=label, species=code,
                         contrast_high_leaf_low_root=round(a, 4),
                         contrast_low_leaf_high_root=round(b, 4),
                         mean_P_leaf_model_high_leaf=float(L.loc[code, "mean_P_leaf_model"]),
                         mean_P_root_model_high_leaf=float(L.loc[code, "mean_P_root_model"]),
                         mean_P_leaf_model_high_root=float(R.loc[code, "mean_P_leaf_model"]),
                         mean_P_root_model_high_root=float(R.loc[code, "mean_P_root_model"]),
                         permutation_p=(float(perm.loc[code]) if perm is not None else None)))
        if perm is not None and float(perm.loc[code]) >= 0.001:
            ax.text(0.405, i, "n.s.", fontsize=4.4, va="center", ha="right", color=MID)
    ax.axvline(0, color=MID, lw=0.5, ls=(0, (2, 2)), zorder=1)
    ax.set_ylim(len(ORDER) - 0.5, -0.5)
    ax.set_xlim(-0.42, 0.42)
    ax.set_xticks([-0.3, 0.0, 0.3])
    ax.set_xticklabels(["-0.3", "0", "0.3"], fontsize=FS_S)
    ax.set_yticks([])
    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(axis="x", length=1.5, width=0.5, pad=1)
    ax.set_xlabel("P(high): leaf - root model", fontsize=FS, labelpad=1)
    ax.text(0.0, -0.215, "root higher", fontsize=FS_S, ha="left", va="top", color=MID,
            transform=ax.transAxes)
    ax.text(1.0, -0.215, "leaf higher", fontsize=FS_S, ha="right", va="top", color=MID,
            transform=ax.transAxes)


def flipped_performance(ax, perf, kind, vals, panel):
    """How well each model actually classifies the swapped genes, scored against the classes of its
    own tissue. The two label sets are exact complements there, so 0.5 is chance."""
    for i, code in enumerate(ORDER):
        s = perf[(perf.model_type == kind) & (perf.species == code)].set_index("model_tissue")
        lf, rt = float(s.loc["leaf", "auROC"]), float(s.loc["root", "auROC"])
        ax.plot([min(lf, rt), max(lf, rt)], [i, i], color=MID, lw=0.5, solid_capstyle="butt", zorder=2)
        ax.scatter([lf], [i], s=9, marker="s", color=DARK, zorder=4, linewidths=0)
        ax.scatter([rt], [i], s=8, marker="s", facecolor="white", edgecolor=DARK, linewidths=0.6, zorder=4)
        vals.append(dict(panel=panel, model_type=kind, species=code,
                         auROC_leaf_model=round(lf, 3), auROC_root_model=round(rt, 3),
                         accuracy_leaf_model=round(float(s.loc["leaf", "accuracy"]), 3),
                         accuracy_root_model=round(float(s.loc["root", "accuracy"]), 3),
                         MCC_leaf_model=round(float(s.loc["leaf", "MCC"]), 3),
                         MCC_root_model=round(float(s.loc["root", "MCC"]), 3),
                         n_swapped=int(s.loc["leaf", "n_flipped"])))
    ax.axvline(0.5, color=MID, lw=0.5, ls=(0, (2, 2)), zorder=1)
    ax.set_ylim(len(ORDER) - 0.5, -0.5)
    ax.set_xlim(0.35, 0.90)
    ax.set_xticks([0.4, 0.6, 0.8])
    ax.set_xticklabels(["0.4", "0.6", "0.8"], fontsize=FS_S)
    ax.set_yticks([])
    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(axis="x", length=1.5, width=0.5, pad=1)
    ax.set_xlabel("auROC on swapped genes", fontsize=FS, labelpad=1)
    ax.text(0.5, -0.215, "chance", fontsize=FS_S, ha="center", va="top", color=MID,
            transform=ax.transAxes)


def main():
    F.setup_rc()
    ssr, msr, sw, swm, comp, perf, perm = load()
    vals = []
    bh = len(ORDER) * CELL
    FW_mm = FIG_W
    FH_mm = 4 * (HEAD + bh) + FOOT_SSR + FOOT_MSR + FOOT_COMP + FOOT_CONTRAST + 3 * GAP
    FW, FH = FW_mm * MM, FH_mm * MM
    fig = plt.figure(figsize=(FW, FH))
    ax_at = lambda x0, w, y0, h: fig.add_axes([x0 / FW_mm, y0 / FH_mm, w / FW_mm, h / FH_mm])
    letter = lambda x, top, s: fig.text(x / FW_mm, (top + 3.0) / FH_mm, s, fontsize=FS_L,
                                        fontweight="bold", va="bottom")
    title = lambda x, top, s: fig.text(x / FW_mm, (top + 7.0) / FH_mm, s, fontsize=FS + 0.5,
                                       fontweight="bold", va="bottom")

    # --- single-species models -------------------------------------------------------
    top = FH_mm - HEAD
    y0 = top - bh
    for (tr, xa, xb), tis, lets in (((L_TREE, L_A, L_B), "leaf", ("a", "b", "c")),
                                    ((R_TREE, R_A, R_B), "root", (None, "d", "e"))):
        C2.make_tree_with_n(fold_counts(tis))(ax_at(*tr, y0, bh), F.TOPO_ALL, ORDER)
        for (x, w), met, xl, lab in ((xa, "auROC", "auROC", lets[1]), (xb, "accuracy", "Accuracy", lets[2])):
            ax = ax_at(x, w, y0, bh)
            ssr_panel(ax, ssr, tis, met, xl, vals, lab)
            header(fig, ax, FW_mm, FH_mm, x, w, top + 0.3)
            letter(x, top, lab)
        if lets[0]:
            letter(tr[0], top, lets[0])
        title(tr[0], top, f"Single-species models (SSR), {tis}")
    F.draw_legend(fig, KEY_SSR, L_A[0], y0 - 6.5, FW, FH)

    # --- multi-species models --------------------------------------------------------
    top = y0 - FOOT_SSR - GAP - HEAD
    y0 = top - bh
    for (tr, xa, xb), tis, lets in (((L_TREE, L_A, L_B), "leaf", ("f", "g")),
                                    ((R_TREE, R_A, R_B), "root", ("h", "i"))):
        F.draw_tree(ax_at(*tr, y0, bh), F.TOPO_ALL, ORDER)
        for (x, w), met, xl, lab in ((xa, "auROC", "auROC", lets[0]), (xb, "accuracy", "Accuracy", lets[1])):
            ax = ax_at(x, w, y0, bh)
            msr_panel(ax, msr, tis, met, xl, vals, lab)
            header(fig, ax, FW_mm, FH_mm, x, w, top + 0.3)
            letter(x, top, lab)
        title(tr[0], top, f"Multi-species models (MSR), {tis} (row: held-out species)")
    F.draw_legend(fig, KEY_MSR, L_A[0], y0 - 6.5, FW, FH)

    # --- class composition and the flipped genes --------------------------------------
    top = y0 - FOOT_MSR - GAP - HEAD
    y0 = top - bh
    F.draw_tree(ax_at(*L_TREE, y0, bh), F.TOPO_ALL, ORDER)
    F.draw_tree(ax_at(*R_TREE, y0, bh), F.TOPO_ALL, ORDER)
    composition(ax_at(L_A[0], WIDE, y0, bh), comp, vals)
    flipped_bars(ax_at(R_A[0], WIDE, y0, bh), sw, vals)
    letter(L_A[0], top, "j")
    letter(R_A[0], top, "k")
    title(L_TREE[0], top, "Expression class of the same gene in leaf and in root")
    F.draw_legend(fig, KEY_COMP, L_A[0], y0 - 6.5, FW, FH)

    # --- what the models predict for the flipped genes --------------------------------
    top = y0 - FOOT_COMP - GAP - HEAD
    y0 = top - bh
    F.draw_tree(ax_at(*L_TREE, y0, bh), F.TOPO_ALL, ORDER)
    F.draw_tree(ax_at(*R_TREE, y0, bh), F.TOPO_ALL, ORDER)
    ps = perm[perm.model == "SSR"].set_index("species").p_perm.reindex(ORDER)
    pm = perm[perm.model == "MSR"].set_index("species").p_perm.reindex(ORDER)
    contrast(ax_at(*L_C, y0, bh), sw, vals, "l", "SSR", ps)
    flipped_performance(ax_at(*L_D, y0, bh), perf, "SSR", vals, "m")
    contrast(ax_at(*R_C, y0, bh), swm, vals, "n", "MSR (all 13 species)", pm)
    flipped_performance(ax_at(*R_D, y0, bh), perf, "MSR (all 13 species)", vals, "o")
    for x, lab in ((L_C[0], "l"), (L_D[0], "m"), (R_C[0], "n"), (R_D[0], "o")):
        letter(x, top, lab)
    title(L_TREE[0], top, "Genes that swap expression class: single-species (l, m) and "
                          "multi-species (n, o) models")
    F.draw_legend(fig, KEY_CONTRAST, L_C[0], y0 - 12.0, FW, FH)

    out = F.OUT / "Supplementary_Figure_2_tissue_n"
    fig.savefig(str(out) + ".pdf")
    fig.savefig(str(out) + ".png", dpi=600)
    plt.close(fig)
    pd.DataFrame(vals).to_csv(str(out) + "_values.tsv", sep="\t", index=False)
    write_legend(str(out) + "_legend.txt", ssr, msr, sw, swm, comp, perf, perm)
    print(f"Saved {out} ({FW_mm:.0f} x {FH_mm:.0f} mm), {len(vals)} value rows")


def med(perf, kind, tissue, col):
    return perf[(perf.model_type == kind) & (perf.model_tissue == tissue)][col].median()


def write_legend(path, ssr, msr, sw, swm, comp, perf, perm):
    own = ssr[ssr.model_species == ssr.target_species]
    oth = ssr[ssr.model_species != ssr.target_species]
    md = lambda x: x.MCC.median()
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(
            "Supplementary Figure 2 | The species boundary costs more than the tissue boundary, and the "
            "models tell the two tissues apart on the genes that swap expression class.\n"
            "Every deployed single-species (SSR) model and every multi-species (MSR) model trained on all "
            "thirteen species was scored against the expression classes of both tissues. No new predictions "
            "were needed: the models score sequence, so only the class labels change. All panels use the "
            "genes a species has predictions for in both tissues. Metrics use the balanced-draw procedure "
            "of Figure 2 (100 seeded draws of the class-balanced gene set, threshold 0.5). Sample sizes: "
            "b-i, n = 1,814-20,370 class-balanced genes per species and tissue for every point (exact n per "
            "species in the Source Data); the SSR models rest on the number of chromosome folds n printed beside the cladograms of the "
            f"single-species block (leaf: n = {C2.n_sentence(fold_counts('leaf'))}; root: n = {C2.n_sentence(fold_counts('root'))}) and the grey clouds on n = 12 models. j, n printed at the right of "
            "each bar. k-m, n printed in k.\n"
            "a, The thirteen species and their phylogenetic relationships, as in Figure 2a; the cladogram "
            "gives the row order of all panels; n beside the two single-species cladograms is the number of chromosome folds per species (leaf, root).\n"
            "b-e, Single-species models: auROC (b, d) and accuracy (c, e) of each species' deployed leaf "
            "model (b, c) and root model (d, e) on its own species, scored against the classes of its own "
            "tissue (black dot) and against the classes of the other tissue (open diamond). The two numbers "
            "at the right of each row are the own-tissue value and the loss incurred when the other "
            "tissue's classes are used. The grey cloud and square are the same model's own-tissue "
            "performance on the twelve other species and its median, which place the tissue cost on the "
            "scale of the species cost.\n"
            "f-i, The same contrast for the multi-species models trained on all thirteen species, each "
            "scored on the species it held out: own tissue (black dot) against the other tissue of the same "
            "species (open diamond), with value and loss printed.\n"
            f"Median MCC for a single-species model was {md(own[own.comparison=='same tissue']):.3f} on its "
            f"own species and tissue against {md(own[own.comparison=='cross-tissue']):.3f} on the other "
            f"tissue, while a model applied to another species in the matching tissue reached only "
            f"{md(oth[oth.comparison=='same tissue']):.3f}: the species boundary costs about four times "
            "what the tissue boundary costs, so a user whose tissue is unavailable is better served by the "
            "other tissue of their own species than by the matching tissue of another species. The "
            "multi-species models lost "
            f"{md(msr[msr.comparison=='same tissue'])-md(msr[msr.comparison=='cross-tissue']):.3f} MCC on "
            f"their held-out species, against {md(oth[oth.comparison=='same tissue'])-md(oth[oth.comparison=='cross-tissue']):.3f} "
            "for a single-species model on a species it has not seen.\n"
            "j, How the expression class of a gene relates between the two tissues, as a percentage of the "
            "genes classed in at least one tissue (total printed at the right). A gene is classed as highly "
            "expressed, lowly expressed or intermediate in each tissue, so it either keeps a stable class in "
            "both, swaps between the high and the low class, or moves from high or from low to intermediate. "
            f"Pooled over the thirteen species, {comp.pct_stable.median():.1f}% of genes were stable, "
            f"{comp.pct_high_to_intermediate.median():.1f}% moved from high to intermediate, "
            f"{comp.pct_low_to_intermediate.median():.1f}% from low to intermediate and only "
            f"{comp.pct_swapped.median():.1f}% swapped class (per-species medians). That small swapping "
            "fraction is why the open diamonds in b-i sit close to the black dots: a model scored against "
            "the other tissue is mostly being asked the same question again. Genes that move to the "
            "intermediate range carry no class in the other tissue and cannot enter l-o.\n"
            "k, The genes that swap expression class, split by direction: high in leaf and low in root "
            "(filled, left) and low in leaf and high in root (open, right). The two groups are the same "
            "genes seen from either side, so a gene high in leaf only is by construction low in root only.\n"
            "l, n, What the models predict for the genes that swap expression class: the difference in mean "
            "predicted probability of high expression between the leaf and the root model of the same "
            "kind, for the single-species models (l) and the multi-species models (n). Filled circles, "
            "genes high in leaf and low in root; open circles, the opposite. A filled circle right of the "
            "dashed line means the leaf model calls those genes more highly expressed than the root model "
            f"does, which is the expected result. Medians {sw['HL'].contrast_leaf_minus_root.median():+.3f} "
            f"and {sw['LH'].contrast_leaf_minus_root.median():+.3f} for the single-species models and "
            f"{swm['HL'].contrast_leaf_minus_root.median():+.3f} and "
            f"{swm['LH'].contrast_leaf_minus_root.median():+.3f} for the multi-species models. The "
            "separation between the two gene groups was tested per species by permuting which swapping "
            "genes belong to which group (10,000 permutations); it exceeded chance in all thirteen species "
            "for the single-species models and in twelve for the multi-species models, the exception being "
            "marked n.s.\n"
            "m, o, How well the same models classify those genes, scored against the classes of their own "
            "tissue: filled squares, leaf model; open squares, root model. On this gene set the two label "
            "sets are exact complements, so 0.5 (dashed line) is chance. Single-species models reached a "
            f"median auROC of {med(perf,'SSR','leaf','auROC'):.3f} in leaf and {med(perf,'SSR','root','auROC'):.3f} "
            f"in root (accuracy {med(perf,'SSR','leaf','accuracy'):.3f} and {med(perf,'SSR','root','accuracy'):.3f}; "
            f"MCC {med(perf,'SSR','leaf','MCC'):+.3f} and {med(perf,'SSR','root','MCC'):+.3f}), above chance in "
            "all thirteen species but far below the auROC these models reach on the gene set as a whole. "
            f"Multi-species models reached {med(perf,'MSR (all 13 species)','leaf','auROC'):.3f} and "
            f"{med(perf,'MSR (all 13 species)','root','auROC'):.3f}, the root models close to chance. The "
            "direction of the effect is therefore systematic while its resolution for individual genes is "
            "poor, and the models should not be used to call differential expression gene by gene. Note "
            f"that C. quinoa contributes few such genes ({int(sw['HL'].loc['Cqui','n_genes'])} and "
            f"{int(sw['LH'].loc['Cqui','n_genes'])}), because its predictions cover scaffolds rather than "
            "full chromosomes.\n"
            "Source data are provided as Source data files.\n"
        )


if __name__ == "__main__":
    main()
