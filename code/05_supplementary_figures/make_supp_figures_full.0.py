"""
Purpose:
    Earlier rebuild script for Supplementary Figures 5 (including *ZmUBI1* swaps) and 6 from
    all-model intermediate tables, plus full-panel bar charts across all deployed models.

Manuscript (Peleke et al., 2026):
    Predecessor to `make_supp_figures_full.15.py` (final published Figure 5 layout uses `.15.py`).

Usage & I/O (Python 3.11 + `requirements_analysis.txt`):
    python3 code/05_supplementary_figures/make_supp_figures_full.0.py \
        <results_dir> <scratch_dir> <out_dir>

How to change or extend:
    Use `make_supp_figures_full.15.py` for the final Supplementary Figure 5 (`5a–o`) layout.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, FancyArrow
from scipy.stats import ttest_1samp, ttest_rel, ttest_ind, mannwhitneyu

ORANGE, NAVY, MAUVE, PURPLE = "#E87028", "#003D64", "#BB7FB3", "#581878"
PEACH, ROSE, LORANGE, LBLUE = "#FED5B2", "#F299A0", "#FFA269", "#B8D8E8"
GREY = "#8A8A8A"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "savefig.dpi": 350, "figure.dpi": 110,
})

FOUR = ["Arabidopsis thaliana", "Solanum lycopersicum", "Sorghum bicolor", "Zea mays"]
SHORT = ["(1) A. thaliana", "(2) S. lycopersicum", "(3) S. bicolor", "(4) Z. mays"]
REGIONS = ["gUR", "gTUR", "gUR+gTUR"]
REG_COL = {"gUR": NAVY, "gTUR": ORANGE, "gUR+gTUR": MAUVE}
MEF1 = "AT5G52630"
BASE_COL = {"A": "#2E8B57", "C": NAVY, "G": "#222222", "T": "#C0392B", "N": GREY}
CONS_LO, CONS_HI = 1185, 1230


def star(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s."


def bars(ax, groups, series, vals, colors, ylab, rot=20):
    """vals[(group, series)] -> per-gene array. Bars mean, errorbars SD, points overlaid."""
    ng, ns = len(groups), len(series)
    w = 0.8 / ns
    rows = []
    rng = np.random.default_rng(7)
    for j, s in enumerate(series):
        for i, g in enumerate(groups):
            v = np.asarray(vals[(g, s)], float)
            v = v[~np.isnan(v)]
            if v.size == 0:
                continue
            xx = i - 0.4 + w * (j + 0.5)
            ax.bar(xx, v.mean(), w * 0.86, yerr=v.std(ddof=1) if v.size > 1 else 0,
                   color=colors[j], edgecolor="black", linewidth=0.5,
                   error_kw=dict(lw=0.7, capsize=2))
            ax.scatter(xx + rng.uniform(-w * 0.22, w * 0.22, v.size), v, s=3.5,
                       color="black", alpha=0.6, linewidths=0, zorder=3)
            t, p = ttest_1samp(v, 0.0) if v.size > 1 else (np.nan, np.nan)
            rows.append((g, s, v.size, float(v.mean()),
                         float(v.std(ddof=1)) if v.size > 1 else np.nan, float(t), float(p)))
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(range(ng))
    ax.set_xticklabels(groups, rotation=rot, ha="right" if rot else "center")
    ax.set_ylabel(ylab)
    ax.set_xlim(-0.6, ng - 0.4)
    return rows


# ── panel 5a: the extraction window ────────────────────────────────────────
def panel_window(ax):
    segs = [("gUR", 1, 1000, NAVY), ("gTUR", 1001, 1500, ORANGE),
            ("N pad", 1501, 1520, GREY), ("gTDR", 1521, 2020, MAUVE),
            ("gDR", 2021, 3020, PURPLE)]
    for name, a, b, c in segs:
        ax.add_patch(Rectangle((a, 0.55), b - a + 1, 0.3, facecolor=c,
                               edgecolor="black", linewidth=0.6))
        ax.text((a + b) / 2, 0.70, name, ha="center", va="center", color="white",
                fontsize=7, fontweight="bold")
        ax.text((a + b) / 2, 0.88, f"{a}-{b}", ha="center", va="bottom", fontsize=6)
    ax.annotate("", xy=(1500, 1.02), xytext=(1521, 1.02),
                arrowprops=dict(arrowstyle="-|>", lw=1.0, color="black"))
    ax.text(1510, 1.06, "TSS / transcript", ha="center", fontsize=6.5)
    for name, a, b, c in segs[:2]:
        ax.add_patch(Rectangle((a, 0.08), b - a + 1, 0.22, facecolor=c, alpha=0.35,
                               edgecolor="black", linewidth=0.6, linestyle="--"))
        ax.annotate("", xy=((a + b) / 2, 0.52), xytext=((a + b) / 2, 0.32),
                    arrowprops=dict(arrowstyle="-|>", lw=1.0, color=c))
    ax.text(750, 0.19, "donor promoter region substituted in place "
                       "(OsACT1, AtACT1, PpACT1, CaMV 35S, ZmUBI1)",
            ha="center", va="center", fontsize=7)
    ax.set_xlim(-40, 3060)
    ax.set_ylim(0, 1.2)
    ax.set_yticks([])
    ax.set_xticks([1, 1000, 1500, 2020, 3020])
    ax.set_xlabel("position in the 3020 bp deepCRE input window (bp)")
    for s in ax.spines.values():
        s.set_visible(False)


# ── Supplementary Figure 5 ─────────────────────────────────────────────────
def suppfig5(sw, extra, sal, meta, donors, gen, out_dir, tag):
    pre = "2024 " if gen == "2024" else "r3 "
    d = sw[sw.generation == gen]
    e = extra[extra.generation == gen]
    stats = {}

    def pull(donor, region, species, tissue, src=None):
        s = (src if src is not None else d)
        s = s[(s.variant == f"{donor} {region}") & (s.model == f"{pre}{species} ({tissue})")]
        return s.set_index("target").delta.sort_index().to_numpy() * 100

    fig = plt.figure(figsize=(13.5, 17.5))
    gs = fig.add_gridspec(6, 4, height_ratios=[0.8, 1.25, 1.25, 1.1, 1.1, 1.25],
                          hspace=0.75, wspace=0.42)

    ax = fig.add_subplot(gs[0, :])
    panel_window(ax)
    ax.set_title("a)", loc="left", fontweight="bold")

    # b, c: OsACT1 leaf/root  |  d: AtACT1 - OsACT1 leaf
    for col, tis, lab in ((0, "leaf", "b"), (1, "root", "c")):
        ax = fig.add_subplot(gs[1, col])
        v = {(sh, r): pull("Os_ACT1", r, sp, tis) for sp, sh in zip(FOUR, SHORT) for r in REGIONS}
        stats[lab] = bars(ax, SHORT, REGIONS, v, [REG_COL[r] for r in REGIONS],
                          f"OsACT1 swap, change in P (%), {tis}")
        ax.set_title(f"{lab})", loc="left", fontweight="bold")

    ax = fig.add_subplot(gs[1, 2])
    v = {(sh, r): pull("At_ACT1", r, sp, "leaf") - pull("Os_ACT1", r, sp, "leaf")
         for sp, sh in zip(FOUR, SHORT) for r in REGIONS}
    stats["d"] = bars(ax, SHORT, REGIONS, v, [REG_COL[r] for r in REGIONS],
                      "AtACT1 minus OsACT1, change in P (%), leaf")
    ax.set_title("d)", loc="left", fontweight="bold")

    # o: ZmUBI1, leaf and root, three regions
    ax = fig.add_subplot(gs[1, 3])
    v = {(sh, r): pull("ZmUBI1", r, sp, "leaf", src=e) for sp, sh in zip(FOUR, SHORT)
         for r in REGIONS}
    stats["o"] = bars(ax, SHORT, REGIONS, v, [REG_COL[r] for r in REGIONS],
                      "ZmUBI1 swap, change in P (%), leaf")
    ax.set_title("o)  ZmUBI1 (new donor)", loc="left", fontweight="bold")

    # e: PpACT1 gTUR  |  f: CaMV 35S gUR  |  o2: ZmUBI1 root
    ax = fig.add_subplot(gs[2, 0])
    v = {(sh, t): pull("Pp_ACT1", "gTUR", sp, t) for sp, sh in zip(FOUR, SHORT)
         for t in ("leaf", "root")}
    stats["e"] = bars(ax, SHORT, ["leaf", "root"], v, [NAVY, LBLUE],
                      "PpACT1 gTUR swap, change in P (%)")
    ax.set_title("e)", loc="left", fontweight="bold")

    ax = fig.add_subplot(gs[2, 1])
    v = {(sh, t): pull("CaMV35S", "gUR", sp, t) for sp, sh in zip(FOUR, SHORT)
         for t in ("leaf", "root")}
    stats["f"] = bars(ax, SHORT, ["leaf", "root"], v, [ORANGE, PEACH],
                      "CaMV 35S gUR swap, change in P (%)")
    ax.set_title("f)", loc="left", fontweight="bold")

    ax = fig.add_subplot(gs[2, 2])
    v = {(sh, r): pull("ZmUBI1", r, sp, "root", src=e) for sp, sh in zip(FOUR, SHORT)
         for r in REGIONS}
    stats["o2"] = bars(ax, SHORT, REGIONS, v, [REG_COL[r] for r in REGIONS],
                       "ZmUBI1 swap, change in P (%), root")
    ax.set_title("o')  ZmUBI1, root", loc="left", fontweight="bold")

    ax = fig.add_subplot(gs[2, 3])
    ax.axis("off")
    ax.legend(handles=[Line2D([], [], marker="s", ls="", ms=7, mfc=REG_COL[r], mec="black",
                              label=f"donor region {r}") for r in REGIONS]
              + [Line2D([], [], marker="s", ls="", ms=7, mfc=c, mec="black", label=l)
                 for c, l in ((NAVY, "leaf model (e)"), (LBLUE, "root model (e)"),
                              (ORANGE, "leaf model (f)"), (PEACH, "root model (f)"))]
              + [Line2D([], [], marker="o", ls="", ms=3, color="black",
                        label="individual gene (n = 9)")],
              loc="center", frameon=False)

    # g-l: MEF1 probabilities and saliency for the three published donors
    sal_models = [m for m in (f"{pre}Arabidopsis thaliana (leaf)", f"{pre}Zea mays (leaf)")]
    donors_gl = [("Os_ACT1", "OsACT1", "g", "h"), ("At_ACT1", "AtACT1", "i", "j"),
                 ("Pp_ACT1", "PpACT1", "k", "l"), ("ZmUBI1", "ZmUBI1", "p", "q")]
    for row, col0, (dk, dn, lp, ls_) in ((3, 0, donors_gl[0]), (3, 2, donors_gl[1]),
                                         (4, 0, donors_gl[2]), (4, 2, donors_gl[3])):
        axp = fig.add_subplot(gs[row, col0])
        src = e if dk == "ZmUBI1" else d
        nat, swp, labs = [], [], []
        for m in sal_models:
            base = sw[(sw.model == m) & (sw.target == MEF1)
                      & (sw.variant == "Os_ACT1 gUR")]
            p_nat = float((base.P - base.delta).iloc[0]) if len(base) else np.nan
            s = src[(src.model == m) & (src.target == MEF1) & (src.variant == f"{dk} gTUR")]
            nat.append(p_nat)
            swp.append(float(s.P.iloc[0]) if len(s) else np.nan)
            labs.append(m.replace(pre, "").replace(" (leaf)", "\nleaf"))
        x = np.arange(len(labs))
        axp.bar(x - 0.19, nat, 0.36, color=GREY, edgecolor="black", lw=0.5, label="native")
        axp.bar(x + 0.19, swp, 0.36, color=ORANGE, edgecolor="black", lw=0.5,
                label=f"{dn} gTUR swap")
        axp.axhline(0.5, color="black", lw=0.6, ls="--")
        axp.set_xticks(x)
        axp.set_xticklabels(labs)
        axp.set_ylabel(f"MEF1 ({MEF1})\npredicted probability")
        axp.set_ylim(0, 1)
        axp.legend(frameon=False, loc="upper left")
        axp.set_title(f"{lp})", loc="left", fontweight="bold")
        stats[lp] = [(labs[i], dn, 1, swp[i] - nat[i], np.nan, np.nan, np.nan)
                     for i in range(len(labs))]

        axs = fig.add_subplot(gs[row, col0 + 1])
        m = sal_models[1] if dk in ("Os_ACT1", "ZmUBI1") else sal_models[0]
        kn, ks = f"{m}||MEF1 native", f"{m}||MEF1 {dk} gTUR"
        if kn in sal and ks in sal:
            a = sal[kn].sum(axis=1)
            b = sal[ks].sum(axis=1)
            x = np.arange(950, 1560)
            axs.plot(x, a[950:1560], color=GREY, lw=0.7, label="native")
            axs.plot(x, b[950:1560], color=ORANGE, lw=0.7, label=f"{dn} gTUR")
            axs.axvspan(CONS_LO, CONS_HI, color=MAUVE, alpha=0.18, lw=0)
            axs.axvline(1000.5, color="black", lw=0.6, ls=":")
            axs.set_xlabel("window position (bp)")
            axs.set_ylabel("contribution score")
            axs.legend(frameon=False, fontsize=6, loc="upper left")
            axs.set_title(f"{ls_})  {m}", loc="left", fontweight="bold", fontsize=7.5)
        else:
            axs.axis("off")
            axs.text(0.5, 0.5, "saliency not computed", ha="center", fontsize=7)

    # m: the conserved block
    axm = fig.add_subplot(gs[5, :2])
    cons_pos = [int(v) for v in str(meta.set_index("key").loc["conserved_positions",
                                                              "value"]).split(",") if v]
    seqs3 = {"OsACT1": donors["Os_ACT1"], "AtACT1": donors["At_ACT1"],
             "PpACT1": donors["Pp_ACT1"], "ZmUBI1": donors["ZmUBI1"]}
    names = list(seqs3)
    for r, nm in enumerate(names):
        for j, pos in enumerate(range(CONS_LO, CONS_HI + 1)):
            b = seqs3[nm][pos - 1001]
            cons = pos in cons_pos
            axm.text(j, -r, b, ha="center", va="center", fontsize=6.5,
                     color=BASE_COL.get(b, GREY),
                     fontweight="bold" if cons else "normal",
                     alpha=1.0 if cons else 0.45)
    m = f"{pre}Zea mays (leaf)"
    kn = f"{m}||MEF1 Os_ACT1 gTUR"
    if kn in sal:
        sc = sal[kn].sum(axis=1)[CONS_LO - 1:CONS_HI]
        sc = sc / (np.abs(sc).max() or 1)
        axm.bar(np.arange(len(sc)), sc * 0.9, 0.8, bottom=-len(names) - 0.4,
                color=[ORANGE if v > 0 else NAVY for v in sc], linewidth=0)
    axm.set_yticks(list(range(0, -len(names), -1)) + [-len(names) - 0.4])
    axm.set_yticklabels(names + ["OsACT1 swap\ncontribution"])
    axm.set_xticks(np.arange(0, CONS_HI - CONS_LO + 1, 5))
    axm.set_xticklabels(range(CONS_LO, CONS_HI + 1, 5))
    axm.set_xlabel("window position (bp)")
    axm.set_xlim(-1, CONS_HI - CONS_LO + 1)
    axm.set_ylim(-len(names) - 1.6, 0.8)
    for s in axm.spines.values():
        s.set_visible(False)
    axm.set_title(f"m)  conserved block {CONS_LO}-{CONS_HI}; "
                  f"{len(cons_pos)} columns identical in OsACT1/AtACT1/PpACT1 (bold)",
                  loc="left", fontweight="bold", fontsize=7.5)

    # n: conserved-site mutation, Zm(leaf)
    axn = fig.add_subplot(gs[5, 2:])
    p5n = extra[(extra.analysis.isin(["panel5n", "panel5n_block"]))
                & (extra.generation == gen)]
    m = f"{pre}Zea mays (leaf)"

    def grab(var):
        s = p5n[(p5n.model == m) & (p5n.variant == var)]
        # "OsACT1 gTUR" is present in both the panel5n and the panel5n_block table
        return s.drop_duplicates("target").set_index("target").delta.sort_index()

    series = [("OsACT1 gTUR swap", grab("OsACT1 gTUR"), ORANGE),
              (f"{len(cons_pos)} identical\ncolumns mutated",
               grab("OsACT1 gTUR, conserved sites mutated"), PEACH),
              (f"whole block\n{CONS_LO}-{CONS_HI} mutated",
               grab(f"OsACT1 gTUR, block {CONS_LO}-{CONS_HI} mutated"), MAUVE)]
    series = [(lab, s, c) for lab, s, c in series if len(s)]
    common = series[0][1].index
    for _, s, _c in series[1:]:
        common = common.intersection(s.index)
    raw = [s.loc[common].to_numpy() * 100 for _, s, _c in series]
    # the swap effect is 40-70 % while the mutations move it by about one point,
    # so the panel shows the paired difference from the unmutated swap
    arrs = [v - raw[0] for v in raw[1:]]
    series = series[1:]
    bp = axn.boxplot(arrs, widths=0.55, patch_artist=True,
                     medianprops=dict(color="black", lw=1.0),
                     flierprops=dict(marker="", ls="none"))
    for patch, (_, _s, c) in zip(bp["boxes"], series):
        patch.set_facecolor(c)
        patch.set_edgecolor("black")
        patch.set_linewidth(0.6)
    rng = np.random.default_rng(3)
    for i, v in enumerate(arrs, start=1):
        axn.scatter(i + rng.uniform(-0.1, 0.1, v.size), v, s=9, color="black",
                    alpha=0.7, linewidths=0, zorder=3)
    axn.axhline(0, color="black", lw=0.6)
    axn.set_xticks(range(1, len(arrs) + 1))
    axn.set_xticklabels([lab for lab, _s, _c in series])
    axn.set_ylabel("change in the swap effect\nafter mutation (percentage points)")
    stats["n"] = [("OsACT1 gTUR swap, unmutated", m, raw[0].size, float(raw[0].mean()),
                   float(raw[0].std(ddof=1)), np.nan, np.nan)]
    for (lab, _s, _c), v, r in zip(series, arrs, raw[1:]):
        t, p = ttest_rel(r, raw[0])
        stats["n"].append((lab.replace("\n", " "), m, v.size, float(r.mean()),
                           float(r.std(ddof=1)), float(t), float(p)))
        stats["n"].append((f"  -> paired difference", m, v.size, float(v.mean()),
                           float(v.std(ddof=1)), np.nan, np.nan))
    axn.set_title("n)", loc="left", fontweight="bold")

    png = os.path.join(out_dir, f"SuppFig5_full_{tag}.png")
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)

    with open(png.replace(".png", "_legend.txt"), "w") as fh:
        fh.write(f"Supplementary Figure 5, rebuilt in full, {tag} models "
                 f"(20261109 rerun)\n\n")
        fh.write("a) Segmentation of the 3020 bp deepCRE input window. The gUR and gTUR "
                 "segments are the ones replaced by donor promoter sequence.\n")
        fh.write("b-f, o) Nine Arabidopsis thaliana genes per bar (n = 9); bars are the mean "
                 "change in predicted probability of high expression, error bars the standard "
                 "deviation, points individual genes. Two-sided one-sample t-test of the nine "
                 "per-gene changes against zero.\n")
        fh.write("o, o') ZmUBI1 is a new donor, not in the published figure: the maize "
                 "polyubiquitin-1 locus of Zea mays GCF_902167145.1, NC_050100.1 minus strand, "
                 "1000 bp upstream of the transcription start site as gUR and the first 500 bp "
                 "of the transcript as gTUR. The clone S94464 (Christensen et al. 1992) matches "
                 "this locus at 99.9 % over 1994 bp, differing by one C in a homopolymer at "
                 "clone position 1013.\n")
        fh.write("g, i, k, p) Predicted probability of MEF1 (AT5G52630) natively and after the "
                 "donor gTUR swap. h, j, l, q) DeepSHAP contribution scores summed over the four "
                 "bases, window positions 950-1560; the dotted line is the gUR/gTUR boundary and "
                 "the shaded band the conserved block of panel m.\n")
        fh.write(f"m) Window positions {CONS_LO}-{CONS_HI} of the four donors; bold columns are "
                 f"identical in OsACT1, AtACT1 and PpACT1. The bar track is the contribution "
                 f"score of the OsACT1 gTUR swap under the Zea mays leaf model, scaled to its "
                 f"own maximum.\n")
        fh.write(f"n) Change in predicted probability of the OsACT1 gTUR swap, of the same swap "
                 f"with only the columns identical across the three donors transverted, and of "
                 f"the same swap with every base of window positions {CONS_LO}-{CONS_HI} "
                 f"transverted; same nine genes, lines pair the genes. Two-sided paired t-test "
                 f"of each mutated construct against the unmutated swap. Note that strict "
                 f"three-way column identity in this block is only {len(cons_pos)} of "
                 f"{CONS_HI - CONS_LO + 1} positions, which is chance level for unaligned "
                 f"sequence, so the whole-block mutation is the informative comparison.\n\n")
        for lab in ("b", "c", "d", "o", "e", "f", "o2", "g", "i", "k", "p", "n"):
            if lab not in stats:
                continue
            fh.write(f"panel {lab})\n")
            for g_, s_, n_, m_, sd_, t_, p_ in stats[lab]:
                line = f"  {g_:28s} {s_:12s} n={n_:2d}  mean {m_:+8.3f}"
                if not np.isnan(sd_):
                    line += f"  SD {sd_:6.3f}"
                if not np.isnan(t_):
                    line += f"  t = {t_:+7.3f}  p = {p_:.3g}  {star(p_)}"
                fh.write(line + "\n")
            fh.write("\n")
    return png


# ── swaps as bar charts, whole model panel ─────────────────────────────────
def swaps_allmodels_bars(sw, extra, gen, out_dir, tag):
    d = pd.concat([sw[sw.generation == gen],
                   extra[(extra.analysis == "swap") & (extra.generation == gen)]])
    d = d[d.variant != "reference"].copy()
    d["donor"] = d.variant.str.split(" ").str[0]
    d["region"] = d.variant.str.split(" ").str[1]
    d["tissue"] = d.model.str.extract(r"\((leaf|root)\)")
    d["sp"] = d.model.str.replace(r"^(2024|r3) ", "", regex=True).str.replace(
        r" \((leaf|root)\)", "", regex=True)
    donors = ["Os_ACT1", "At_ACT1", "CaMV35S", "ZmUBI1", "Pp_ACT1"]
    names = {"Os_ACT1": "OsACT1", "At_ACT1": "AtACT1", "CaMV35S": "CaMV 35S",
             "ZmUBI1": "ZmUBI1", "Pp_ACT1": "PpACT1"}
    sps = sorted(d.sp.unique())

    fig, axes = plt.subplots(len(donors), len(REGIONS),
                             figsize=(5.2 * len(REGIONS), 2.9 * len(donors)),
                             squeeze=False)
    rows = []
    for i, dn in enumerate(donors):
        for j, reg in enumerate(REGIONS):
            ax = axes[i][j]
            s = d[(d.donor == dn) & (d.region == reg)]
            if s.empty:
                ax.axis("off")
                continue
            v = {}
            for sp in sps:
                for tis in ("leaf", "root"):
                    v[(sp, tis)] = (s[(s.sp == sp) & (s.tissue == tis)]
                                    .set_index("target").delta.sort_index().to_numpy() * 100)
            r = bars(ax, sps, ["leaf", "root"], v, [NAVY, LBLUE],
                     "change in P (%)", rot=90)
            rows += [(dn, reg) + t for t in r]
            ax.set_title(f"{names[dn]} {reg}", loc="left", fontsize=8, fontweight="bold")
            if i or j:
                ax.set_ylabel("")
    axes[0][0].legend(handles=[Line2D([], [], marker="s", ls="", ms=7, mfc=c, mec="black",
                                      label=l) for c, l in ((NAVY, "leaf model"),
                                                            (LBLUE, "root model"))],
                      frameon=False, loc="upper left")
    fig.tight_layout()
    png = os.path.join(out_dir, f"SuppFig5_swaps_allmodels_bars_{tag}.png")
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)

    st = pd.DataFrame(rows, columns=["donor", "region", "species_or_set", "tissue", "n",
                                     "mean_change_pct", "sd", "t", "p"])
    st.to_csv(png.replace(".png", "_stats.tsv"), sep="\t", index=False)
    with open(png.replace(".png", "_legend.txt"), "w") as fh:
        fh.write(f"Supplementary Figure 5b-f in the published bar format, extended to every "
                 f"{tag} deepCRE model\n")
        fh.write("Each bar is the mean change in predicted probability of high expression over "
                 "the same nine Arabidopsis thaliana genes used in the published figure (n = 9); "
                 "error bars are the standard deviation and points the individual genes. Rows: "
                 "donor promoter; columns: donor region; navy leaf model, light blue root model. "
                 "Statistics: two-sided one-sample t-test of the nine per-gene changes against "
                 f"zero, reported per bar in {os.path.basename(png).replace('.png', '_stats.tsv')}"
                 f" ({len(st)} bars).\n")
    return png, st


# ── Supplementary Figure 6 ─────────────────────────────────────────────────
def suppfig6(allt, per, summ, extra6f, d5, sal, gen, out_dir, tag):
    pre = "2024 " if gen == "2024" else "r3 "
    at_leaf = f"{pre}Arabidopsis thaliana (leaf)"
    at_root = f"{pre}Arabidopsis thaliana (root)"
    stats = {}

    fig = plt.figure(figsize=(13.5, 13.0))
    gs = fig.add_gridspec(4, 2, height_ratios=[1.25, 1.0, 1.0, 1.25],
                          hspace=0.62, wspace=0.30)

    # a) the 672 haplotypes
    ax = fig.add_subplot(gs[0, 0])
    dd = per[per.model == at_leaf].sort_values("P").reset_index(drop=True)
    is_w = dd.group == "Hap-W box"
    ax.scatter(dd.index[~is_w], dd.P[~is_w], s=5, color=NAVY, linewidths=0, alpha=0.75)
    ax.scatter(dd.index[is_w], dd.P[is_w], s=16, color=ORANGE, edgecolor="black",
               linewidths=0.3, zorder=3)
    tib = dd.index[dd.accession == "Tibet"][0]
    ax.scatter([tib], [dd.P[tib]], s=46, facecolor="none", edgecolor=PURPLE, lw=1.0, zorder=4)
    ax.annotate("Tibet", (tib, dd.P[tib]), textcoords="offset points", xytext=(-4, 8),
                ha="right", fontsize=7, color=PURPLE)
    ax.set_xlabel("RAP2.12 promoter haplotype, ranked by predicted probability")
    ax.set_ylabel("predicted probability of high expression")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=4, mfc=NAVY, mec="none",
                              label="Hap-WT box (n = 660)"),
                       Line2D([], [], marker="o", ls="", ms=5, mfc=ORANGE, mec="black",
                              mew=0.3, label="Hap-W box (n = 12)")],
              loc="lower right", frameon=False)
    ax.set_title(f"a)  {at_leaf}", loc="left", fontweight="bold")
    pw, pwt = dd.P[is_w].to_numpy(), dd.P[~is_w].to_numpy()
    t, p = ttest_ind(pwt, pw, equal_var=False)
    u, _ = mannwhitneyu(pw, pwt, alternative="two-sided")
    stats["a"] = [(at_leaf, "Hap-WT", pwt.size, pwt.mean(), pwt.std(ddof=1), np.nan, np.nan),
                  (at_leaf, "Hap-W", pw.size, pw.mean(), pw.std(ddof=1), float(t), float(p)),
                  ("rank separation AUC", "", 0, 1 - u / (pw.size * pwt.size),
                   np.nan, np.nan, np.nan),
                  ("Hap-WT below the Hap-W maximum", "", int((pwt < pw.max()).sum()),
                   np.nan, np.nan, np.nan, np.nan)]

    # b) variant sets
    ax = fig.add_subplot(gs[0, 1])
    rp = allt[allt.analysis == "rap2.12"]
    order = [v for v in rp.variant.unique() if v != "reference"]
    short = [v.replace("I-Cat-0 ", "IP-Cat-0 ").replace(" of the alignment", "")
             .replace("(W group, 17 SNPs)", "17 SNPs") for v in order]
    x = np.arange(len(order))
    for k, (m, c, lab) in enumerate(((at_leaf, NAVY, "leaf"), (at_root, LBLUE, "root"))):
        vals = [float(rp[(rp.model == m) & (rp.variant == v)].delta.iloc[0]) * 100
                if len(rp[(rp.model == m) & (rp.variant == v)]) else np.nan for v in order]
        ax.bar(x - 0.19 + 0.38 * k, vals, 0.36, color=c, edgecolor="black", lw=0.5, label=lab)
        stats.setdefault("b", []).extend(
            [(short[i], lab, 1, vals[i], np.nan, np.nan, np.nan) for i in range(len(order))])
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(short, rotation=25, ha="right")
    ax.set_ylabel("change in P vs Col-0 reference (%)")
    ax.legend(frameon=False, title="A. thaliana model")
    ax.set_title("b)", loc="left", fontweight="bold")

    # c, d) saliency for the 41 SNPs and the nine contributors
    pos = d5[d5.HAPLOTYPE.astype(str).str.startswith("I-Cat")].POS.astype(int)
    nine = d5[(d5.HAPLOTYPE.astype(str).str.startswith("I-Cat"))
              & (d5.CONTRIBUTOR.astype(str) == "yes")].POS.astype(int)
    all41 = [k for k in sal if k.startswith(f"{at_leaf}||RAP2.12 all")]
    for lab, key, marks, ttl in (("c", all41[0] if all41 else None, pos,
                                  "all alignment SNPs"),
                                 ("d", f"{at_leaf}||RAP2.12 nine contributors", nine,
                                  "nine contributors")):
        ax = fig.add_subplot(gs[1 if lab == "c" else 2, :])
        kref = f"{at_leaf}||RAP2.12 reference"
        if key and key in sal and kref in sal:
            a, b = sal[kref].sum(axis=1), sal[key].sum(axis=1)
            xs = np.arange(1000, 1520)
            ax.plot(xs, a[1000:1520], color=GREY, lw=0.8, label="Col-0 reference")
            ax.plot(xs, b[1000:1520], color=ORANGE, lw=0.8, label=ttl)
            lo = 20135013 - 1
            ax.set_xlabel("window position (bp), gTUR of RAP2.12")
            ax.set_ylabel("contribution score")
            ax.legend(frameon=False, loc="upper left")
            d_ = np.abs(b - a)
            top = np.argsort(d_[1000:1520])[-len(marks):] + 1000
            ax.scatter(top, b[top], s=12, color=PURPLE, zorder=4, linewidths=0)
            stats[lab] = [(ttl, at_leaf, len(marks),
                           float(d_[1000:1520].sum()), np.nan, np.nan, np.nan)]
        else:
            ax.axis("off")
            ax.text(0.5, 0.5, "saliency not computed", ha="center")
        ax.set_title(f"{lab})  {ttl}, {at_leaf}", loc="left", fontweight="bold")

    # e) the 250 bp alignment block for representative accessions
    ax = fig.add_subplot(gs[3, 0])
    reps = ["Col-0", "Sichuan", "Tibet", "IP-Cat-0", "IP-Her-12", "IP-Rel-0"]
    reps = [r for r in reps if r in set(per.accession)]
    pv = per[per.model == at_leaf].set_index("accession")
    ax.barh(range(len(reps)), [pv.loc[r, "P"] for r in reps],
            color=[ORANGE if pv.loc[r, "group"] == "Hap-W box" else NAVY for r in reps],
            edgecolor="black", lw=0.5)
    ax.set_yticks(range(len(reps)))
    ax.set_yticklabels([f"{r}\n{pv.loc[r, 'group']}" for r in reps], fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel("predicted probability of high expression")
    ax.set_title("e)  representative accessions", loc="left", fontweight="bold")
    stats["e"] = [(r, str(pv.loc[r, "group"]), 1, float(pv.loc[r, "P"]),
                   np.nan, np.nan, np.nan) for r in reps]

    # f) nine contributors vs the five published control substitutions
    ax = fig.add_subplot(gs[3, 1])
    f6 = extra6f[extra6f.model == at_leaf]
    real = f6[f6.variant.str.startswith("nine")]
    ctrl = f6[f6.variant.str.startswith("control")]
    rv = float(real.delta.iloc[0]) * 100 if len(real) else np.nan
    cv = ctrl.delta.to_numpy(float) * 100
    ax.bar([0], [rv], 0.6, color=ORANGE, edgecolor="black", lw=0.5)
    ax.bar(np.arange(1, len(cv) + 1), cv, 0.6, color=PEACH, edgecolor="black", lw=0.5)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(range(len(cv) + 1))
    ax.set_xticklabels(["nine\ncontributors"] + [f"control {i}" for i in range(1, len(cv) + 1)],
                       rotation=25, ha="right")
    ax.set_ylabel("change in P vs Col-0 reference (%)")
    ax.set_title(f"f)  {at_leaf}", loc="left", fontweight="bold")
    if len(cv) > 1:
        t, p = ttest_1samp(cv, rv)
        stats["f"] = [("nine contributors", at_leaf, 1, rv, np.nan, np.nan, np.nan),
                      ("five control sets", at_leaf, len(cv), float(cv.mean()),
                       float(cv.std(ddof=1)), float(t), float(p))]

    png = os.path.join(out_dir, f"SuppFig6_full_{tag}.png")
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)

    with open(png.replace(".png", "_legend.txt"), "w") as fh:
        fh.write(f"Supplementary Figure 6, rebuilt in full, {tag} models (20261109 rerun). "
                 f"Every panel uses the Arabidopsis thaliana leaf model, as the published "
                 f"figure and Supplementary Note 3 do; panel b additionally shows the root "
                 f"model for comparison.\n\n")
        fh.write("a) Predicted probabilities of the 672 RAP2.12 promoter haplotypes of Lou et "
                 "al. (2022), ranked. Navy Hap-WT box (n = 660), orange Hap-W box (n = 12); the "
                 "Tibet accession, a Hap-W haplotype, is circled. Two-sided Welch's t-test "
                 "between the group means, and the Mann-Whitney AUC as a rank-based separation "
                 "measure that the 660-vs-12 imbalance cannot inflate.\n")
        fh.write("b) Change in predicted probability for each variant set relative to the Col-0 "
                 "reference sequence (single sequences, no replication, so no test).\n")
        fh.write("c, d) DeepSHAP contribution scores summed over the four bases across the gTUR "
                 "region, reference against the variant sequence; purple points mark the "
                 "positions of largest change.\n")
        fh.write("e) Predicted probabilities of representative accessions of the two haplotype "
                 "groups. The transcription-factor binding sites of the published panel come "
                 "from PlantPAN 4.0 and are not reproducible from the archived files, so that "
                 "track is not redrawn here.\n")
        fh.write("f) The nine IP-Cat-0 contributor SNPs against the five control substitution "
                 "sets of Supplementary Data 5 (columns ALTCtrMUT-1 to -5), which replace the "
                 "same nine positions with different bases. Two-sided one-sample t-test of the "
                 "five controls against the value of the real nine.\n\n")
        for lab in ("a", "b", "c", "d", "e", "f"):
            if lab not in stats:
                continue
            fh.write(f"panel {lab})\n")
            for g_, s_, n_, m_, sd_, t_, p_ in stats[lab]:
                line = f"  {g_:34s} {s_:30s} n={n_:3d}  value {m_:+9.4f}"
                if not np.isnan(sd_):
                    line += f"  SD {sd_:7.4f}"
                if not np.isnan(t_):
                    line += f"  t = {t_:+8.3f}  p = {p_:.3g}  {star(p_)}"
                fh.write(line + "\n")
            fh.write("\n")
    return png


def read_fasta(path):
    d, n = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            n = line[1:]
            d[n] = ""
        elif n:
            d[n] += line.upper()
    return d


def main(res_dir, scratch, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    sw = pd.read_csv(os.path.join(res_dir, "supp_allmodels_20261109_swaps_long.tsv"), sep="\t")
    allt = pd.read_csv(os.path.join(res_dir, "supp_allmodels_20261109_all.tsv"), sep="\t")
    per = pd.read_csv(os.path.join(res_dir, "lou672_20261109_per_accession.tsv"), sep="\t")
    summ = pd.read_csv(os.path.join(res_dir, "lou672_20261109_summary.tsv"), sep="\t")
    extra = pd.concat([
        pd.read_csv(os.path.join(scratch, "suppextra_20261109_predictions.tsv"), sep="\t"),
        pd.read_csv(os.path.join(scratch, "panel5n_block_20261109.tsv"), sep="\t"),
    ], ignore_index=True)
    extra6f = pd.read_csv(os.path.join(scratch, "panel6f_published_20261109.tsv"), sep="\t")
    meta = pd.read_csv(os.path.join(scratch, "suppextra_20261109_meta.tsv"), sep="\t")
    d5 = pd.read_csv(os.path.join(scratch, "suppdata5_published.tsv"), sep="\t")
    sal = dict(np.load(os.path.join(scratch, "suppextra_20261109_saliency.npz"),
                       allow_pickle=True))

    fa = read_fasta(os.path.join(
        scratch, "tutorial/File_2_deepCRE_tutorial_promoter_sequences_"
                 "Supplementary-file-2_dCRE-Pelekeetal2025.txt"))
    zm = read_fasta(os.path.join(scratch, "..", "zmubi1", "ZmUBI1_donor.fa"))
    donors = {
        "Os_ACT1": fa["gTUR_Osativa_OsACT1_KP100426-PIG2_5UTR500BP"],
        "At_ACT1": fa["gTUR_AthalianaTARI10_ACT1-AT2G37620_5UTR500BP"],
        "Pp_ACT1": fa["gTUR_Ppatens_PpACT1_Physcomitrium-patens-v6.1_5UTR500BP"],
        "ZmUBI1": zm["gTUR_Zmays_ZmUBI1_NC_050100.1_5UTR500BP"],
    }

    for gen, tag in (("round-3", "round3models"), ("2024", "2024models")):
        print("written", suppfig5(sw, extra, sal, meta, donors, gen, out_dir, tag))
        print("written", swaps_allmodels_bars(sw, extra, gen, out_dir, tag)[0])
        print("written", suppfig6(allt, per, summ, extra6f, d5, sal, gen, out_dir, tag))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
