"""
Purpose:
    Final rebuild script rendering Supplementary Figure 5 (`5a–o`: window schematic, cross-species
    swaps incl. *ZmUBI1*, *MEF1* DeepSHAP saliency tracks, published *OsACT1*/*AtACT1*/*PpACT1*
    alignment in `5n`, and fragment-only box plots in `5o`), panel bar summaries, and companion
    statistics/legend `.txt` files across both 2024 and 2026 model generations.

Manuscript (Peleke et al., 2026):
    Generates Supplementary Figure 5 offline using pre-computed tables in `code/05_supplementary_figures/`.

Usage & I/O (Python 3.11 + `requirements_analysis.txt`):
    python3 code/05_supplementary_figures/make_supp_figures_full.15.py \
        <results_dir> <scratch_dir> <out_dir>
    - Expects `suppextra_20261109_*`, `aligned_20261109_*`, `fragment_20261109_*`, `peak_20261109_*`,
      `lou672_20261109_*`, `supp_allmodels_20261109_*`, and `tutorial/` files.

How to change or extend:
    Adjust species order (`FOUR`, `SHORT`), palette constants (`SITE_*`), or figure dimensions
    at the top of the script and pass a custom `<out_dir>`.
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

# Palette taken from the deepCRE web tool itself, so the figures and the site
# read as one system:
#   app.py COLOR_PALETTE = ["#4F1787", "#EB3678"]   High / Low expression
#   lib/ui/mutation.py   palette = ["#6b7280", "#33BBC5"]  original / mutated
#   lib/ui/saliency.py   base_scale  green / cornflowerblue / darkorange / red
#   lib/ui/saliency.py   _REGION_SPANS  gUR grey, gTUR red, gTDR cornflowerblue,
#                                       gDR grey
#   app.py CSS accent gradient #FF4B4B -> #FF8E53
SITE_HIGH = "#4F1787"      # deep purple, "High"
SITE_LOW = "#EB3678"       # pink, "Low"
SITE_MUT = "#33BBC5"       # teal, "mutated"
SITE_ORIG = "#6b7280"      # grey, "original"
SITE_ACCENT = "#FF8E53"    # button gradient, used here for callouts
SITE_HIGH_L = "#B9A3D4"    # light tints for the paired light/dark series
SITE_LOW_L = "#F6A6C1"
SITE_MUT_L = "#A3DEE3"
SITE_PEAK = "#E0A526"      # mustard, for the aligned high-contribution regions

NAVY, ORANGE, MAUVE, PURPLE = SITE_HIGH, SITE_LOW, SITE_MUT, SITE_ACCENT
PEACH, ROSE, LORANGE, LBLUE = SITE_LOW_L, SITE_LOW_L, SITE_ACCENT, SITE_HIGH_L
GREY = SITE_ORIG

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "savefig.dpi": 350, "figure.dpi": 110,
})

# ordered by relatedness from Arabidopsis outwards: Brassicaceae, Caryophyllales,
# asterids, then the grasses (BOP before PACMAD)
FOUR = ["Arabidopsis thaliana", "Camelina sativa", "Beta vulgaris", "Solanum lycopersicum",
        "Oryza sativa", "Sorghum bicolor", "Zea mays"]
SHORT = ["A. thaliana", "C. sativa", "B. vulgaris", "S. lycopersicum",
         "O. sativa", "S. bicolor", "Z. mays"]
REGIONS = ["gUR", "gTUR", "gUR+gTUR"]
REG_COL = {"gUR": NAVY, "gTUR": ORANGE, "gUR+gTUR": MAUVE}
MEF1 = "AT5G52630"
# which model each donor's saliency panel is shown with
SAL_MODEL_OF = {"Os_ACT1": "PREZea mays (leaf)",
                "At_ACT1": "PREZea mays (leaf)",
                "Pp_ACT1": "PREZea mays (leaf)",
                "ZmUBI1": "PREZea mays (leaf)"}
# the web tool's own base scale (lib/ui/saliency.py base_scale)
BASE_COL = {"A": "#008000", "C": "#6495ED", "G": "#FF8C00", "T": "#FF0000",
            "N": GREY, "-": "#D1D5DB"}
CONS_LO, CONS_HI = 1185, 1230


def disp(m):
    """Model name as shown inside the plot; the round-3 prefix is dropped."""
    return m.replace("r3 ", "")


def star(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s."


def bars(ax, groups, series, vals, colors, ylab, rot=35):
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
    ax.set_xticklabels(groups, rotation=rot, ha="right" if rot else "center",
                       fontsize=6.3)
    ax.set_ylabel(ylab)
    ax.set_xlim(-0.6, ng - 0.4)
    return rows


# ── panel 5a: the extraction window ────────────────────────────────────────
def panel_window(ax):
    # grey tones; the web tool's region hues are not needed to read the schematic
    segs = [("gUR", 1, 1000, "#374151"), ("gTUR", 1001, 1500, "#6b7280"),
            ("N pad", 1501, 1520, "#e5e7eb"), ("gTDR", 1521, 2020, "#9ca3af"),
            ("gDR", 2021, 3020, "#cbd0d7")]
    for name, a, b, c in segs:
        ax.add_patch(Rectangle((a, 0.55), b - a + 1, 0.3, facecolor=c,
                               edgecolor="black", linewidth=0.6))
        lum = sum(int(c[i:i + 2], 16) for i in (1, 3, 5)) / 3
        ax.text((a + b) / 2, 0.70, name, ha="center", va="center",
                color="black" if lum > 150 else "white",
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
def suppfig5(sw, extra, sal, aln, frag, peaks, pk, gen, out_dir, tag):
    """Panels a-o. ZmUBI1 sits in sequence as g), so the published g-n shift to h-o."""
    pre = "2024 " if gen == "2024" else "r3 "
    d = sw[sw.generation == gen]
    e = extra[extra.generation == gen]
    stats, donor_peaks = {}, {}

    def pull(donor, region, species, tissue, src=None):
        s = (src if src is not None else d)
        s = s[(s.variant == f"{donor} {region}") & (s.model == f"{pre}{species} ({tissue})")]
        return s.set_index("target").delta.sort_index().to_numpy() * 100

    fig = plt.figure(figsize=(14.0, 17.0))
    gs = fig.add_gridspec(7, 12, height_ratios=[0.7, 0.10, 1.35, 1.35, 1.05, 1.05, 1.3],
                          hspace=0.85, wspace=1.25)
    C3 = [slice(0, 4), slice(4, 8), slice(8, 12)]       # three panels across
    C4 = [slice(0, 3), slice(3, 6), slice(6, 9), slice(9, 12)]

    ax = fig.add_subplot(gs[0, :])
    panel_window(ax)
    ax.set_title("a)", loc="left", fontweight="bold")

    # ── b, c: OsACT1 leaf / root ; d: AtACT1 - OsACT1 leaf ──────────────────
    for col, tis, lab in ((0, "leaf", "b"), (1, "root", "c")):
        ax = fig.add_subplot(gs[2, C3[col]])
        v = {(sh, r): pull("Os_ACT1", r, sp, tis) for sp, sh in zip(FOUR, SHORT) for r in REGIONS}
        stats[lab] = bars(ax, SHORT, REGIONS, v, [REG_COL[r] for r in REGIONS],
                          f"OsACT1 swap, change in P (%), {tis}")
        ax.set_title(f"{lab})", loc="left", fontweight="bold")

    ax = fig.add_subplot(gs[2, C3[2]])
    v = {(sh, r): pull("At_ACT1", r, sp, "leaf") - pull("Os_ACT1", r, sp, "leaf")
         for sp, sh in zip(FOUR, SHORT) for r in REGIONS}
    stats["d"] = bars(ax, SHORT, REGIONS, v, [REG_COL[r] for r in REGIONS],
                      "AtACT1 minus OsACT1, change in P (%), leaf")
    ax.set_title("d)", loc="left", fontweight="bold")

    ax = fig.add_subplot(gs[1, :])
    ax.axis("off")
    ax.legend(handles=[Line2D([], [], marker="s", ls="", ms=8, mfc=REG_COL[r], mec="black",
                              mew=0.4, label=f"donor region {r}") for r in REGIONS]
                     + [Line2D([], [], marker="o", ls="", ms=3.5, color="black",
                               label="individual gene (n = 9)")],
              ncol=4, loc="center", frameon=False, fontsize=8,
              handletextpad=0.6, columnspacing=2.2)

    # ── e: PpACT1 gTUR ; f: CaMV 35S gUR ; g: ZmUBI1 gUR ────────────────────
    for col, (donor, region, lab, cols, src) in enumerate((
            ("Pp_ACT1", "gTUR", "e", (SITE_HIGH, SITE_HIGH_L), None),
            ("CaMV35S", "gUR", "f", (SITE_LOW, SITE_LOW_L), None),
            ("ZmUBI1", "gTUR", "g", (SITE_MUT, SITE_MUT_L), e))):
        ax = fig.add_subplot(gs[3, C3[col]])
        v = {(sh, t): pull(donor, region, sp, t, src=src)
             for sp, sh in zip(FOUR, SHORT) for t in ("leaf", "root")}
        name = {"Pp_ACT1": "PpACT1", "CaMV35S": "CaMV 35S", "ZmUBI1": "ZmUBI1"}[donor]
        stats[lab] = bars(ax, SHORT, ["leaf", "root"], v, list(cols),
                          f"{name} {region} swap, change in P (%)")
        ax.set_title(f"{lab})  {name} {region}", loc="left", fontweight="bold")
        ax.legend(handles=[Line2D([], [], marker="s", ls="", ms=6, mfc=c, mec="black",
                                  mew=0.4, label=l)
                           for c, l in zip(cols, ("leaf model", "root model"))],
                  frameon=False, fontsize=6, loc="upper left",
                  handletextpad=0.4, labelspacing=0.25)

    # ── h-m: MEF1 probability and saliency for the three published donors ───
    sal_models = [f"{pre}Arabidopsis thaliana (leaf)", f"{pre}Zea mays (leaf)"]
    panels = [("Os_ACT1", "OsACT1", "h", "i"), ("At_ACT1", "AtACT1", "j", "k"),
              ("Pp_ACT1", "PpACT1", "l", "m")]
    for (row, col0), (dk, dn, lp, ls_) in zip(((4, 0), (4, 2), (5, 0)), panels):
        axp = fig.add_subplot(gs[row, C4[col0]])
        mm = SAL_MODEL_OF[dk].replace("PRE", pre)
        base = sw[(sw.model == mm) & (sw.target == MEF1) & (sw.variant == "Os_ACT1 gUR")]
        nat = float((base.P - base.delta).iloc[0])
        s = d[(d.model == mm) & (d.target == MEF1) & (d.variant == f"{dk} gTUR")]
        swp = float(s.P.iloc[0])
        axp.bar([0], [nat], 0.52, color=SITE_ORIG, edgecolor="black", lw=0.5, label="native")
        axp.bar([1], [swp], 0.52, color=SITE_MUT, edgecolor="black", lw=0.5,
                label=f"{dn} gTUR swap")
        axp.axhline(0.5, color="black", lw=0.6, ls="--")
        axp.set_xticks([0, 1])
        axp.set_xticklabels(["native", f"{dn}\ngTUR swap"])
        axp.set_xlim(-0.6, 1.6)
        axp.set_ylabel(f"MEF1 ({MEF1})\npredicted probability")
        axp.set_ylim(0, 1)
        axp.set_title(f"{lp})  {disp(mm)}", loc="left", fontweight="bold", fontsize=7.5)
        stats[lp] = [("native", mm, 1, nat, np.nan, np.nan, np.nan),
                     (f"{dn} gTUR swap", mm, 1, swp, np.nan, np.nan, np.nan),
                     ("change", mm, 1, swp - nat, np.nan, np.nan, np.nan)]

        axs = fig.add_subplot(gs[row, C4[col0 + 1]])
        m = SAL_MODEL_OF[dk].replace("PRE", pre)
        kn, ks = f"{m}||MEF1 native", f"{m}||MEF1 {dk} gTUR"
        a, b = sal[kn].sum(axis=1), sal[ks].sum(axis=1)
        xx = np.arange(950, 1560)
        axs.plot(xx, a[950:1560], color=SITE_ORIG, lw=0.7, label="native")
        axs.plot(xx, b[950:1560], color=SITE_MUT, lw=0.7, label=f"{dn} gTUR")
        g = np.abs(b[1000:1500])
        cs = np.convolve(g, np.ones(CONS_HI - CONS_LO + 1), "valid")
        i = int(np.argmax(cs))
        plo, phi = 1001 + i, 1000 + i + (CONS_HI - CONS_LO + 1)
        donor_peaks[dn] = (plo, phi, m)
        axs.axvspan(plo, phi, color=SITE_PEAK, alpha=0.32, lw=0)
        axs.axvline(1000.5, color="black", lw=0.6, ls=":")
        axs.set_xlabel("window position (bp)")
        axs.set_ylabel("contribution score")
        axs.legend(frameon=False, fontsize=6, loc="upper left")
        axs.set_title(f"{ls_})  {dn} peak {plo}-{phi}\n{disp(m)}", loc="left",
                      fontweight="bold", fontsize=7.0)

    ax = fig.add_subplot(gs[5, 6:])
    ax.axis("off")
    ax.legend(handles=[Line2D([], [], marker="s", ls="", ms=8, mfc=SITE_ORIG, mec="black",
                              mew=0.4, label="native sequence"),
                       Line2D([], [], marker="s", ls="", ms=8, mfc=SITE_MUT, mec="black",
                              mew=0.4, label="donor gTUR swap (saliency)"),
                       Line2D([], [], marker="s", ls="", ms=8, mfc=SITE_PEAK, mec="none",
                              alpha=0.55, label="aligned high contribution regions"),
                       Line2D([], [], ls=":", color="black", lw=1.0,
                              label="gUR / gTUR boundary")],
              loc="center left", frameon=False, fontsize=7.5,
              handletextpad=0.6, labelspacing=0.5)

    # ── n: the published three-donor alignment, Arabidopsis leaf ────────────
    axm = fig.add_subplot(gs[6, 0:6])
    m = f"{pre}Zea mays (leaf)"
    LO, HI = CONS_LO, CONS_HI
    names = ["OsACT1", "AtACT1", "PpACT1"]
    cols = [i for i in range(len(aln)) if aln.os_pos[i] and LO <= aln.os_pos[i] <= HI]
    cols = list(range(min(cols), max(cols) + 1))
    cons_cols = [i for i in cols if aln.conserved[i]]
    cons_pos = [aln.os_pos[i] for i in cons_cols]
    for r, nm in enumerate(names):
        for j, i in enumerate(cols):
            bch = aln.rows[nm][i]
            axm.text(j, -r, bch, ha="center", va="center", fontsize=6.0,
                     color=BASE_COL.get(bch, GREY),
                     fontweight="bold" if i in cons_cols else "normal",
                     alpha=1.0 if i in cons_cols else 0.40)
    for i in cons_cols:
        axm.axvspan(cols.index(i) - 0.5, cols.index(i) + 0.5, ymin=0.30,
                    color=SITE_ACCENT, alpha=0.18, lw=0)
    tr = sal[f"{m}||MEF1 Os_ACT1 gTUR"].sum(axis=1)
    sc = np.array([tr[aln.os_pos[i] - 1] if aln.os_pos[i] else 0.0 for i in cols])
    sc = sc / (np.abs(sc).max() or 1)
    axm.bar(np.arange(len(sc)), sc * 0.9, 0.8, bottom=-len(names) - 0.4,
            color=[SITE_MUT if v > 0 else SITE_HIGH for v in sc], linewidth=0)
    axm.set_yticks(list(range(0, -len(names), -1)) + [-len(names) - 0.4])
    axm.set_yticklabels(names + [f"contribution\n{disp(m)}"])
    tick = [j for j, i in enumerate(cols) if aln.os_pos[i] and aln.os_pos[i] % 5 == 0]
    axm.set_xticks(tick)
    axm.set_xticklabels([aln.os_pos[cols[j]] for j in tick])
    axm.set_xlabel("window position in OsACT1 (bp)")
    axm.set_xlim(-1, len(cols))
    axm.set_ylim(-len(names) - 1.6, 0.8)
    for s_ in axm.spines.values():
        s_.set_visible(False)
    n_all = aln.n_conserved
    obs, exp, bp_p = aln.enrichment(LO, HI)
    axm.set_title(f"n)  donor alignment, window {LO}-{HI}",
                  loc="left", fontweight="bold", fontsize=8)

    # ── o: the fragment-only experiment, Arabidopsis leaf and root ──────────
    axn = fig.add_subplot(gs[6, 6:])
    fr = frag[frag.generation == gen]

    def fseries(mm, var):
        s_ = fr[(fr.model == mm) & (fr.variant == var)]
        return s_.set_index("target").delta.sort_index()

    def fcontrol(mm):
        s_ = fr[(fr.model == mm) & (fr.variant.str.startswith("control"))]
        return s_.groupby("target").delta.mean().sort_index()

    at_l, at_r = f"{pre}Zea mays (leaf)", f"{pre}Zea mays (root)"
    labels, arrs, colors = [], [], []
    stats["o"] = []
    for lab, mm in (("leaf", at_l), ("root", at_r)):
        intact, degen, ctrl = (fseries(mm, f"fragment {CONS_LO}-{CONS_HI} intact"),
                               fseries(mm, "conserved sites degenerated"), fcontrol(mm))
        common = intact.index.intersection(degen.index).intersection(ctrl.index)
        I, D, C = (intact.loc[common].to_numpy() * 100, degen.loc[common].to_numpy() * 100,
                   ctrl.loc[common].to_numpy() * 100)
        for nm, v, c in (("fragment\nintact", I, SITE_ORIG),
                         ("conserved sites\ndegenerated", D, SITE_MUT)):
            labels.append(f"{lab}\n{nm}")
            arrs.append(v)
            colors.append(c)
        stats["o"].append((f"{lab}: fragment intact", mm, I.size, float(I.mean()),
                           float(I.std(ddof=1)), np.nan, np.nan))
        for nm, v in (("conserved sites degenerated", D), ("control sites degenerated", C)):
            t, p = ttest_rel(v, I)
            stats["o"].append((f"{lab}: {nm}", mm, v.size, float(v.mean()),
                               float(v.std(ddof=1)), float(t), float(p)))
        t, p = ttest_rel(D, C)
        stats["o"].append((f"{lab}: degenerated vs control", mm, D.size,
                           float((D - C).mean()), float((D - C).std(ddof=1)),
                           float(t), float(p)))

    bp = axn.boxplot(arrs, widths=0.6, patch_artist=True,
                     medianprops=dict(color="black", lw=1.0),
                     flierprops=dict(marker="", ls="none"))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_edgecolor("black")
        patch.set_linewidth(0.6)
    rng = np.random.default_rng(3)
    for i, v in enumerate(arrs, start=1):
        axn.scatter(i + rng.uniform(-0.11, 0.11, v.size), v, s=9, color="black",
                    alpha=0.7, linewidths=0, zorder=3)
    axn.axhline(0, color="black", lw=0.6)
    axn.axvline(2.5, color=GREY, lw=0.6, ls=":")
    axn.set_xticks(range(1, len(arrs) + 1))
    axn.set_xticklabels(labels, fontsize=6.2)
    axn.set_ylabel("change in P vs the native gene (%)")

    # significance brackets: intact vs degenerated, and degenerated vs control
    top = max(v.max() for v in arrs)
    span = top - min(v.min() for v in arrs)
    for k, base_i in enumerate((0, 2)):
        pv = [r for r in stats["o"]
              if r[0].endswith("conserved sites degenerated")][k][6]
        x1, x2 = base_i + 1, base_i + 2
        y = top + span * 0.09
        axn.plot([x1, x1, x2, x2], [y, y + span * 0.03, y + span * 0.03, y],
                 lw=0.7, color="black")
        axn.text((x1 + x2) / 2, y + span * 0.04, f"{star(pv)}  p = {pv:.3g}",
                 ha="center", va="bottom", fontsize=6.5)
    axn.set_ylim(min(v.min() for v in arrs) - span * 0.08, top + span * 0.30)
    axn.set_title("o)", loc="left", fontweight="bold")

    png = os.path.join(out_dir, f"SuppFig5_full_{tag}.png")
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)

    with open(png.replace(".png", "_legend.txt"), "w") as fh:
        fh.write(f"Supplementary Figure 5, rebuilt in full, {tag} models (20261109 rerun). "
                 f"ZmUBI1 is a donor added in this rerun and sits in sequence as panel g), so "
                 f"the published panels g-n appear here as h-o.\n\n")
        fh.write("a) Segmentation of the 3020 bp deepCRE input window; gUR and gTUR are the "
                 "segments replaced by donor promoter sequence. Donor coordinates were checked "
                 "against the genes they come from: the AtACT1 donor is identical to window "
                 "positions 1-1500 of AT2G37620 extracted by the toolkit itself (0 mismatches), "
                 "the CaMV 35S donor is the 1500 bp immediately upstream of the Hyg(R) start in "
                 "pGZ14.1021 (GenBank KP100426, exact), and for every donor gURgTUR[1:1000] "
                 "equals its gUR file and gURgTUR[1001:1500] its gTUR file.\n")
        fh.write("b-g) Species run left to right by relatedness from Arabidopsis "
                 "thaliana: Brassicaceae, Caryophyllales, asterids, then the grasses. "
                 "C. sativa, B. vulgaris and O. sativa are added in this rerun; the published "
                 "figure had only A. thaliana, S. lycopersicum, S. bicolor and Z. mays. "
                 "Nine Arabidopsis thaliana genes per bar (n = 9); bars are the mean change "
                 "in predicted probability of high expression, error bars the standard "
                 "deviation, points individual genes. Two-sided one-sample t-test of the nine "
                 "per-gene changes against zero.\n")
        fh.write("g) ZmUBI1 gTUR: the maize polyubiquitin-1 locus of Zea mays "
                 "GCF_902167145.1, NC_050100.1 minus strand; the gTUR is the first 500 bp of "
                 "the transcript, i.e. the 5' UTR exon plus the start of intron 1. Its gUR and "
                 "gUR+gTUR swaps are in SuppFig5_swaps_allmodels_bars_*.png.\n")
        fh.write("h, j, l) Predicted probability of MEF1 (AT5G52630) natively and after the "
                 "donor gTUR swap. i, k, m) DeepSHAP contribution scores summed over the four "
                 "bases, window positions 950-1560; the dotted line is the gUR/gTUR boundary "
                 "and the mustard band the aligned high contribution region, the 46 bp "
                 "window of largest summed |contribution| inside that donor's gTUR, each shown "
                 "with the model named in the panel title. Regions: "
                 + ", ".join(f"{k} {v[0]}-{v[1]} ({v[2]})" for k, v in sorted(donor_peaks.items()))
                 + ".\n")
        fh.write(f"n) The published OsACT1/AtACT1/PpACT1 alignment (Supplementary file 3, "
                 f"{aln.n_cols} columns; each donor's ungapped {aln.n_ref} bp occupies window "
                 f"positions 1151-1400 of its own gTUR), shown over the block {LO}-{HI}. A dash "
                 f"is a gap. Bold, highlighted columns are identical in all three donors: {obs} "
                 f"of the {n_all} such columns in the whole alignment fall in this 46 bp block, "
                 f"against {exp:.1f} expected (one-sided binomial p = {bp_p:.3g}), and no 46 bp "
                 f"window of the aligned region holds more. The bar track is the DeepSHAP "
                 f"contribution score of the OsACT1 gTUR swap under {m}, scaled to its own "
                 f"maximum.\n")
        fh.write(f"o) Only window positions {CONS_LO}-{CONS_HI} of the OsACT1 gTUR — the 46 bp "
                 f"block of panel n — are pasted into the same positions of each test gene's "
                 f"own window; everything else stays native, so the change is attributable to "
                 f"those 46 bases. Three constructs per gene: the fragment intact, the fragment "
                 f"with all {len(cons_pos)} sites conserved across OsACT1/AtACT1/PpACT1 "
                 f"transverted (A<->C, G<->T). Nine genes; the leaf and the root model named "
                 f"below. Brackets give the two-sided paired t-test of the intact fragment "
                 f"against the conserved-site degeneration. A control transverting the same "
                 f"number of NON-conserved sites of the same fragment (mean of 20 seeded draws "
                 f"per gene; drawn sites in fragment_20261109_control_sites.tsv) was run "
                 f"alongside and is reported below but no longer plotted.\n\n")
        for lab in ("b", "c", "d", "e", "f", "g", "h", "j", "l", "o"):
            if lab not in stats:
                continue
            fh.write(f"panel {lab})\n")
            for g_, s_, n_, m_, sd_, t_, p_ in stats[lab]:
                line = f"  {g_:34s} {s_:34s} n={n_:2d}  mean {m_:+8.3f}"
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


class DonorAlignment:
    """The published OsACT1/AtACT1/PpACT1 alignment (Supplementary file 3).

    339 columns; each donor's ungapped 250 bp sits at offset 150 of its own gTUR,
    i.e. window positions 1151-1400. os_pos[i] is the OsACT1 window position of
    alignment column i, or None where OsACT1 carries a gap.
    """

    ORDER = ["OsACT1", "AtACT1", "PpACT1"]

    def __init__(self, path, donor_gtur, ref_offset=150):
        raw, name = {}, None
        for line in open(path):
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].replace("250", "")
                raw[name] = ""
            elif name:
                raw[name] += line.upper()
        self.rows = {k: raw[k] for k in self.ORDER}
        n = {len(v) for v in self.rows.values()}
        assert len(n) == 1, f"ragged alignment: {n}"
        for k in self.ORDER:
            assert donor_gtur[k].find(self.rows[k].replace("-", "")) == ref_offset, (
                f"{k} does not sit at gTUR offset {ref_offset}")
        self.n_cols = n.pop()
        self.os_pos, j = [], 0
        for c in self.rows["OsACT1"]:
            if c == "-":
                self.os_pos.append(None)
            else:
                self.os_pos.append(1001 + ref_offset + j)
                j += 1
        self.n_ref = j
        self.conserved = [
            self.rows["OsACT1"][i] == self.rows["AtACT1"][i] == self.rows["PpACT1"][i]
            != "-" for i in range(self.n_cols)]
        self.n_conserved = sum(self.conserved)

    def __len__(self):
        return self.n_cols

    def enrichment(self, lo, hi):
        from scipy.stats import binomtest
        obs = sum(1 for i in range(self.n_cols)
                  if self.conserved[i] and self.os_pos[i] and lo <= self.os_pos[i] <= hi)
        frac = (hi - lo + 1) / self.n_ref
        return obs, self.n_conserved * frac, binomtest(
            obs, self.n_conserved, frac, alternative="greater").pvalue


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
        pd.read_csv(os.path.join(scratch, "aligned_20261109_predictions.tsv"), sep="\t"),
    ], ignore_index=True)
    pk = pd.read_csv(os.path.join(scratch, "peak_20261109_predictions.tsv"), sep="\t")
    pkt = pd.read_csv(os.path.join(scratch, "peak_20261109_peaks.tsv"), sep="\t")
    peaks = {r.model: (int(r.peak_lo), int(r.peak_hi)) for _, r in pkt.iterrows()}
    extra6f = pd.read_csv(os.path.join(scratch, "panel6f_published_20261109.tsv"), sep="\t")
    d5 = pd.read_csv(os.path.join(scratch, "suppdata5_published.tsv"), sep="\t")
    fa = read_fasta(os.path.join(
        scratch, "tutorial/File_2_deepCRE_tutorial_promoter_sequences_"
                 "Supplementary-file-2_dCRE-Pelekeetal2025.txt"))
    donor_gtur = {
        "OsACT1": fa["gTUR_Osativa_OsACT1_KP100426-PIG2_5UTR500BP"],
        "AtACT1": fa["gTUR_AthalianaTARI10_ACT1-AT2G37620_5UTR500BP"],
        "PpACT1": fa["gTUR_Ppatens_PpACT1_Physcomitrium-patens-v6.1_5UTR500BP"]}
    aln = DonorAlignment(
        os.path.join(scratch, "donor_alignment_Supplementary-file-3.txt"), donor_gtur)
    frag = pd.read_csv(os.path.join(scratch, "fragment_20261109_predictions.tsv"), sep="\t")
    sal = dict(np.load(os.path.join(scratch, "suppextra_20261109_saliency.npz"),
                       allow_pickle=True))

    for gen, tag in (("round-3", "round3models"), ("2024", "2024models")):
        print("written", suppfig5(sw, extra, sal, aln, frag, peaks, pk, gen, out_dir, tag))
        print("written", swaps_allmodels_bars(sw, extra, gen, out_dir, tag)[0])
        print("written", suppfig6(allt, per, summ, extra6f, d5, sal, gen, out_dir, tag))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
