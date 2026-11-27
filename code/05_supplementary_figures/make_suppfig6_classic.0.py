"""
Purpose:
    Earlier Supplementary Figure 6 script rebuilding panels `6a–g` in the published layout
    (haplotype strip plot, reference vs. *IP-Cat-0* bar, saliency maps, accession SNP track,
    and 5-control comparison) across 2024 and 2026 deepCRE models.

Manuscript (Peleke et al., 2026):
    Predecessor to `make_suppfig6_classic.14.py` (which updates panel `6g` to the 1,000-draw
    permutation histogram).

Usage & I/O (Python 3.11 + `requirements_analysis.txt`):
    python3 code/05_supplementary_figures/make_suppfig6_classic.0.py \
        <results_dir> <scratch_dir> <out_dir>

How to change or extend:
    Use `make_suppfig6_classic.14.py` for the final permutation-based panel `6g` figure.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import ttest_ind, mannwhitneyu

# the deepCRE web tool's palette
SITE_HIGH, SITE_LOW = "#4F1787", "#EB3678"      # High / Low expression
SITE_MUT, SITE_ORIG = "#33BBC5", "#6b7280"      # mutated / original
SITE_ACCENT = "#FF8E53"
BASE_COL = {"A": "#008000", "C": "#6495ED", "G": "#FF8C00", "T": "#FF0000",
            "N": "#9ca3af", "-": "#D1D5DB"}
# lib/ui/saliency.py _REGION_SPANS
REGIONS = [(1, 1000, "#6b7280"), (1001, 1500, "#FF0000"),
           (1521, 2020, "#6495ED"), (2021, 3020, "#6b7280")]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "savefig.dpi": 350, "figure.dpi": 110,
})

GENE = "AT1G53910"
GENE_START = 20134843                      # verified against the published panel e labels
REPS = ["Col-0", "Sichuan", "BuckhornPass", "Utrecht", "Lan-0",
        "Tibet", "Tanz-1", "IP-Rel-0", "IP-Cat-0"]


def star(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s."


def disp(m):
    return m.replace("r3 ", "")


def win(pos):
    """genomic position -> 1-based position in the 3020 bp window (plus strand)."""
    return pos - (GENE_START - 1000)


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


def region_bands(ax):
    for a, b, c in REGIONS:
        ax.axvspan(a, b, color=c, alpha=0.10, lw=0)
    ax.text(1000, 1.02, "TSS", transform=ax.get_xaxis_transform(),
            ha="center", va="bottom", fontsize=6.5)
    ax.text(1520, 1.02, "TTS", transform=ax.get_xaxis_transform(),
            ha="center", va="bottom", fontsize=6.5)


def main(res_dir, scratch, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    per = pd.read_csv(os.path.join(res_dir, "lou672_20261109_per_accession.tsv"), sep="\t")
    allt = pd.read_csv(os.path.join(res_dir, "supp_allmodels_20261109_all.tsv"), sep="\t")
    f6 = pd.read_csv(os.path.join(scratch, "panel6f_published_20261109.tsv"), sep="\t")
    icat = pd.read_csv(os.path.join(scratch, "tutorial/Tutorial_table3_icat.csv"))
    sal = dict(np.load(os.path.join(scratch, "suppextra_20261109_saliency.npz"),
                       allow_pickle=True))
    sal.update(dict(np.load(os.path.join(scratch, "saliency_icat_20261109.npz"),
                            allow_pickle=True)))

    import gzip
    vcf = []
    with gzip.open(os.path.join(
            scratch, "tutorial/File_3_deepCRE_tutorial_rap12-2-variant_"
                     "Supplementary-file-3_dCRE_Peleketal2025.vcf.gz"), "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                vcf.append(int(line.split("\t")[1]))
    vcf_pos = sorted(set(vcf))

    icat_pos = list(icat.POS.astype(int))
    nine = list(icat[icat.CONTRIBUTOR == "yes"].POS.astype(int))

    # the Lou alignment, projected onto Col-0 window coordinates
    aln = read_fasta(os.path.join(scratch, "lou2022_data_s1.fas"))
    col = aln["Col-0"]
    keep = [i for i, c in enumerate(col) if c != "-"]
    geno = {a: "".join(aln[a][i] if aln[a][i] != "-" else "N" for i in keep) for a in aln}
    # the block's first base, verified by requiring Col-0 to carry the REF allele
    # of every IP-Cat-0 SNP at the implied offset
    first = None
    for cand in range(1050, 1300):
        if all(0 <= win(pos) - cand < len(geno["Col-0"])
               and geno["Col-0"][win(pos) - cand] == ref
               for pos, ref in zip(icat.POS.astype(int), icat.REF.astype(str)))\
                and all(0 <= win(pos) - cand < len(geno["Col-0"]) for pos in icat_pos):
            first = cand
            break
    assert first is not None, "the Lou alignment block could not be placed in the window"
    print(f"Lou alignment block: window {first}-{first + len(geno['Col-0']) - 1}; "
          f"Col-0 carries the REF allele at all {len(icat_pos)} IP-Cat-0 SNPs")
    for gen, tag in (("2024", "2024models"), ("round-3", "round3models")):
        pre = "2024 " if gen == "2024" else "r3 "
        m = f"{pre}Arabidopsis thaliana (leaf)"
        stats = {}

        fig = plt.figure(figsize=(13.5, 13.5))
        gs = fig.add_gridspec(4, 12, height_ratios=[1.15, 1.15, 1.45, 1.0],
                              hspace=0.70, wspace=1.3)

        # ── a  the 672 haplotypes ───────────────────────────────────────────
        ax = fig.add_subplot(gs[0, 0:4])
        d = per[per.model == m]
        rng = np.random.default_rng(11)
        for grp, c, s_, z in (("Hap-WT box", SITE_HIGH, 9, 2),
                              ("Hap-W box", SITE_LOW, 20, 3)):
            v = d[d.group == grp].P.to_numpy()
            ax.scatter(rng.uniform(-0.28, 0.28, v.size), v, s=s_, color=c,
                       edgecolor="black" if grp.endswith("W box") else "none",
                       linewidths=0.3, zorder=z, label=f"{grp.split('-')[1]} (n = {v.size})")
        for acc, dx, ha in (("IP-Cat-0", 0.36, "left"), ("Tibet", -0.36, "right")):
            y = float(d[d.accession == acc].P.iloc[0])
            ax.annotate(acc, (0, y), xytext=(dx, y), textcoords="data", ha=ha,
                        va="center", fontsize=6.5, color=SITE_ACCENT,
                        arrowprops=dict(arrowstyle="-", lw=0.6, color=SITE_ACCENT))
        ax.set_xlim(-0.8, 0.8)
        ax.set_xticks([])
        ax.set_ylabel("probability of high gene expression")
        ax.legend(title="haplotype", frameon=False, loc="center left", fontsize=6.5,
                  title_fontsize=6.5)
        ax.set_title(f"a)  {disp(m)}", loc="left", fontweight="bold")
        pw = d[d.group == "Hap-W box"].P.to_numpy()
        pwt = d[d.group == "Hap-WT box"].P.to_numpy()
        t, p = ttest_ind(pwt, pw, equal_var=False)
        u, _ = mannwhitneyu(pw, pwt, alternative="two-sided")
        stats["a"] = [("Hap-WT", pwt.size, pwt.mean(), pwt.min(), pwt.max(), np.nan, np.nan),
                      ("Hap-W", pw.size, pw.mean(), pw.min(), pw.max(), float(t), float(p)),
                      ("AUC", 0, 1 - u / (pw.size * pwt.size), np.nan, np.nan, np.nan, np.nan)]

        # ── b  reference against the IP-Cat-0 haplotype ─────────────────────
        ax = fig.add_subplot(gs[1, 0:4])
        rp = allt[(allt.analysis == "rap2.12") & (allt.model == m)]
        ref = float(rp[rp.variant == "reference"].P.iloc[0])
        mut = float(rp[rp.variant.str.startswith("I-Cat-0 (W group")].P.iloc[0])
        ax.bar([0], [ref], 0.55, color=SITE_ORIG, edgecolor="black", lw=0.5)
        ax.bar([1], [mut], 0.55, color=SITE_MUT, edgecolor="black", lw=0.5)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([GENE, f"{GENE}:\nIP-Cat-0 mutated"])
        ax.set_xlim(-0.6, 1.6)
        ax.set_ylabel("probability of high gene expression")
        ax.set_title("b)", loc="left", fontweight="bold")
        stats["b"] = [("reference", 1, ref, np.nan, np.nan, np.nan, np.nan),
                      ("IP-Cat-0, 17 SNPs", 1, mut, np.nan, np.nan, np.nan, np.nan),
                      ("change (%)", 1, 100 * (mut - ref), np.nan, np.nan, np.nan, np.nan)]

        # ── c, d  saliency maps ─────────────────────────────────────────────
        key_all = [k for k in sal if k.startswith(f"{m}||RAP2.12 all")]
        for lab, row, key, marks, ttl in (
                ("c", 0, key_all[0] if key_all else None, vcf_pos,
                 f"all {len(vcf_pos)} alignment SNPs"),
                ("d", 1, f"{m}||RAP2.12 nine contributors", nine, "nine key SNPs")):
            ax = fig.add_subplot(gs[row, 4:12])
            kref = f"{m}||RAP2.12 reference"
            a = sal[kref].sum(axis=1)
            b = sal[key].sum(axis=1)
            xs = np.arange(1, 3021)
            region_bands(ax)
            for pos in marks:
                ax.axvline(win(pos), color="#444444", lw=0.45, ls=(0, (2, 2)), zorder=1)
            ax.plot(xs, a, color=SITE_ORIG, lw=0.45, label=GENE, zorder=4)
            ax.plot(xs, b, color=SITE_MUT, lw=0.45, label=f"{GENE}: mutated", zorder=5)
            ax.axhline(0, color="black", lw=0.6, ls=(0, (1, 2)))
            ax.set_xlim(1, 3020)
            ax.set_xlabel("nucleotide position")
            ax.set_ylabel("saliency score")
            ax.legend(frameon=False, loc="upper right", fontsize=6.5)
            ax.set_title(f"{lab})  {ttl}", loc="left", fontweight="bold")
            stats[lab] = [(ttl, len(marks), float(np.abs(b - a).sum()),
                           np.nan, np.nan, np.nan, np.nan)]

        # ── e  the SNP track of the representative accessions ───────────────
        ax = fig.add_subplot(gs[2, :])
        reps = [a for a in REPS if a in geno]
        cols = [p for p in icat_pos if 0 <= win(p) - first < len(geno["Col-0"])]
        grp = dict(zip(per[per.model == m].accession, per[per.model == m].group))
        for r, acc in enumerate(reps):
            for j, pos in enumerate(cols):
                b_ = geno[acc][win(pos) - first]
                ref_b = geno["Col-0"][win(pos) - first]
                if b_ != ref_b:
                    ax.add_patch(plt.Rectangle((j - 0.42, -r - 0.42), 0.84, 0.84,
                                               facecolor=BASE_COL.get(b_, "#9ca3af"),
                                               edgecolor="none"))
                ax.text(j, -r, b_, ha="center", va="center", fontsize=5.4,
                        color="white" if b_ != ref_b else "#b0b4ba",
                        fontweight="bold" if b_ != ref_b else "normal")
        ax.set_yticks([-r for r in range(len(reps))])
        ax.set_yticklabels(reps, fontsize=7)
        for tick, acc in zip(ax.get_yticklabels(), reps):
            tick.set_color(SITE_LOW if grp.get(acc) == "Hap-W box" else SITE_HIGH)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels([f"{p} ({win(p)})" + ("  ***" if p in nine else "")
                            for p in cols], rotation=90, fontsize=5.6)
        for tick, p in zip(ax.get_xticklabels(), cols):
            if p in nine:
                tick.set_fontweight("bold")
        ax.set_xlim(-0.8, len(cols) - 0.2)
        ax.set_ylim(-len(reps) + 0.4, 0.7)
        for s_ in ax.spines.values():
            s_.set_visible(False)
        ax.set_title("e)  chromosome 1 position (window position); "
                     "bold with *** = one of the nine key SNPs",
                     loc="left", fontweight="bold", fontsize=7.5)

        # ── f  the nine contributors against the five controls ──────────────
        ax = fig.add_subplot(gs[3, 0:6])
        g = f6[f6.model == m]
        real = float(g[g.variant.str.startswith("nine")].delta.iloc[0]) * 100
        ctrl = g[g.variant.str.startswith("control")].delta.to_numpy(float) * 100
        ax.bar([0], [real], 0.6, color=SITE_ORIG, edgecolor="black", lw=0.5)
        ax.bar(np.arange(1, len(ctrl) + 1), ctrl, 0.6, color="#c7cacf",
               edgecolor="black", lw=0.5)
        ax.axhline(0, color="black", lw=0.6)
        ax.set_xticks(range(len(ctrl) + 1))
        ax.set_xticklabels(["Col-0\nnine SNPs"] + [f"random\nSNPs {i}"
                                                   for i in range(1, len(ctrl) + 1)],
                           fontsize=6.5)
        ax.set_ylabel("change in predicted probability (%)")
        ax.set_title("f)", loc="left", fontweight="bold")
        stats["f"] = [("nine contributors", 1, real, np.nan, np.nan, np.nan, np.nan)] + [
            (f"control {i}", 1, float(v), np.nan, np.nan, np.nan, np.nan)
            for i, v in enumerate(ctrl, 1)]

        ax = fig.add_subplot(gs[3, 6:])
        ax.axis("off")
        ax.legend(handles=[Line2D([], [], marker="s", ls="", ms=8, mfc=SITE_HIGH,
                                  mec="none", label="Hap-WT box accession"),
                           Line2D([], [], marker="s", ls="", ms=8, mfc=SITE_LOW,
                                  mec="black", mew=0.4, label="Hap-W box accession"),
                           Line2D([], [], marker="s", ls="", ms=8, mfc=SITE_ORIG,
                                  mec="black", mew=0.4, label=f"{GENE} reference"),
                           Line2D([], [], marker="s", ls="", ms=8, mfc=SITE_MUT,
                                  mec="black", mew=0.4, label=f"{GENE} mutated")]
                          + [Line2D([], [], marker="s", ls="", ms=8, mfc=BASE_COL[b],
                                    mec="none", label=f"substituted base {b}")
                             for b in "ACGT"],
                  ncol=2, loc="center", frameon=False, fontsize=7.5,
                  handletextpad=0.6, columnspacing=1.6, labelspacing=0.5)

        png = os.path.join(out_dir, f"SuppFig6_classic_{tag}.png")
        fig.savefig(png, bbox_inches="tight")
        plt.close(fig)

        with open(png.replace(".png", "_legend.txt"), "w") as fh:
            fh.write(f"Supplementary Figure 6, rebuilt in the layout of the published figure, "
                     f"{tag} (20261109 rerun). All panels use the Arabidopsis thaliana leaf "
                     f"model, as Supplementary Note 3 specifies.\n\n")
            fh.write("a) Predicted probabilities of the 672 RAP2.12 promoter haplotypes of Lou "
                     "et al. (2022) as a strip plot. Statistics: two-sided Welch's t-test "
                     "between the group means and the Mann-Whitney AUC.\n")
            fh.write(f"b) Predicted probability of the Col-0 reference window against the "
                     f"IP-Cat-0 haplotype (all 17 of its SNPs).\n")
            fh.write("c, d) DeepSHAP saliency summed over the four bases across the whole "
                     "3020 bp window, reference against the mutated sequence. Dashed lines "
                     "mark the substituted positions; the shaded bands are the regions the web "
                     "tool itself draws (gUR and gDR grey, gTUR red, gTDR blue).\n")
            fh.write("e) Bases of the representative accessions at the IP-Cat-0 SNP positions, "
                     "read from the Lou et al. alignment. A coloured cell marks a base "
                     "differing from Col-0. Accession labels are coloured by haplotype group. "
                     "The transcription-factor binding sites of the published panel come from "
                     "PlantPAN 4.0 and are in none of the archived files, so that track is not "
                     "redrawn and has to be carried over.\n")
            fh.write("f) The nine key SNPs against the five control substitution sets of "
                     "Supplementary Data 5 (ALTCtrMUT-1 to -5), which replace the same nine "
                     "positions with different bases.\n\n")
            for lab in ("a", "b", "c", "d", "f"):
                fh.write(f"panel {lab})\n")
                for row in stats.get(lab, []):
                    name, n_, v, lo, hi, t_, p_ = row
                    line = f"  {name:28s} n={n_:4d}  value {v:+10.4f}"
                    if not np.isnan(lo):
                        line += f"  range {lo:.4f}-{hi:.4f}"
                    if not np.isnan(t_):
                        line += f"  t = {t_:+8.3f}  p = {p_:.3g}  {star(p_)}"
                    fh.write(line + "\n")
                fh.write("\n")
        print("written", png)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
