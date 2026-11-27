"""
Purpose:
    Final rebuild script rendering Supplementary Figure 6 (`6a–g`) on TAIR10 windows across 2024
    and 2026 *A. thaliana* leaf models, with panel `6g` plotting the 1,000-permutation random-allele
    null distribution against the observed 9-SNP drop and `6f` displaying *RAP2.12* haplotypes and
    PlantPAN 4.0 TFBS intervals (transcribed in `TFBS` from the published figure).

Manuscript (Peleke et al., 2026):
    Produces Supplementary Figure 6 and `SuppFig6_classic_*_legend.txt`.

Usage & I/O (Python 3.11 + `requirements_analysis.txt`):
    python3 code/05_supplementary_figures/make_suppfig6_classic.14.py \
        <results_dir> <scratch_dir> <out_dir>

How to change or extend:
    Pass alternative `<results_dir>` / `<scratch_dir>` tables or edit `TFBS` / `REPS` at the top
    of the script.
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
GENE_START = 20134843                      # TAIR10; verified against the published panel e labels
# Inference is run on TAIR10 for every model generation, as the toolkit is used.
# The v3 models were trained on Helixer gene models, which place this gene's start
# 250 bp downstream; scored on that window the v3 model predicts RAP2.12 high
# (0.8311 leaf, 0.9339 root) instead of low. See
# helixerv3_20261109_variants.tsv and rap212_helixer_v3.0.py for that check.
REPS = ["Col-0", "Sichuan", "BuckhornPass", "Utrecht", "Lan-0",
        "Tibet", "Tanz-1", "IP-Rel-0", "IP-Cat-0"]

# Predicted transcription-factor binding regions over the RAP2.12 gTUR, as window
# position ranges. PlantPAN 4.0 output is in none of the archived files, so these
# are transcribed from the published figure and are the only part of this panel
# that is not recomputed here.
TFBS = [((1168, 1172), "PIF5"),
        ((1228, 1241), "ABF4, DNAJ1, MYC3"),
        ((1251, 1257), "ELF3"),
        ((1298, 1302), "EIN6, WRKY"),
        ((1317, 1321), "DNAJ2"),
        ((1336, 1345), "SOC1, AP1"),
        ((1388, 1395), "ERF115")]


def star(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s."


def disp(m):
    """Model name for the plot, with the deepCRE version spelled out."""
    if m.startswith("r3 "):
        return m[3:] + ", deepCRE v3"
    return m.replace("2024 ", "") + ", deepCRE 2024"


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
    """Two plain boundary lines; the tinted region bands were distracting."""
    for x_, lab_ in ((1000.5, "TSS"),):
        ax.axvline(x_, color="black", lw=0.6, ls=":", zorder=1)
        ax.text(x_, 0.965, lab_, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=6.5, color="#4b5563")


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
    published_nine = list(icat[icat.CONTRIBUTOR == "yes"].POS.astype(int))
    best10 = pd.read_csv(os.path.join(scratch, "best10_20261109_positions.tsv"), sep="\t")
    best10_pos = list(best10.genomic_position.astype(int))
    b10 = pd.read_csv(os.path.join(scratch, "best10_20261109_predictions.tsv"), sep="\t")
    sal.update(dict(np.load(os.path.join(scratch, "best10_20261109_saliency.npz"),
                            allow_pickle=True)))

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
        # the deepCRE v3 figure uses the ten SNPs the exhaustive search found optimal
        # for that model; the 2024 figure keeps the published nine
        if gen == "round-3":
            nine, key_lab = best10_pos, "ten SNPs, deepCRE v3 optimum"
            key_sal, key_src = f"{m}||RAP2.12 best ten (v3)", b10
            key_var, ctrl_pre = "best ten (v3)", "control"
        else:
            nine, key_lab = published_nine, "nine key SNPs"
            key_sal, key_src = f"{m}||RAP2.12 nine contributors", f6
            key_var, ctrl_pre = "nine", "control"

        fig = plt.figure(figsize=(13.5, 14.2))
        gs = fig.add_gridspec(5, 12, height_ratios=[1.0, 1.0, 1.05, 1.55, 0.85],
                              hspace=0.62, wspace=1.3)

        # ── a  the 672 haplotypes ───────────────────────────────────────────
        sub = gs[0, 0:4].subgridspec(1, 2, width_ratios=[2.2, 1.0], wspace=0.10)
        d = per[per.model == m].set_index("accession")
        uniq = {}
        for a, sq in geno.items():
            uniq.setdefault(sq, []).append(a)
        hap = []
        for sq, accs in uniq.items():
            grp = {str(d.loc[a, "group"]) for a in accs}
            hap.append(dict(n=len(accs), P=float(d.loc[accs[0], "P"]),
                            group="Hap-W box" if grp == {"Hap-W box"} else "Hap-WT box",
                            members=accs))
        hap = pd.DataFrame(hap)
        axes_a, rng = [], np.random.default_rng(5)
        for j, (grp, c) in enumerate((("Hap-WT box", SITE_HIGH), ("Hap-W box", SITE_LOW))):
            axa = fig.add_subplot(sub[0, j], sharey=axes_a[0] if axes_a else None)
            axes_a.append(axa)
            h = hap[hap.group == grp].sort_values("P")
            xj = rng.uniform(-0.30, 0.30, len(h))
            axa.scatter(xj, h.P, s=18 + 26 * np.log10(h.n + 1), color=c,
                        edgecolor="black", linewidths=0.4, zorder=3)
            for acc in ("Tibet", "IP-Cat-0"):
                sel = h[[acc in mm for mm in h.members]]
                if len(sel):
                    k = list(h.index).index(sel.index[0])
                    axa.scatter([xj[k]], [sel.P.iloc[0]], s=150, facecolor="none",
                                edgecolor=SITE_ACCENT, linewidths=1.2, zorder=4)
                    axa.annotate(acc, (xj[k], sel.P.iloc[0]), textcoords="offset points",
                                 xytext=(0, 13), ha="center", fontsize=7,
                                 color=SITE_ACCENT, fontweight="bold")
            axa.set_xlim(-0.75, 0.75)
            axa.set_xticks([])
            axa.set_xlabel(f"{grp}\n{len(h)} sequences, {int(h.n.sum())} accessions",
                           fontsize=7)
            for sp_ in ("top", "right"):
                axa.spines[sp_].set_visible(False)
            if j:
                axa.spines["left"].set_visible(False)
                axa.tick_params(left=False, labelleft=False)
        axes_a[0].set_ylabel("probability of high gene expression")
        axes_a[0].set_title(f"a)  distinct RAP2.12 promoter sequences\n{disp(m)}",
                            loc="left", fontweight="bold", fontsize=8)
        pw = per[(per.model == m) & (per.group == "Hap-W box")].P.to_numpy()
        pwt = per[(per.model == m) & (per.group == "Hap-WT box")].P.to_numpy()
        t, p = ttest_ind(pwt, pw, equal_var=False)
        u, _ = mannwhitneyu(pw, pwt, alternative="two-sided")
        stats["a"] = [("Hap-WT sequences", int((hap.group == "Hap-WT box").sum()),
                       float(hap[hap.group == "Hap-WT box"].P.mean()),
                       float(hap[hap.group == "Hap-WT box"].P.min()),
                       float(hap[hap.group == "Hap-WT box"].P.max()), np.nan, np.nan),
                      ("Hap-W sequences", int((hap.group == "Hap-W box").sum()),
                       float(hap[hap.group == "Hap-W box"].P.mean()),
                       float(hap[hap.group == "Hap-W box"].P.min()),
                       float(hap[hap.group == "Hap-W box"].P.max()), np.nan, np.nan),
                      ("accession-level Welch t-test", pw.size, float(pw.mean()),
                       np.nan, np.nan, float(t), float(p)),
                      ("AUC", 0, 1 - u / (pw.size * pwt.size), np.nan, np.nan,
                       np.nan, np.nan)]

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
                ("d", 1, key_sal, nine, key_lab)):
            ax = fig.add_subplot(gs[row, 4:12])
            C_LO, C_HI = 950, 1560          # crop: gUR/gTUR boundary plus the variant block
            kref = f"{m}||RAP2.12 reference"
            # the web tool's "Average Saliency Map": the mean over the four bases.
            # Only the observed base is non-zero, so this is the summed score / 4,
            # which is the scale the published panels are drawn on.
            a = sal[kref].mean(axis=1)
            b = sal[key].mean(axis=1)
            xs = np.arange(C_LO, C_HI + 1)
            region_bands(ax)
            for pos in marks:
                if C_LO <= win(pos) <= C_HI:
                    ax.axvline(win(pos), color="#444444", lw=0.5, ls=(0, (2, 2)), zorder=1)
            ax.plot(xs, a[C_LO - 1:C_HI], color=SITE_ORIG, lw=0.6, label=GENE, zorder=4)
            ax.plot(xs, b[C_LO - 1:C_HI], color=SITE_MUT, lw=0.6,
                    label=f"{GENE}: mutated", zorder=5)
            ax.axhline(0, color="black", lw=0.6, ls=(0, (1, 2)))
            ax.set_xlim(C_LO, C_HI)
            ax.set_xlabel("nucleotide position")
            ax.set_ylabel("saliency score")
            ax.legend(frameon=False, loc="lower right", fontsize=6.5, ncol=2)
            ax.set_title(f"{lab})  {ttl}, window {C_LO}-{C_HI}",
                         loc="left", fontweight="bold")
            stats[lab] = [(ttl, len(marks), float(np.abs(b - a).sum()),
                           np.nan, np.nan, np.nan, np.nan)]

        # ── e  the gTUR zoom: the variants sit on the saliency tips ─────────
        ax = fig.add_subplot(gs[2, :])
        a_ = sal[f"{m}||RAP2.12 reference"].mean(axis=1)
        b_ = sal[key_sal].mean(axis=1)
        lo_, hi_ = 1140, 1410
        xs = np.arange(lo_, hi_ + 1)
        for pos in nine:
            ax.axvline(win(pos), color="#444444", lw=0.6, ls=(0, (2, 2)), zorder=1)
        ax.plot(xs, a_[lo_ - 1:hi_], color=SITE_ORIG, lw=0.8, label=GENE, zorder=3)
        ax.plot(xs, b_[lo_ - 1:hi_], color=SITE_MUT, lw=0.8, label=f"{GENE}: mutated",
                zorder=4)
        ax.scatter([win(p_) for p_ in nine], [b_[win(p_) - 1] for p_ in nine], s=22,
                   color=SITE_ACCENT, edgecolor="black", linewidths=0.4, zorder=6)
        for p_ in nine:
            ax.text(win(p_), 1.015, str(win(p_)), transform=ax.get_xaxis_transform(),
                    rotation=90, ha="center", va="bottom", fontsize=5.8, color=SITE_ACCENT)
        ax.axhline(0, color="black", lw=0.6, ls=(0, (1, 2)))
        ax.set_xlim(lo_, hi_)
        ax.set_xlabel("nucleotide position")
        ax.set_ylabel("saliency score")
        ax.legend(frameon=False, loc="upper left", fontsize=6.5, ncol=2)
        ax.set_title(f"e)  gTUR {lo_}-{hi_}, the {len(nine)} SNPs on the saliency profile",
                     loc="left", fontweight="bold", pad=18)
        tips = sum(1 for p_ in nine
                   if abs(b_[win(p_) - 1]) >= np.quantile(np.abs(b_[lo_ - 1:hi_]), 0.75))
        stats["e"] = [(f"{key_lab}: in the top quartile of |saliency| over {lo_}-{hi_}",
                       len(nine), float(tips), np.nan, np.nan, np.nan, np.nan)]

        # ── f  the SNP track of the representative accessions ───────────────
        ax = fig.add_subplot(gs[3, :])
        reps = [a for a in REPS if a in geno]
        nblock = len(geno["Col-0"])
        cols = sorted({p for p in range(GENE_START - 1000 + first,
                                        GENE_START - 1000 + first + nblock)
                       if 0 <= win(p) - first < nblock
                       and any(geno[a][win(p) - first] != geno["Col-0"][win(p) - first]
                               and geno[a][win(p) - first] in "ACGT" for a in reps)})
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
                            + (" †" if p in published_nine and p not in nine else "")
                            for p in cols], rotation=90, fontsize=5.6)
        for tick, p in zip(ax.get_xticklabels(), cols):
            if p in nine:
                tick.set_fontweight("bold")
        # transcription-factor binding regions and which key SNPs fall inside them
        yb = -len(reps) - 0.25
        in_tfbs, tf_of = [], {}
        for (lo_, hi_), name in TFBS:
            idx = [j for j, p_ in enumerate(cols) if lo_ <= win(p_) <= hi_]
            if not idx:
                continue
            ax.add_patch(plt.Rectangle((min(idx) - 0.45, yb - 0.30),
                                       max(idx) - min(idx) + 0.9, 0.60,
                                       facecolor="#e7e9ec", edgecolor="#9ca3af",
                                       linewidth=0.5))
            ax.text((min(idx) + max(idx)) / 2, yb, name, ha="center", va="center",
                    fontsize=5.8)
            for j in idx:
                tf_of[cols[j]] = name
                if cols[j] in nine:
                    in_tfbs.append(cols[j])
        n_hit = len(set(in_tfbs))
        ax.set_xlim(-0.8, len(cols) - 0.2)
        ax.set_ylim(yb - 0.55, 0.7)
        for s_ in ax.spines.values():
            s_.set_visible(False)
        stats["f"] = [(f"key SNP {p_} (window {win(p_)})", 1,
                       1.0 if p_ in tf_of else 0.0, np.nan, np.nan, np.nan, np.nan)
                      for p_ in nine]
        ax.set_title(f"f)  haplotype variants and predicted binding regions — "
                     f"{n_hit} of the {len(nine)} SNPs ({key_lab}) fall inside one"
                     + ("; † = published nine only" if set(nine) != set(published_nine)
                        else ""),
                     loc="left", fontweight="bold", fontsize=8)

        # ── f  the nine contributors against the five controls ──────────────
        # ── g  the nine key SNPs against random substitutions at the same sites ──
        # Supplementary Note 3 states the claim this panel has to carry: five
        # artificial sequences in which the nine key SNPs were substituted with
        # random bases "failed to recreate the observed drop in predicted
        # probabilities". That is one null — the same sites, other alleles — so
        # the panel shows that one null, with 1000 draws instead of five.
        ax = fig.add_subplot(gs[4, 0:6])
        frame = "helixer" if "helixer" in tag else "tair10"
        gen_slug = "2024" if gen == "2024" else "r3"
        npz = os.path.join(scratch,
                           f"perm_20261110_{frame}_{gen_slug}_leaf_published_9.npz")
        z = np.load(npz)
        a_null = z["allele_null"]
        obs = float(z["observed"])
        pv = (int((a_null <= obs).sum()) + 1) / (len(a_null) + 1)
        lo = min(float(a_null.min()), obs) - 0.6
        hi = max(float(a_null.max()), obs) + 0.6
        ax.hist(a_null, bins=np.linspace(lo, hi, 46), color=SITE_ORIG, alpha=0.65,
                edgecolor="white", linewidth=0.4)
        ax.axvline(obs, color=SITE_MUT, lw=1.8, zorder=5)
        right = obs > 0.5 * (lo + hi)
        ax.annotate(f"nine key SNPs, {obs:+.2f} %", xy=(obs, 1.0),
                    xycoords=("data", "axes fraction"),
                    xytext=(-8 if right else 8, -4), textcoords="offset points",
                    va="top", ha="right" if right else "left", fontsize=7.5,
                    color=SITE_MUT, fontweight="bold")
        ax.text(0.985, 0.94, f"{len(a_null)} random substitution sets\n"
                             f"at the same nine sites\n"
                             f"p = {pv:.3f} ({star(pv)})",
                transform=ax.transAxes, ha="right", va="top", fontsize=7)
        ax.set_xlim(lo, hi)
        ax.set_xlabel("change in predicted probability (%)")
        ax.set_ylabel("random substitution sets")
        ax.set_title("g)  the nine key SNPs against random substitutions "
                     "at the same sites", loc="left", fontweight="bold", pad=1)
        stats["g"] = [("observed, published nine", 1, obs, np.nan, np.nan, np.nan, np.nan),
                      ("random alleles, same sites", len(a_null), float(a_null.mean()),
                       float(a_null.min()), float(a_null.max()), np.nan, pv)]
        if "position_null" in z:                 # reported, not plotted
            p_null = z["position_null"]
            pv2 = (int((p_null <= obs).sum()) + 1) / (len(p_null) + 1)
            stats["g"].append(("random sites (not plotted)", len(p_null),
                               float(p_null.mean()), float(p_null.min()),
                               float(p_null.max()), np.nan, pv2))

        ax = fig.add_subplot(gs[4, 6:])
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
                     f"model, as Supplementary Note 3 specifies. Windows are extracted with "
                     f"the TAIR10 annotation for every model generation, as the toolkit is "
                     f"used.\n")
            if tag.startswith("round3"):
                fh.write("Caveat on this generation: the v3 models were trained on Helixer "
                         "gene models, which place the RAP2.12 start 250 bp downstream of the "
                         "TAIR10 one, so only 859 of the 3020 window positions coincide. On "
                         "the TAIR10 window shown here the v3 leaf model predicts RAP2.12 low "
                         "(P = 0.0747) although the gene is highly expressed in both tissues "
                         "(log10 mean TPM 1.71 leaf, 2.25 root; deepPHY label high). Scored "
                         "on a Helixer-anchored window the same model predicts it high "
                         "(0.8311 leaf, 0.9339 root). The low call is therefore a consequence "
                         "of the annotation the window comes from, not of the sequence.\n")
            fh.write("\n")
            fh.write("a) The 672 accessions of Lou et al. (2022) carry only a few distinct "
                     "RAP2.12 promoter sequences; each is drawn once, at its predicted "
                     "probability and at the number of accessions that carry it (log scale). "
                     "Identical sequences necessarily receive identical predictions, so this "
                     "shows what the model does with the sequence variation that exists. Statistics: two-sided Welch's t-test "
                     "between the group means and the Mann-Whitney AUC.\n")
            fh.write(f"b) Predicted probability of the Col-0 reference window against the "
                     f"IP-Cat-0 haplotype (all 17 of its SNPs).\n")
            fh.write("c, d) DeepSHAP saliency averaged over the four bases across the whole "
                     "3020 bp window, reference against the mutated sequence. Dashed lines "
                     "mark the substituted positions; the shaded bands are the regions the web "
                     "tool itself draws (gUR and gDR grey, gTUR red, gTDR blue).\n")
            fh.write("e) The same saliency profiles zoomed to the gTUR window that carries "
                     "the variants, with the nine key SNP positions as dashed lines and the "
                     "value of the mutated profile at each of them marked. The count below is "
                     "how many of the nine sit in the top quartile of |saliency| over that "
                     "window.\n")
            fh.write("f) Bases of the representative accessions at every position of the "
                     "aligned block where at least one of them differs from Col-0, "
                     "read from the Lou et al. alignment. A coloured cell marks a base "
                     "differing from Col-0. Accession labels are coloured by haplotype group. "
                     "The transcription-factor binding sites of the published panel come from "
                     "PlantPAN 4.0 and are in none of the archived files; the boxes drawn "
                     "here are transcribed from the published figure and are the only element "
                     "of this figure that is not recomputed. Panel e reports how many of the "
                     "nine key SNPs fall inside one (1 = inside, 0 = outside, listed below).\n")
            fh.write("g) The control of Supplementary Note 3, scaled up: 1000 instead of "
                     "five artificial sequences in which the nine key SNPs of Lou et al. "
                     "are substituted with a random base that is neither the reference nor "
                     "the IP-Cat-0 allele. Histogram, the 1000 random substitution sets; "
                     "line, the observed change for the nine key SNPs. p is one-sided "
                     "empirical, (#null <= observed + 1)/1001, so its floor is 0.001. A "
                     "second null that also randomises the positions within the aligned "
                     "250 bp block was computed but is not plotted; its p is listed below "
                     "and in perm_20261110_*_summary.tsv.\n\n")
            for lab in ("a", "b", "c", "d", "e", "f", "g"):
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
