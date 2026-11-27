"""
Purpose:
    Export the multi-sheet Source Data Excel workbook for Figure 2 (leaf), Supplementary Figure 1
    (root), and Supplementary Figure 2 (cross-tissue evaluation), directly importing UPGMA / Mantel
    distance calculations from `make_pub_figures_split.7.py` so Excel values match plotted figures.

Manuscript (Peleke et al., 2026):
    Produces Source Data Excel files for Figure 2 (`a–s`, excluding tree-only cladograms `a,f,j,o`),
    Supplementary Figure 1 (via `tissue` column), and Supplementary Figure 2 (`b–o`).

Usage & I/O (run from repo root):
    python3 code/04_figures/build_source_data_fig2.5.py \
        <S2_dir> <fig_script_dir> <suppfig2_values.tsv> <out.xlsx>
    Example:
    python3 code/04_figures/build_source_data_fig2.5.py \
        results/evaluation code/04_figures <out/SuppFig_tissue_v7_values.tsv> out.xlsx

How to change or extend:
    Update sheet mapping or panel letters inside `build()` if figure panels are reordered;
    inherits phylogeny order `ORDER` automatically from `make_pub_figures_split.7.py`.
"""
import csv
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

SP = {"Atha": "Arabidopsis thaliana", "Bvul": "Beta vulgaris", "Bole": "Brassica oleracea",
      "Csat": "Camelina sativa", "Cqui": "Chenopodium quinoa", "Cucu": "Cucumis sativus",
      "Dcar": "Daucus carota", "Gmax": "Glycine max", "Hvul": "Hordeum vulgare",
      "Osat": "Oryza sativa", "Slyc": "Solanum lycopersicum", "Sbic": "Sorghum bicolor",
      "Zmay": "Zea mays"}
ORDER = ["Atha", "Csat", "Bole", "Cucu", "Gmax", "Dcar", "Slyc", "Cqui", "Bvul", "Osat", "Hvul", "Sbic", "Zmay"]
SETS = {"all": ("all 13 species", "i"), "dicots": ("eudicots only", "m / n"), "monocots": ("grasses only", "r / s")}
DOT_PANEL = {"all": ("g", "h"), "dicots": ("k", "l"), "monocots": ("p", "q")}

LEGEND = (
    "Source Data for Figure 2. deepCRE model performance for single-species reference (SSR) and multi-species "
    "reference (MSR) training across 13 flowering plant species, leaf and root. Panels a, f, j and o are cladograms and "
    "carry no data. Sheets: Fig2b_c_folds (every chromosome fold behind the grey dots of b and c), Fig2b_c_summary "
    "(the web-tool model value, fold mean and s.d. drawn in b and c; n_folds = the n printed in a), Fig2d_SSR_crosspred (the two heatmaps of d), "
    "Fig2e_UPGMA (the Euclidean distances between the 52-value cross-prediction profiles clustered in e, with the "
    "resulting leaf order and the Mantel statistics), Fig2_MSR_dots (the black dots, open diamonds and grey dots of "
    "g, h, k, l, p and q), Fig2_MSR_crosspred (the MSR heatmaps i, m, n, r and s). Labels: deepPHY logic on log10 "
    "mean TPM; SSR validation genes: paralog-filtered list v8; seed 2026. SSR fold metrics are computed on the "
    "paralog-filtered validation genes of the held-out chromosome; all MSR and all cross-species values are computed "
    "on all low- and high-expression genes of the target species. Both are balanced by random downsampling of the "
    "larger class (mean of 100 draws, threshold 0.5)."
)


def load_fig_module(fig_dir):
    path = os.path.join(fig_dir, "make_pub_figures_split.7.py")
    spec = importlib.util.spec_from_file_location("figv7", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["figv7"] = mod
    spec.loader.exec_module(mod)
    return mod


def read_matrix(path, value_name, **fixed):
    rows = []
    with open(path, encoding="utf-8") as f:
        r = csv.reader(f, delimiter="\t")
        cols = next(r)[1:]
        for line in r:
            if not line or not line[0]:
                continue
            for c, v in zip(cols, line[1:]):
                if v not in ("", "nan", "NA"):
                    rows.append({**fixed, "model_species": SP[line[0]], "target_species": SP[c],
                                 value_name: round(float(v), 4)})
    return rows


def suppfig2_sheets(values_tsv):
    v = pd.read_csv(values_tsv, sep="\t")
    v["species"] = v.species.map(SP)
    pick = lambda ps: v[v.panel.isin(ps)].dropna(axis=1, how="all")
    return [("SuppFig2b_e_SSR_models", pick(["b", "c", "d", "e"])),
            ("SuppFig2f_i_MSR_models", pick(["f", "g", "h", "i"])),
            ("SuppFig2j_class_composition", pick(["j"])),
            ("SuppFig2k_flipped_gene_counts", pick(["k"])),
            ("SuppFig2l_n_prediction_contrast", pick(["l", "n"])),
            ("SuppFig2m_o_flipped_performance", pick(["m", "o"]))]


def build(s2, fig_dir, sf2_values, out):
    F = load_fig_module(fig_dir)

    pf = pd.read_csv(os.path.join(s2, "final", "S2v8_per_fold.tsv"), sep="\t")
    pf = pf[pf.status.isin(["converged", "converged_rescue"])].copy()
    pf["_s"] = pf.code.map(ORDER.index)
    pf = pf.sort_values(["tissue", "_s", "chrom"])
    folds = pd.DataFrame({"panel": "b, c", "tissue": pf.tissue, "species": pf.code.map(SP),
                          "fold_validation_chromosome": pf.chrom, "model_id": pf.model,
                          "auROC": pf.auROC_v8.round(4), "accuracy": pf.accuracy_v8.round(4),
                          "MCC": pf.MCC_v8.round(4), "n_validation_balanced": pf.n_v8_balanced.astype(int),
                          "deployed_on_web_tool": pf.deployed.astype(bool)})

    summ = []
    for t in ("leaf", "root"):
        for code in ORDER:
            g = pf[(pf.tissue == t) & (pf.code == code)]
            dep = g[g.deployed]
            assert len(dep) == 1, (t, code)
            for panel, col, name in (("b", "auROC_v8", "auROC"), ("c", "accuracy_v8", "accuracy")):
                v = g[col].to_numpy(float)
                summ.append(dict(panel=panel, tissue=t, species=SP[code], metric=name,
                                 web_tool_model_value=round(float(dep[col].iloc[0]), 4),
                                 fold_mean=round(float(v.mean()), 4),
                                 fold_sd=round(float(v.std(ddof=1)), 4) if len(v) > 1 else 0.0,
                                 fold_min=round(float(v.min()), 4), fold_max=round(float(v.max()), 4),
                                 n_folds=len(v)))
    summary = pd.DataFrame(summ)

    d_rows, e_rows, mantel_rows = [], [], []
    for t in ("leaf", "root"):
        mats = {}
        for cls, lab in ((0, "lowly expressed"), (1, "highly expressed")):
            p = os.path.join(s2, "final", f"crossspecies_mean_{t}_class{cls}.tsv")
            mats[cls] = p
            d_rows += read_matrix(p, "mean_predicted_probability", panel="d", tissue=t, gene_class=lab)
        m0, _ = F.read_matrix_tsv(F.Path(mats[0]))
        m1, _ = F.read_matrix_tsv(F.Path(mats[1]))
        a0 = F.reindex_matrix(m0, ORDER, ORDER)
        a1 = F.reindex_matrix(m1, ORDER, ORDER)
        Xp = np.hstack([a0, a1, a0.T, a1.T])
        n = len(ORDER)
        D = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                D[i, j] = float(np.linalg.norm(Xp[i] - Xp[j]))
        for i in range(n):
            for j in range(n):
                e_rows.append(dict(panel="e", tissue=t, species_1=SP[ORDER[i]], species_2=SP[ORDER[j]],
                                   euclidean_distance=round(float(D[i, j]), 4)))
        # the layout order and the Mantel statistics are taken from the function that draws panel e,
        # so the source data cannot disagree with the figure or its legend
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig_tmp, ax_tmp = plt.subplots()
        info = F.draw_perf_tree(ax_tmp, m0, m1)
        plt.close(fig_tmp)
        mantel_rows.append(dict(panel="e", tissue=t, n_species=n, n_species_pairs=n * (n - 1) // 2,
                                profile_length=4 * n, mantel_r=round(float(info["mantel_r"]), 4),
                                mantel_p=round(float(info["mantel_p"]), 4),
                                upgma_leaf_order=", ".join(info["order"])))
    fig2d = pd.DataFrame(d_rows)
    fig2e = pd.DataFrame(e_rows)
    fig2e_stats = pd.DataFrame(mantel_rows)

    dots, heat = [], []
    for t in ("leaf", "root"):
        full = F.read_tsv_dicts(F.Path(os.path.join(s2, "msr", f"msr_eval_{t}_full.tsv")))
        pfr = F.read_tsv_dicts(F.Path(os.path.join(s2, "final", "S2v8_per_fold.tsv")))
        pfr = [x for x in pfr if x["status"] in ("converged", "converged_rescue")]
        sm, oth = F.compute_msr_summary(full, pfr, t)
        for st, (label, _) in SETS.items():
            members = [c for c in ORDER if sm[c].get(f"{st}_kind") == "held_out"]
            pa, pb = DOT_PANEL[st]
            for code in members:
                s = sm[code]
                o_auc = [float(x["auROC"]) for x in oth if x["target"] == code]
                o_acc = [float(x["accuracy"]) for x in oth if x["target"] == code]
                for panel, metric, mval, own, others in (
                        (pa, "auROC", s[f"{st}_auROC"], s["SSR_deployed_own_auROC"], o_auc),
                        (pb, "accuracy", s[f"{st}_accuracy"], s["SSR_deployed_own_accuracy"], o_acc)):
                    dots.append(dict(panel=panel, tissue=t, msr_training_set=label, held_out_species=SP[code],
                                     metric=metric, MSR_on_held_out_species=round(mval, 4),
                                     SSR_web_tool_model_same_species=round(own, 4),
                                     SSR_other_species_median=round(float(np.median(others)), 4),
                                     SSR_other_species_min=round(float(np.min(others)), 4),
                                     SSR_other_species_max=round(float(np.max(others)), 4),
                                     n_SSR_other_species=len(others)))
            for cls, lab in ((0, "lowly expressed"), (1, "highly expressed")):
                heat += read_matrix(os.path.join(s2, "msr", f"suppfig_msr_mean_{t}_{st}_class{cls}.tsv"),
                                    "mean_predicted_probability", panel=SETS[st][1], tissue=t,
                                    msr_training_set=label, gene_class=lab)
    fig2_dots = pd.DataFrame(dots)
    fig2_heat = pd.DataFrame(heat).rename(columns={"model_species": "held_out_species_of_model"})

    sheets = [("Fig2b_c_folds", folds), ("Fig2b_c_summary", summary), ("Fig2d_SSR_crosspred", fig2d),
              ("Fig2e_UPGMA", fig2e), ("Fig2e_statistics", fig2e_stats), ("Fig2_MSR_dots", fig2_dots),
              ("Fig2_MSR_crosspred", fig2_heat)] + suppfig2_sheets(sf2_values)
    with pd.ExcelWriter(out, engine="openpyxl") as w:
        pd.DataFrame([[LEGEND]]).to_excel(w, sheet_name="legend", index=False, header=False)
        for name, df in sheets:
            df.to_excel(w, sheet_name=name, index=False)
            print(f"  {name}: {len(df)} rows x {len(df.columns)} cols")
    print("written:", out)
    for _, r in fig2e_stats.iterrows():
        print(f"  {r.tissue}: Mantel r={r.mantel_r} p={r.mantel_p} | order {r.upgma_leaf_order}")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    build(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
