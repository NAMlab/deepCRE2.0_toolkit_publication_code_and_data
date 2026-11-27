"""
Purpose:
    Evaluate all 52 multi-species (MSR: `all`, `dicots`, `monocots`) models on every target
    species (`held_out`, `seen`, `outside`) using 100 per-(model, target)-seeded class-balanced
    draws (threshold 0.5) and compare against deployed single-species (SSR) models.

Manuscript (Peleke et al., 2026):
    Generates `results/evaluation/msr/msr_eval_{leaf,root}_{full,summary}.tsv` plotted in
    Figure 2f–s, Supplementary Figure 1f–s, and included in Supplementary Data 1.

Usage & I/O (run from repo root):
    python3 code/03_evaluation/msr_eval_dPm.2.py \
        <tissue> <msr_pred_dir> <ssr_pred_dir> <S2v8_per_species.tsv> <out_prefix>
    - `<msr_pred_dir>`: MSR cross-prediction CSVs (`EXTERNAL/predictions/`)
    - `<ssr_pred_dir>`: deployed SSR cross-prediction CSVs (merged across folders if split)
    - Writes `<out_prefix>_full.tsv` and `<out_prefix>_summary.tsv`

How to change or extend:
    Update `ORDER`, `DICOTS`, `MONOCOTS`, or `SPECIES` mappings for new species/clade groups,
    or adjust `N_DRAWS` (default 100) at the top of the module.
"""
import glob
import os
import re
import sys
import zlib

import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, roc_auc_score

ORDER = ["Atha", "Csat", "Bole", "Cucu", "Gmax", "Dcar", "Slyc", "Cqui", "Bvul", "Osat", "Hvul", "Sbic", "Zmay"]
SP = {"Atha": "arabidopsis_thaliana", "Bvul": "beta_vulgaris", "Bole": "brassica_oleracea", "Csat": "camelina_sativa",
      "Cqui": "chenopodium_quinoa", "Cucu": "cucumis_sativus", "Dcar": "daucus_carota", "Gmax": "glycine_max",
      "Osat": "oryza_sativa", "Slyc": "solanum_lycopersicum", "Sbic": "sorghum_bicolor", "Zmay": "zea_mays",
      "Hvul": "hordeum_vulgare"}
SETS = {"all": set(ORDER), "dicots": {"Atha", "Bvul", "Bole", "Csat", "Cqui", "Cucu", "Dcar", "Gmax", "Slyc"},
        "monocots": {"Osat", "Sbic", "Zmay", "Hvul"}}
LAB = "data/"
MSR_RE = re.compile(r"^MSR_dPm_(all|dicots|monocots)_(leaf|root)_seed2026_([A-Z][a-z]{3})_dPm_(leaf|root)_msr_train_models_\d{6}_\d{6}\.h5$")
N_DRAWS = 100


def balanced(y, p, seed):
    rng = np.random.default_rng(seed)
    lo, hi = np.where(y == 0)[0], np.where(y == 1)[0]
    k = min(len(lo), len(hi))
    res = []
    for _ in range(N_DRAWS):
        idx = np.concatenate([rng.choice(lo, k, replace=False), rng.choice(hi, k, replace=False)])
        yy, pp = y[idx], p[idx]
        c = (pp > 0.5).astype(int)
        res.append(((c == yy).mean(), roc_auc_score(yy, pp), matthews_corrcoef(yy, c)))
    a = np.array(res).mean(0)
    return dict(n_balanced=2 * k, accuracy=a[0], auROC=a[1], MCC=a[2])


def score(pred_dir, pattern, tissue, parse):
    rows = []
    for code in ORDER:
        f = sorted(glob.glob(os.path.join(pred_dir, pattern.format(code=code, tissue=tissue))))
        if not f:
            continue
        d = None
        for ff in f:
            x = pd.read_csv(ff)
            x = x[["genes"] + [c for c in x.columns if c.endswith(".h5")]]
            if d is None:
                d = x
            else:
                new = [c for c in x.columns if c.endswith(".h5") and c not in d.columns]
                n0 = len(d)
                d = d.merge(x[["genes"] + new], on="genes", how="inner")
                assert len(d) == n0, (code, ff, n0, len(d))   # same genes extracted in every file
        lab = pd.read_csv(os.path.join(LAB, SP[code], f"{SP[code]}_{tissue}_log10_mean_tpm_Q0.75_labeled_dPHYlogic.csv"),
                          usecols=["gene_id", "target"]).set_index("gene_id").target
        d["y"] = d.genes.map(lab)
        d = d[d.y.isin([0, 1])]
        y = d.y.to_numpy().astype(int)
        for col in [c for c in d.columns if c.endswith(".h5")]:
            info = parse(col)
            if info is None:
                continue
            r = dict(target=code, model=col, **info)
            r.update(balanced(y, d[col].to_numpy(), zlib.crc32(f"{col}|{code}".encode())))
            rows.append(r)
    return pd.DataFrame(rows)


def main(tissue, msr_dir, ssr_dir, s2_species, prefix):
    def parse_msr(col):
        m = MSR_RE.match(col)
        if not m:
            return None
        return dict(kind="MSR", set=m.group(1), held_out=m.group(3))

    def parse_ssr(col):
        return dict(kind="SSR_deployed", set="single", held_out=col[:4])

    msr = score(msr_dir, "MSR_{tissue}_deepcre_crosspredict_{code}_{tissue}.csv", tissue, parse_msr)
    msr["status"] = np.where(msr.target == msr.held_out, "held_out",
                             np.where([t in SETS[s] for t, s in zip(msr.target, msr.set)], "seen", "outside"))
    ssr = score(ssr_dir, "deployed_{tissue}_*_deepcre_crosspredict_{code}_{tissue}.csv", tissue, parse_ssr)
    ssr["status"] = np.where(ssr.target == ssr.held_out, "own_species", "other_species")
    full = pd.concat([msr, ssr], ignore_index=True)
    full["tissue"] = tissue
    full.to_csv(prefix + "_full.tsv", sep="\t", index=False)
    s2 = pd.read_csv(s2_species, sep="\t")
    s2 = s2[s2.tissue == tissue]
    summ = []
    for code in ORDER:
        r = dict(tissue=tissue, target=code)
        own = s2[s2.species.str.startswith({"Atha": "Arabidopsis", "Csat": "Camelina", "Bole": "Brassica", "Cucu": "Cucumis",
                                             "Gmax": "Glycine", "Dcar": "Daucus", "Slyc": "Solanum", "Cqui": "Chenopodium",
                                             "Bvul": "Beta", "Osat": "Oryza", "Sbic": "Sorghum", "Zmay": "Zea", "Hvul": "Hordeum"}[code])]
        r["SSR_own_species_fold_mean_MCC"] = float(own.MCC_v8_mean.iloc[0]) if len(own) else np.nan
        o = ssr[(ssr.target == code) & (ssr.status == "other_species")].MCC
        r["SSR_other_species_median_MCC"], r["SSR_best_other_species_MCC"] = o.median(), o.max()
        for st in ["all", "dicots", "monocots"]:
            h = msr[(msr.target == code) & (msr.set == st) & (msr.status == "held_out")]
            out = msr[(msr.target == code) & (msr.set == st) & (msr.status == "outside")].MCC
            r[f"MSR_{st}_held_out_MCC"] = float(h.MCC.iloc[0]) if len(h) else np.nan
            r[f"MSR_{st}_held_out_accuracy"] = float(h.accuracy.iloc[0]) if len(h) else np.nan
            r[f"MSR_{st}_held_out_auROC"] = float(h.auROC.iloc[0]) if len(h) else np.nan
            r[f"MSR_{st}_outside_median_MCC"] = out.median() if len(out) else np.nan
        summ.append(r)
    summ = pd.DataFrame(summ)
    summ.round(4).to_csv(prefix + "_summary.tsv", sep="\t", index=False)
    cols = ["target", "SSR_own_species_fold_mean_MCC", "SSR_other_species_median_MCC", "SSR_best_other_species_MCC",
            "MSR_all_held_out_MCC", "MSR_dicots_held_out_MCC", "MSR_monocots_held_out_MCC", "MSR_dicots_outside_median_MCC",
            "MSR_monocots_outside_median_MCC"]
    print(summ[cols].round(3).to_string(index=False))
    print("medians:", summ[cols[1:]].median().round(3).to_dict())


if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    main(*sys.argv[1:])
