"""
Purpose:
    Compute per-fold, per-species, and per-tissue validation metrics for all final SSR models:
    (a) 100 CRC32-seeded class-balanced draws on `v8` paralog-filtered held-out-chromosome genes
    (`accuracy_v8`, `auROC_v8`, `MCC_v8`) and (b) unbalanced metrics across all held-out
    chromosome low/high genes; selects the deployed fold per species x tissue by max `auROC_v8`.

Manuscript (Peleke et al., 2026):
    Produces `S2v8_per_fold.tsv`, `S2v8_per_species.tsv`, `S2v8_per_tissue.tsv`, and
    `deployed_models_{leaf,root}.tsv` under `results/evaluation/final/`.

Usage & I/O (run from repo root):
    python3 code/03_evaluation/s2v8_validation_metrics.2.py \
        <final_models.tsv> <label_dir> <out_dir> <pred_dir> <v8.pickle> [<v8.pickle> ...]

How to change or extend:
    Adjust `N_DRAWS`, `ORDER`, or fold selection criteria (`idxmax` on `auROC_v8`), or pass
    additional validation pickles cleanly as positional trailing arguments.
"""
import glob
import zlib
import os
import pickle
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score, matthews_corrcoef, roc_auc_score

N_DRAWS = 100
SP = {"Atha": "arabidopsis_thaliana", "Bvul": "beta_vulgaris", "Bole": "brassica_oleracea", "Csat": "camelina_sativa",
      "Cqui": "chenopodium_quinoa", "Cucu": "cucumis_sativus", "Dcar": "daucus_carota", "Gmax": "glycine_max",
      "Osat": "oryza_sativa", "Slyc": "solanum_lycopersicum", "Sbic": "sorghum_bicolor", "Zmay": "zea_mays",
      "Hvul": "hordeum_vulgare"}
CHROM_RE = r"_(N[CW]_\d+\.\d+)_\d+$"


def basic(y, p):
    c = (p > 0.5).astype(int)
    return dict(accuracy=float((c == y).mean()), auROC=float(roc_auc_score(y, p)), MCC=float(matthews_corrcoef(y, c)),
                balanced_accuracy=float(balanced_accuracy_score(y, c)))


def main(final_tsv, label_dir, out_dir, pred_dir, pickles):
    os.makedirs(out_dir, exist_ok=True)
    vgenes = {}
    for pk in pickles:
        vgenes.update({k: set(v) for k, v in pickle.load(open(pk, "rb")).items()})
    fm = pd.read_csv(final_tsv, sep="\t")
    rng = np.random.default_rng(2026)
    rows = []
    for _, r in fm[fm.status != "missing"].iterrows():
        base = os.path.basename(r.file)[:-3]
        hits = glob.glob(os.path.join(pred_dir, f"{base}_deepcre_crosspredict_{r.code}_S2v8*.csv"))
        out = dict(tissue=r.tissue, code=r.code, species=SP[r.code].replace("_", " ").capitalize(), chrom=r.chrom,
                   status=r.status, model=os.path.basename(r.file))
        if not hits:
            rows.append(dict(out, note="no prediction file")); continue
        d = pd.read_csv(sorted(hits)[-1])
        col = [c for c in d.columns if c.endswith(".h5")][0]
        sp = SP[r.code]
        lab = pd.read_csv(os.path.join(label_dir, sp, f"{sp}_{r.tissue}_log10_mean_tpm_Q0.75_labeled_dPHYlogic.csv"),
                          usecols=["gene_id", "target"]).set_index("gene_id").target
        d["chrom"] = d.genes.str.extract(CHROM_RE)[0]
        d["y"] = d.genes.map(lab)
        held = d[(d.chrom == r.chrom) & d.y.isin([0, 1])]
        y, p = held.y.to_numpy().astype(int), held[col].to_numpy()
        out.update(n_low_all=int((y == 0).sum()), n_high_all=int((y == 1).sum()))
        out.update({k + "_all": v for k, v in basic(y, p).items()})
        v = held[held.genes.isin(vgenes[r.code.lower()])]
        yv, pv = v.y.to_numpy().astype(int), v[col].to_numpy()
        lo, hi = np.where(yv == 0)[0], np.where(yv == 1)[0]
        k = min(len(lo), len(hi))
        out.update(n_low_v8=len(lo), n_high_v8=len(hi), n_v8_balanced=2 * k)
        if k > 0:
            rng = np.random.default_rng(zlib.crc32(os.path.basename(r.file).encode()))
            res = []
            for _ in range(N_DRAWS):
                idx = np.concatenate([rng.choice(lo, k, replace=False), rng.choice(hi, k, replace=False)])
                res.append(basic(yv[idx], pv[idx]))
            res = pd.DataFrame(res)
            for c in ("accuracy", "auROC", "MCC"):
                out[c + "_v8"] = res[c].mean()
                out[c + "_v8_sd_draws"] = res[c].std()
        out["pred_sd"] = float(np.std(p))
        rows.append(out)
    f = pd.DataFrame(rows)
    ok = f[f.status.isin(["converged", "converged_rescue"]) & f.auROC_v8.notna()].copy()
    ok["deployed"] = False
    ok.loc[ok.groupby(["tissue", "code"]).auROC_v8.idxmax(), "deployed"] = True
    f = f.merge(ok[["model", "deployed"]], on="model", how="left")
    f["deployed"] = f.deployed.fillna(False)
    f.to_csv(os.path.join(out_dir, "S2v8_per_fold.tsv"), sep="\t", index=False)

    mets = ["accuracy_v8", "auROC_v8", "MCC_v8", "MCC_all", "auROC_all", "balanced_accuracy_all"]
    g = ok.groupby(["tissue", "species"])
    s = g[mets].agg(["mean", "std"]).round(3)
    s.columns = [f"{a}_{b}" for a, b in s.columns]
    s["n_folds"] = g.size()
    s["n_folds_not_converged"] = f[f.status == "not_converged"].groupby(["tissue", "species"]).size()
    s["median_n_v8_balanced"] = g.n_v8_balanced.median()
    s["median_n_all"] = g.apply(lambda x: (x.n_low_all + x.n_high_all).median())
    dep = ok[ok.deployed].set_index(["tissue", "species"])
    for c in ["chrom", "accuracy_v8", "auROC_v8", "MCC_v8", "MCC_all", "auROC_all"]:
        s["deployed_" + c] = dep[c].round(3) if c != "chrom" else dep[c]
    s = s.reset_index()
    s.to_csv(os.path.join(out_dir, "S2v8_per_species.tsv"), sep="\t", index=False)
    t = s.groupby("tissue")[[m + "_mean" for m in mets]].agg(["mean", "std", "min", "max"]).round(3)
    t.to_csv(os.path.join(out_dir, "S2v8_per_tissue.tsv"), sep="\t")
    print(f.status.value_counts().to_string())
    print(t.T.to_string())


if __name__ == "__main__":
    if len(sys.argv) < 6:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5:])
