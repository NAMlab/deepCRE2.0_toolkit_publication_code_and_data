"""
Purpose:
    Combine `S2v8_per_fold.tsv` metrics with exact training-set and validation-set gene counts
    (low, high, balanced = `2 * min`) across all converged chromosome folds.

Manuscript (Peleke et al., 2026):
    Generates `results/evaluation/final/SuppTable1_per_fold_counts.tsv` backing Supplementary Table 1
    and the `fold_sizes` sheet of Supplementary Data 1.

Usage & I/O (run from repo root):
    python3 code/03_evaluation/supptable1_fold_counts.0.py \
        <S2v8_per_fold.tsv> <label_dir> <out.tsv> <pred_dir> [<pred_dir> ...]

How to change or extend:
    Update `ORDER` / `SPECIES` mappings to support additional species.
"""
import glob
import os
import sys

import pandas as pd

SP = {"Atha": "arabidopsis_thaliana", "Bvul": "beta_vulgaris", "Bole": "brassica_oleracea", "Csat": "camelina_sativa",
      "Cqui": "chenopodium_quinoa", "Cucu": "cucumis_sativus", "Dcar": "daucus_carota", "Gmax": "glycine_max",
      "Osat": "oryza_sativa", "Slyc": "solanum_lycopersicum", "Sbic": "sorghum_bicolor", "Zmay": "zea_mays",
      "Hvul": "hordeum_vulgare"}
ORDER = ["Atha", "Csat", "Bole", "Cucu", "Gmax", "Dcar", "Slyc", "Cqui", "Bvul", "Osat", "Hvul", "Sbic", "Zmay"]
CHROM_RE = r"_(N[CW]_\d+\.\d+)_\d+$"


def main(per_fold, label_dir, out, pred_dirs):
    pf = pd.read_csv(per_fold, sep="\t")
    pf = pf[pf.status.isin(["converged", "converged_rescue"])]
    genes = {}
    for pdir in pred_dirs:
        for f in sorted(glob.glob(os.path.join(pdir, "*.csv"))):
            code = f.split("crosspredict_")[1][:4]
            code = "Dcar" if code == "Dcau" else code
            tis = os.path.basename(f)[:-4].split("_")[-1]
            if code in SP and (code, tis) not in genes:
                genes[(code, tis)] = pd.read_csv(f, usecols=["genes"]).genes
    rows = []
    for (tissue, code), g in pf.groupby(["tissue", "code"]):
        sp = SP[code]
        lab = pd.read_csv(os.path.join(label_dir, sp, f"{sp}_{tissue}_log10_mean_tpm_Q0.75_labeled_dPHYlogic.csv"),
                          usecols=["gene_id", "target"]).set_index("gene_id").target
        d = pd.DataFrame({"gene": genes[(code, tissue)]})
        d["chrom"] = d.gene.str.extract(CHROM_RE)[0]
        d["y"] = d.gene.map(lab)
        d = d[d.y.isin([0, 1])]
        folds = set(g.chrom)
        for _, r in g.iterrows():
            t = d[d.chrom.isin(folds - {r.chrom})]
            tlo, thi = int((t.y == 0).sum()), int((t.y == 1).sum())
            rows.append(dict(tissue=tissue, species=r.species, code=code, fold_chromosome=r.chrom,
                             train_low=tlo, train_high=thi, train_balanced=2 * min(tlo, thi),
                             val_low=int(r.n_low_v8), val_high=int(r.n_high_v8), val_balanced=int(r.n_v8_balanced),
                             deployed=bool(r.deployed), accuracy=round(r.accuracy_v8, 4), auROC=round(r.auROC_v8, 4),
                             MCC=round(r.MCC_v8, 4)))
    o = pd.DataFrame(rows)
    o["_s"] = o.code.map(ORDER.index)
    o = o.sort_values(["tissue", "_s", "fold_chromosome"]).drop(columns="_s")
    o.to_csv(out, sep="\t", index=False)
    print(o.groupby(["tissue", "code"]).size().unstack(0).to_string())


if __name__ == "__main__":
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])
