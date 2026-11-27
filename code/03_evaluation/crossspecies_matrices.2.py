"""
Purpose:
    Aggregate raw `deepcre_crosspredict.py` prediction CSVs into 13 x 13 matrices of mean
    predicted probability of high expression per model species (rows) x target species (columns),
    split into true low-expression (`class0`) and true high-expression (`class1`) genes.

Manuscript (Peleke et al., 2026):
    Produces `results/evaluation/final/crossspecies_mean_{leaf,root}_class{0,1}.tsv` plotted
    in Figure 2d–e and Supplementary Figure 1d–e (13 species including *H. vulgare*; mean
    probability per class; grass order `Osat, Hvul, Sbic, Zmay`).

Usage & I/O (run from repo root):
    python3 code/03_evaluation/crossspecies_matrices.2.py \
        <out_prefix> <label_dir> <tissue> <models.tsv> <pred_dir> [<pred_dir> ...]
    - `<models.tsv>`: `code, file` (e.g. `results/evaluation/final/deployed_models_leaf.tsv`)
    - `<pred_dir>`: one or more directories holding `EXTERNAL/predictions/` CSVs
    - Writes `<out_prefix>_class0.tsv` and `<out_prefix>_class1.tsv`

How to change or extend:
    Extend `ORDER` and `SPECIES` dictionaries to include additional model/target species, or
    pass an alternative `<models.tsv>` to compare different fold selections.
"""
import glob
import os
import sys

import pandas as pd

ORDER = ["Atha", "Csat", "Bole", "Cucu", "Gmax", "Dcar", "Slyc", "Cqui", "Bvul", "Osat", "Hvul", "Sbic", "Zmay"]
SP = {"Atha": "arabidopsis_thaliana", "Bvul": "beta_vulgaris", "Bole": "brassica_oleracea", "Csat": "camelina_sativa",
      "Cqui": "chenopodium_quinoa", "Cucu": "cucumis_sativus", "Dcar": "daucus_carota", "Gmax": "glycine_max",
      "Osat": "oryza_sativa", "Slyc": "solanum_lycopersicum", "Sbic": "sorghum_bicolor", "Zmay": "zea_mays",
      "Hvul": "hordeum_vulgare"}


def main(prefix, label_dir, tissue, models_tsv, pred_dirs, label_pattern="{sp}_{tissue}_log10_mean_tpm_Q0.75_labeled_dPHYlogic.csv"):
    models = pd.read_csv(models_tsv, sep="\t")
    want = {os.path.basename(f) for f in models.file}
    code_of = {os.path.basename(f): c for f, c in zip(models.file, models.code)}
    rows = []
    for pdir in pred_dirs:
        for f in glob.glob(os.path.join(pdir, "*.csv")):
            target = f.split("crosspredict_")[1][:4]
            target = "Dcar" if target == "Dcau" else target
            if target not in SP:
                continue
            cols = [c for c in pd.read_csv(f, nrows=0).columns if c in want]
            if not cols:
                continue
            d = pd.read_csv(f, usecols=cols + ["genes"])
            sp = SP[target]
            lab = pd.read_csv(os.path.join(label_dir, sp, label_pattern.format(sp=sp, tissue=tissue)), usecols=["gene_id", "target"])
            d["y"] = d.genes.map(lab.set_index("gene_id").target)
            for cls in (0, 1):
                med = d[d.y == cls][cols].mean()
                for c, v in med.items():
                    rows.append(dict(target=target, model_species=code_of[c], model=c, cls=cls, mean=v))
    r = pd.DataFrame(rows).drop_duplicates(["target", "model", "cls"])
    for cls in (0, 1):
        m = (r[r.cls == cls].groupby(["model_species", "target"])["mean"].mean().unstack("target")
             .reindex(index=ORDER, columns=ORDER))
        m.to_csv(f"{prefix}_class{cls}.tsv", sep="\t")
        print(f"class {cls}: {m.notna().sum().sum()} of {m.size} cells")


if __name__ == "__main__":
    if len(sys.argv) < 6:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5:])
