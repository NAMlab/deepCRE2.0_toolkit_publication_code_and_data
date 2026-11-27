"""
Purpose:
    Compute MSR predicted-probability heatmap matrices (mean predicted probability per true
    expression class 0/1; rows = MSR model by held-out species, columns = target species) for
    all-species, dicot-only, and monocot-only MSR model sets across 13 species.

Manuscript (Peleke et al., 2026):
    Generates `results/evaluation/msr/suppfig_msr_mean_{leaf,root}_{all,dicots,monocots}_class{0,1}.tsv`
    plotted in Figure 2i/m/n/r/s and Supplementary Figure 1i/m/n/r/s.

Usage & I/O (run from repo root):
    python3 code/03_evaluation/msr_matrices.0.py <tissue> <msr_pred_dir> <out_prefix>
    - Writes `<out_prefix>_<set>_class<k>.tsv` (`set` in `{all, dicots, monocots}`, `k` in `{0, 1}`)

How to change or extend:
    Modify `ORDER` / `SPECIES` / filename regexes if additional species or MSR subsets are added.
"""
import glob
import os
import re
import sys

import pandas as pd

ORDER = ["Atha", "Csat", "Bole", "Cucu", "Gmax", "Dcar", "Slyc", "Cqui", "Bvul", "Osat", "Hvul", "Sbic", "Zmay"]
SP = {"Atha": "arabidopsis_thaliana", "Bvul": "beta_vulgaris", "Bole": "brassica_oleracea", "Csat": "camelina_sativa",
      "Cqui": "chenopodium_quinoa", "Cucu": "cucumis_sativus", "Dcar": "daucus_carota", "Gmax": "glycine_max",
      "Osat": "oryza_sativa", "Slyc": "solanum_lycopersicum", "Sbic": "sorghum_bicolor", "Zmay": "zea_mays",
      "Hvul": "hordeum_vulgare"}
LAB = "data/"
MSR_RE = re.compile(r"^MSR_dPm_(all|dicots|monocots)_(leaf|root)_seed2026_([A-Z][a-z]{3})_dPm_(leaf|root)_msr_train_models_\d{6}_\d{6}\.h5$")


def main(tissue, pred_dir, prefix):
    rows = []
    for code in ORDER:
        f = glob.glob(os.path.join(pred_dir, f"MSR_{tissue}_deepcre_crosspredict_{code}_{tissue}.csv"))
        assert len(f) == 1, (code, f)
        d = pd.read_csv(f[0])
        lab = pd.read_csv(os.path.join(LAB, SP[code], f"{SP[code]}_{tissue}_log10_mean_tpm_Q0.75_labeled_dPHYlogic.csv"),
                          usecols=["gene_id", "target"]).set_index("gene_id").target
        d["y"] = d.genes.map(lab)
        for col in d.columns:
            m = MSR_RE.match(col)
            if not m:
                continue
            for cls in (0, 1):
                rows.append(dict(set=m.group(1), held_out=m.group(3), target=code, cls=cls, mean=d.loc[d.y == cls, col].mean()))
    r = pd.DataFrame(rows)
    for s in ["all", "dicots", "monocots"]:
        for cls in (0, 1):
            m = (r[(r.set == s) & (r.cls == cls)].pivot_table(index="held_out", columns="target", values="mean")
                 .reindex(index=ORDER, columns=ORDER))
            m.to_csv(f"{prefix}_{s}_class{cls}.tsv", sep="\t")
            print(s, cls, int(m.notna().sum().sum()), "cells")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
