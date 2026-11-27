"""
Purpose:
    Test tissue specificity by scoring each species' deployed leaf and root SSR models
    exclusively on genes whose expression class flips between leaf and root (`high in leaf /
    low in root` vs. `low in leaf / high in root`, where chance auROC is 0.5).

Manuscript (Peleke et al., 2026):
    Generates `results/evaluation/tissue_cross/tissue_switchers.tsv` for Supplementary Figure 2.

Usage & I/O (run from repo root):
    python3 code/03_evaluation/tissue_cross_switchers.1.py <ssr_pred_dir> <out.tsv>
    - Imports shared prediction/label loaders from `tissue_cross_eval.1.py`

How to change or extend:
    Modify shared gene/label helpers in `tissue_cross_eval.1.py` or adjust output columns here.
"""
import glob
import importlib.util
import os
import sys
import zlib

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("tcross", os.path.join(HERE, "tissue_cross_eval.1.py"))
T = importlib.util.module_from_spec(spec)
spec.loader.exec_module(T)
ORDER = T.ORDER
N_DRAWS = 100


def balanced_auroc(y, p, seed):
    rng = np.random.default_rng(seed)
    lo, hi = np.where(y == 0)[0], np.where(y == 1)[0]
    k = min(len(lo), len(hi))
    if k < 10:
        return np.nan, 2 * k
    v = []
    for _ in range(N_DRAWS):
        idx = np.concatenate([rng.choice(lo, k, replace=False), rng.choice(hi, k, replace=False)])
        v.append(roc_auc_score(y[idx], p[idx]))
    return float(np.mean(v)), 2 * k


def main(pred_dir, out):
    rows = []
    for code in ORDER:
        ll, rl = T.labels(code, "leaf"), T.labels(code, "root")
        j = pd.concat([ll.rename("leaf"), rl.rename("root")], axis=1, join="inner")
        j = j[j.leaf.isin([0, 1]) & j.root.isin([0, 1])]
        shared = T.shared_genes(pred_dir, code)
        j = j[j.index.isin(shared)]
        sw = j[j.leaf != j.root]
        r = dict(species=code, n_labelled_both=len(j), n_switching=len(sw),
                 pct_switching=round(100 * len(sw) / len(j), 1))
        for model_tissue in ("leaf", "root"):
            d = T.predictions(pred_dir, code, model_tissue, shared)
            col = [c for c in d.columns if c.startswith(code) and c.endswith(".h5")]
            assert len(col) == 1, (code, model_tissue, col)
            d = d[d.genes.isin(sw.index)]
            y = sw.loc[d.genes, "leaf"].to_numpy().astype(int)
            a, n = balanced_auroc(y, d[col[0]].to_numpy(), zlib.crc32(f"sw|{code}|{model_tissue}".encode()))
            r[f"{model_tissue}_model_auROC_vs_leaf_labels"] = None if np.isnan(a) else round(a, 3)
            r["n_balanced"] = n
        rows.append(r)
        print(f"  {code}: {len(sw)} switching genes of {len(j)} ({r['pct_switching']}%)")
    d = pd.DataFrame(rows)
    d["leaf_minus_root"] = (d.leaf_model_auROC_vs_leaf_labels - d.root_model_auROC_vs_leaf_labels).round(3)
    d.to_csv(out, sep="\t", index=False)
    print("\nOn the genes that switch class between tissues (auROC against the LEAF labels;")
    print("0.5 = no tissue signal, >0.5 = leaf-like, <0.5 = root-like):")
    print(d[["species", "n_switching", "pct_switching", "leaf_model_auROC_vs_leaf_labels",
             "root_model_auROC_vs_leaf_labels", "leaf_minus_root"]].to_string(index=False))
    print(f"\nmedian: leaf model {d.leaf_model_auROC_vs_leaf_labels.median():.3f}, "
          f"root model {d.root_model_auROC_vs_leaf_labels.median():.3f}, "
          f"separation {d.leaf_minus_root.median():.3f}")
    print("written:", out)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
