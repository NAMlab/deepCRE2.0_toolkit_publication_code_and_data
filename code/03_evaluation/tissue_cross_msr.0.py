"""
Purpose:
    Evaluate multi-species (MSR) leaf and root models on each held-out species against both
    leaf and root expression labels and measure leaf-minus-root probability contrasts on genes
    that switch class between tissues (restricted to the shared leaf/root gene set).

Manuscript (Peleke et al., 2026):
    Generates `results/evaluation/tissue_cross/tissue_cross_msr_{full,switch_direction}.tsv`
    for Supplementary Figure 2 (panels f–i, m–n).

Usage & I/O (run from repo root):
    python3 code/03_evaluation/tissue_cross_msr.0.py <msr_pred_dir> <ssr_pred_dir> <out_prefix>
    - Uses `<ssr_pred_dir>` to determine the shared leaf/root gene set per species

How to change or extend:
    Depends on `tissue_cross_eval.1.py` (`T`) in the same folder; adjust `TISSUES` or MSR sets
    at the top if extending to additional tissues or species pools.
"""
import glob
import importlib.util
import os
import re
import sys
import zlib

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("tcross", os.path.join(HERE, "tissue_cross_eval.1.py"))
T = importlib.util.module_from_spec(spec)
spec.loader.exec_module(T)
ORDER, TISSUES = T.ORDER, ("leaf", "root")
MSR_RE = re.compile(r"^MSR_dPm_(all|dicots|monocots)_(leaf|root)_seed2026_([A-Z][a-z]{3})_dPm_(leaf|root)_"
                    r"msr_train_models_\d{6}_\d{6}\.h5$")
SETS = {"all": set(ORDER), "dicots": set(ORDER[:9]), "monocots": set(ORDER[9:])}


def msr_predictions(msr_dir, code, tissue, shared):
    f = sorted(glob.glob(os.path.join(msr_dir, f"MSR_{tissue}_deepcre_crosspredict_{code}_{tissue}.csv")))
    assert f, (code, tissue)
    d = pd.concat([pd.read_csv(x) for x in f], axis=0).drop_duplicates("genes")
    d = d[["genes"] + [c for c in d.columns if c.endswith(".h5")]]
    return d[d.genes.isin(shared)].reset_index(drop=True)


def main(msr_dir, ssr_dir, prefix):
    rows, contrast = [], []
    for code in ORDER:
        shared = T.shared_genes(ssr_dir, code)
        lab = {t: T.labels(code, t) for t in TISSUES}
        j = pd.concat([lab["leaf"].rename("leaf"), lab["root"].rename("root")], axis=1, join="inner")
        j = j[j.leaf.isin([0, 1]) & j.root.isin([0, 1])]
        j = j[j.index.isin(shared)]
        groups = {"high in leaf only": j[(j.leaf == 1) & (j.root == 0)].index,
                  "high in root only": j[(j.leaf == 0) & (j.root == 1)].index}
        p_held = {}
        for model_tissue in TISSUES:
            d = msr_predictions(msr_dir, code, model_tissue, shared)
            cols = [c for c in d.columns if MSR_RE.match(c)]
            assert cols, (code, model_tissue)
            for col in cols:
                m = MSR_RE.match(col)
                sset, held = m.group(1), m.group(3)
                if held != code:                      # only the model that held this species out
                    continue
                p_held[(sset, model_tissue)] = d.set_index("genes")[col]
                for label_tissue in TISSUES:
                    dd = d.assign(y=d.genes.map(lab[label_tissue]))
                    dd = dd[dd.y.isin([0, 1])]
                    r = dict(kind="MSR", msr_set=sset, held_out_species=code, model_tissue=model_tissue,
                             label_tissue=label_tissue,
                             comparison="same tissue" if model_tissue == label_tissue else "cross-tissue",
                             model=col)
                    r.update(T.balanced(dd.y.to_numpy().astype(int), dd[col].to_numpy(),
                                        zlib.crc32(f"{col}|{code}|{label_tissue}".encode())))
                    rows.append(r)
        for sset in sorted({s for s, _ in p_held}):
            if (sset, "leaf") not in p_held or (sset, "root") not in p_held:
                continue
            for g, idx in groups.items():
                pl = p_held[(sset, "leaf")].reindex(idx).dropna()
                pr = p_held[(sset, "root")].reindex(idx).dropna()
                contrast.append(dict(msr_set=sset, held_out_species=code, gene_group=g, n_genes=len(pl),
                                     mean_P_leaf_model=round(float(pl.mean()), 4),
                                     mean_P_root_model=round(float(pr.mean()), 4),
                                     contrast_leaf_minus_root=round(float(pl.mean() - pr.mean()), 4)))
        print(f"  {code}: {len(j)} genes labelled in both and predicted, "
              f"{len(groups['high in leaf only'])} high in leaf only, "
              f"{len(groups['high in root only'])} high in root only")

    f = pd.DataFrame(rows)
    f["_s"] = f.held_out_species.map(ORDER.index)
    f = f.sort_values(["msr_set", "model_tissue", "_s", "label_tissue"]).drop(columns="_s")
    f.to_csv(prefix + "_full.tsv", sep="\t", index=False)
    c = pd.DataFrame(contrast)
    c["_s"] = c.held_out_species.map(ORDER.index)
    c = c.sort_values(["msr_set", "_s", "gene_group"]).drop(columns="_s")
    c.to_csv(prefix + "_switch_direction.tsv", sep="\t", index=False)

    print("\nMSR models on their held-out species, median over species:")
    for sset in ("all", "dicots", "monocots"):
        s = f[f.msr_set == sset]
        if s.empty:
            continue
        a = s[s.comparison == "same tissue"]
        b = s[s.comparison == "cross-tissue"]
        print(f"  {sset:9s} n={len(a):3d}  same tissue MCC {a.MCC.median():.3f} auROC {a.auROC.median():.3f}"
              f" | cross-tissue MCC {b.MCC.median():.3f} auROC {b.auROC.median():.3f}"
              f" | loss {a.MCC.median()-b.MCC.median():+.3f}")
    print("\nGenes high in one tissue only, MSR contrast (leaf model - root model):")
    for sset in ("all", "dicots", "monocots"):
        s = c[c.msr_set == sset]
        if s.empty:
            continue
        for g in ("high in leaf only", "high in root only"):
            v = s[s.gene_group == g]
            print(f"  {sset:9s} {g:20s} median {v.contrast_leaf_minus_root.median():+.3f}  "
                  f"(expected sign in {int((v.contrast_leaf_minus_root > 0).sum() if 'leaf' in g else (v.contrast_leaf_minus_root < 0).sum())}/{len(v)})")
    print("\nwritten:", prefix + "_full.tsv", "and", prefix + "_switch_direction.tsv")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
