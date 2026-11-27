"""
Purpose:
    For genes switching class between leaf and root (`higher in leaf` vs. `higher in root`),
    compute mean predicted probabilities for each species' deployed leaf and root SSR models,
    their difference `P(leaf model) - P(root model)`, and paired 10,000-permutation tests.

Manuscript (Peleke et al., 2026):
    Generates `results/evaluation/tissue_cross/tissue_switch_direction.tsv` and
    `tissue_switch_permutation.tsv` plotted in Supplementary Figure 2k–n.

Usage & I/O (run from repo root):
    python3 code/03_evaluation/tissue_switch_direction.1.py <ssr_pred_dir> <out.tsv>

How to change or extend:
    Adjust permutation count or summary statistics inside `main()`.
"""
import importlib.util
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("tcross", os.path.join(HERE, "tissue_cross_eval.1.py"))
T = importlib.util.module_from_spec(spec)
spec.loader.exec_module(T)
ORDER = T.ORDER


def main(pred_dir, out):
    rows = []
    for code in ORDER:
        ll, rl = T.labels(code, "leaf"), T.labels(code, "root")
        j = pd.concat([ll.rename("leaf"), rl.rename("root")], axis=1, join="inner")
        j = j[j.leaf.isin([0, 1]) & j.root.isin([0, 1])]
        shared = T.shared_genes(pred_dir, code)
        j = j[j.index.isin(shared)]
        groups = {"higher in leaf": j[(j.leaf == 1) & (j.root == 0)].index,
                  "higher in root": j[(j.leaf == 0) & (j.root == 1)].index}
        p = {}
        for mt in ("leaf", "root"):
            d = T.predictions(pred_dir, code, mt, shared)
            col = [c for c in d.columns if c.startswith(code) and c.endswith(".h5")]
            assert len(col) == 1, (code, mt, col)
            p[mt] = d.set_index("genes")[col[0]]
        for g, idx in groups.items():
            pl = p["leaf"].reindex(idx).dropna()
            pr = p["root"].reindex(idx).dropna()
            assert len(pl) == len(pr), (code, g)
            rows.append(dict(species=code, gene_group=g, n_genes=len(idx), n_predicted=len(pl),
                             n_labelled_both_tissues=len(j),
                             pct_of_labelled=round(100 * len(idx) / len(j), 2),
                             mean_P_leaf_model=round(float(pl.mean()), 4),
                             mean_P_root_model=round(float(pr.mean()), 4),
                             contrast_leaf_minus_root=round(float(pl.mean() - pr.mean()), 4)))
        print(f"  {code}: {len(groups['higher in leaf'])} higher in leaf, "
              f"{len(groups['higher in root'])} higher in root, of {len(j)} labelled in both")
    d = pd.DataFrame(rows)
    d["_s"] = d.species.map(ORDER.index)
    d = d.sort_values(["_s", "gene_group"]).drop(columns="_s")
    d.to_csv(out, sep="\t", index=False)
    for g in ("higher in leaf", "higher in root"):
        s = d[d.gene_group == g]
        print(f"median contrast, genes {g}: {s.contrast_leaf_minus_root.median():+.3f}")
    piv = d.pivot(index="species", columns="gene_group", values="contrast_leaf_minus_root").reindex(ORDER)
    ok = int(((piv["higher in leaf"] > 0) & (piv["higher in root"] < 0)).sum())
    sep = (piv["higher in leaf"] - piv["higher in root"])
    print(f"expected sign in both groups: {ok}/{len(ORDER)} species; "
          f"separation between the two groups positive in {int((sep > 0).sum())}/{len(ORDER)}")
    print("written:", out)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
