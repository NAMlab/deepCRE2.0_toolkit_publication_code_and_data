"""
Purpose:
    Categorize every gene scored in both leaf and root into the four inter-tissue label
    transitions (`stable`, `swapped`, `high to intermediate`, `low to intermediate`) per species.

Manuscript (Peleke et al., 2026):
    Generates `results/evaluation/tissue_cross/tissue_label_composition.tsv` plotted in
    Supplementary Figure 2j.

Usage & I/O (run from repo root):
    python3 code/03_evaluation/tissue_label_composition.1.py <ssr_pred_dir> <out.tsv>

How to change or extend:
    Imports gene set helpers from `tissue_cross_eval.1.py`; edit category definitions in `main()`
    to inspect finer expression bins.
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
        shared = T.shared_genes(pred_dir, code)
        lab = {}
        for t in ("leaf", "root"):
            s = T.labels(code, t)
            lab[t] = s[s.index.isin(shared) & s.isin([0, 1])]
        both = pd.concat([lab["leaf"].rename("l"), lab["root"].rename("r")], axis=1, join="inner")
        rows.append(dict(
            species=code, n_predicted_both_tissues=len(shared),
            n_classed_in_both=len(both),
            n_high_in_both=int(((both.l == 1) & (both.r == 1)).sum()),
            n_low_in_both=int(((both.l == 0) & (both.r == 0)).sum()),
            n_high_leaf_low_root=int(((both.l == 1) & (both.r == 0)).sum()),
            n_low_leaf_high_root=int(((both.l == 0) & (both.r == 1)).sum()),
            n_classed_in_leaf_only=len(set(lab["leaf"].index) - set(lab["root"].index)),
            n_classed_in_root_only=len(set(lab["root"].index) - set(lab["leaf"].index)),
            # a gene classed in one tissue only moved to the intermediate range in the other; the
            # manuscript groups these by the class it left, not by the tissue it kept
            n_high_to_intermediate=int(sum(
                len(set(lab[t].index[lab[t] == 1]) - set(lab[o].index)) for t, o in
                (("leaf", "root"), ("root", "leaf")))),
            n_low_to_intermediate=int(sum(
                len(set(lab[t].index[lab[t] == 0]) - set(lab[o].index)) for t, o in
                (("leaf", "root"), ("root", "leaf"))))))
        print(f"  {code}: {len(both)} classed in both, "
              f"{rows[-1]['n_high_leaf_low_root'] + rows[-1]['n_low_leaf_high_root']} swapped, "
              f"{rows[-1]['n_high_to_intermediate']} high->intermediate, "
              f"{rows[-1]['n_low_to_intermediate']} low->intermediate")
    d = pd.DataFrame(rows)
    d["n_classed_in_at_least_one"] = (d.n_classed_in_both + d.n_classed_in_leaf_only
                                      + d.n_classed_in_root_only)
    d["n_swapped"] = d.n_high_leaf_low_root + d.n_low_leaf_high_root
    d["n_stable"] = d.n_high_in_both + d.n_low_in_both
    for c in ("n_stable", "n_swapped", "n_high_to_intermediate", "n_low_to_intermediate"):
        d["pct_" + c[2:]] = (100 * d[c] / d.n_classed_in_at_least_one).round(2)
    d["n_flipped"] = d.n_swapped
    d["pct_flipped_of_classed_in_both"] = (100 * d.n_flipped / d.n_classed_in_both).round(2)
    d["pct_flipped_of_at_least_one"] = (100 * d.n_flipped / d.n_classed_in_at_least_one).round(2)
    d.to_csv(out, sep="\t", index=False)
    tot = d.sum(numeric_only=True)
    print(f"\ntotal: {int(tot.n_classed_in_both):,} classed in both "
          f"({int(tot.n_high_in_both):,} high in both, {int(tot.n_low_in_both):,} low in both, "
          f"{int(tot.n_flipped):,} flipped), {int(tot.n_classed_in_leaf_only):,} leaf only, "
          f"{int(tot.n_classed_in_root_only):,} root only")
    print(f"swapped share: {d.pct_flipped_of_classed_in_both.median():.1f}% of the genes classed in both, "
          f"{d.pct_flipped_of_at_least_one.median():.1f}% of those classed in at least one tissue")
    print("pooled %% of genes classed in at least one tissue: stable %.1f, swapped %.1f, "
          "high->intermediate %.1f, low->intermediate %.1f" % tuple(
              100 * tot[c] / tot.n_classed_in_at_least_one for c in
              ("n_stable", "n_swapped", "n_high_to_intermediate", "n_low_to_intermediate")))
    print("written:", out)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
