"""
Purpose:
    Evaluate all 26 deployed single-species models across both leaf and root expression labels
    on the shared gene universe scored in both tissues' cross-prediction runs (important for
    *B. vulgaris*, whose leaf and root runs differ in chromosome fold count: 8 vs. 9), using
    100 CRC32-seeded balanced draws per model x target x label tissue.

Manuscript (Peleke et al., 2026):
    Generates `results/evaluation/tissue_cross/tissue_cross_{full,own_species}.tsv` (Supplementary Figure 2).
    Also acts as the shared utility module imported via `importlib` by `tissue_*.py` scripts.

Usage & I/O (run from repo root):
    python3 code/03_evaluation/tissue_cross_eval.1.py <ssr_pred_dir> <out_prefix>
    - Writes `<out_prefix>_full.tsv` (`676` rows) and `<out_prefix>_own_species.tsv`

How to change or extend:
    Pass alternative `<ssr_pred_dir>` paths or modify `N_DRAWS` / `ORDER` / `TISSUES` in this
    file (automatically inherited by downstream `tissue_*` scripts).
"""
import glob
import importlib.util
import os
import sys
import zlib

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("msreval", os.path.join(HERE, "msr_eval_dPm.2.py"))
E = importlib.util.module_from_spec(spec)
spec.loader.exec_module(E)

ORDER, SP, LAB, balanced = E.ORDER, E.SP, E.LAB, E.balanced
TISSUES = ("leaf", "root")


def shared_genes(pred_dir, code):
    """Genes a species has deployed-model predictions for in both tissues."""
    return set(_predictions(pred_dir, code, "leaf").genes) & set(_predictions(pred_dir, code, "root").genes)


def predictions(pred_dir, code, model_tissue, shared=None):
    d = _predictions(pred_dir, code, model_tissue)
    if shared is None:
        shared = shared_genes(pred_dir, code)
    return d[d.genes.isin(shared)].reset_index(drop=True)


def _predictions(pred_dir, code, model_tissue):
    """All deployed models of one tissue on the genes of one species, merged over the all12 + Hvul files."""
    files = sorted(glob.glob(os.path.join(
        pred_dir, f"deployed_{model_tissue}_*_deepcre_crosspredict_{code}_{model_tissue}.csv")))
    assert files, (code, model_tissue)
    d = None
    for f in files:
        x = pd.read_csv(f)
        x = x[["genes"] + [c for c in x.columns if c.endswith(".h5")]]
        if d is None:
            d = x
        else:
            new = [c for c in x.columns if c.endswith(".h5") and c not in d.columns]
            n0 = len(d)
            d = d.merge(x[["genes"] + new], on="genes", how="inner")
            assert len(d) == n0, (code, f, n0, len(d))
    return d


def labels(code, tissue):
    p = os.path.join(LAB, SP[code], f"{SP[code]}_{tissue}_log10_mean_tpm_Q0.75_labeled_dPHYlogic.csv")
    return pd.read_csv(p, usecols=["gene_id", "target"]).set_index("gene_id").target


def main(pred_dir, prefix):
    rows = []
    for code in ORDER:
        shared = shared_genes(pred_dir, code)
        for model_tissue in TISSUES:
            d = predictions(pred_dir, code, model_tissue, shared)
            cols = [c for c in d.columns if c.endswith(".h5")]
            assert len(cols) == len(ORDER), (code, model_tissue, len(cols))
            for label_tissue in TISSUES:
                lab = labels(code, label_tissue)
                dd = d.assign(y=d.genes.map(lab))
                dd = dd[dd.y.isin([0, 1])]
                y = dd.y.to_numpy().astype(int)
                for col in cols:
                    ms = col[:4]
                    r = dict(model_species=ms, model_tissue=model_tissue, target_species=code,
                             label_tissue=label_tissue,
                             comparison="same tissue" if model_tissue == label_tissue else "cross-tissue",
                             species="own species" if ms == code else "other species", model=col)
                    r.update(balanced(y, dd[col].to_numpy(),
                                      zlib.crc32(f"{col}|{code}|{label_tissue}".encode())))
                    rows.append(r)
            print(f"  scored {model_tissue} models on {code} ({len(dd)} labelled genes "
                  f"of {len(shared)} shared between tissues)")
    full = pd.DataFrame(rows)
    full["_m"] = full.model_species.map(ORDER.index)
    full["_t"] = full.target_species.map(ORDER.index)
    full = full.sort_values(["model_tissue", "_m", "label_tissue", "_t"]).drop(columns=["_m", "_t"])
    full.to_csv(prefix + "_full.tsv", sep="\t", index=False)

    own = full[full.model_species == full.target_species].copy()
    piv = own.pivot_table(index="model_species", columns=["model_tissue", "label_tissue"],
                          values=["auROC", "MCC"])
    piv = piv.reindex(ORDER)
    piv.round(4).to_csv(prefix + "_own_species.tsv", sep="\t")

    print("\nWithin-species: deployed model of one tissue scored on both tissues' labels")
    print("                 model=leaf              model=root")
    print("species     labels=leaf  labels=root  labels=root  labels=leaf   (auROC / MCC)")
    for c in ORDER:
        g = own[own.model_species == c].set_index(["model_tissue", "label_tissue"])
        f = lambda mt, lt: g.loc[(mt, lt)]
        print(f"{c:10s}  {f('leaf','leaf').auROC:.3f}/{f('leaf','leaf').MCC:.3f}  "
              f"{f('leaf','root').auROC:.3f}/{f('leaf','root').MCC:.3f}  "
              f"{f('root','root').auROC:.3f}/{f('root','root').MCC:.3f}  "
              f"{f('root','leaf').auROC:.3f}/{f('root','leaf').MCC:.3f}")
    for mt in TISSUES:
        s = own[own.model_tissue == mt]
        same = s[s.comparison == "same tissue"]
        cross = s[s.comparison == "cross-tissue"]
        print(f"median {mt} model: same tissue auROC {same.auROC.median():.3f} MCC {same.MCC.median():.3f} | "
              f"cross-tissue auROC {cross.auROC.median():.3f} MCC {cross.MCC.median():.3f}")
    print("\nwritten:", prefix + "_full.tsv", "and", prefix + "_own_species.tsv")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
