"""
Purpose:
    Join shuffled-sequence control restart logs (`..._ctrlSSC.tsv`) with native/rescue training
    restart logs per species x tissue x chromosome fold, emitting per-fold and per-species comparisons.

Manuscript (Peleke et al., 2026):
    Generates `results/evaluation/ssc_control/ssc_control_per_fold.tsv` (231 folds across all 13
    species including *H. vulgare*) and `ssc_control_per_species.tsv` (26 species x tissue rows),
    which enter Supplementary Data 1 (`shuffled_control` sheet) and Methods.

Usage & I/O (run from repo root):
    python3 code/03_evaluation/ssc_control_table.0.py <restart_log_dir> <out_dir>
    - Writes `<out_dir>/ssc_control_per_fold.tsv` and `<out_dir>/ssc_control_per_species.tsv`

How to change or extend:
    Update species order `ORDER` or log filename matching rules inside `main()` to include new
    species or control runs.
"""
import glob
import os
import re
import sys

import pandas as pd

ORDER = ["Atha", "Csat", "Bole", "Cucu", "Gmax", "Dcar", "Slyc", "Cqui", "Bvul", "Osat", "Hvul", "Sbic", "Zmay"]
RE_CTRL = re.compile(r"restart_log_([A-Z][a-z]{3})_(S0X|S1X)0\.75dPm_v8val_seed2026_ctrlSSC\.tsv$")


def read(path):
    return pd.read_csv(path, sep="\t")


def main(log_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    rows = []
    for f in sorted(glob.glob(os.path.join(log_dir, "restart_log_*_ctrlSSC.tsv"))):
        m = RE_CTRL.search(os.path.basename(f))
        if not m:
            continue
        code, st = m.groups()
        tissue = "leaf" if st == "S0X" else "root"
        base = os.path.join(log_dir, f"restart_log_{code}_{st}0.75dPm_v8val_seed2026")
        real = read(base + ".tsv")
        if os.path.exists(base + "_rescue.tsv"):
            resc = read(base + "_rescue.tsv")
            resc = resc[resc.accepted == True]
            real = pd.concat([real[~real.fold.isin(resc.fold)], resc], ignore_index=True)
        # one row per fold: the attempt that was kept (accepted), else the last attempt tried
        real = real.sort_values(["fold", "accepted", "attempt"]).groupby("fold", as_index=False).last()
        c = read(f).sort_values(["fold", "accepted", "attempt"]).groupby("fold", as_index=False).last()
        j = c.merge(real[["fold", "val_accuracy", "val_auROC"]], on="fold", how="left",
                    suffixes=("", "_real"))
        assert j.val_accuracy_real.notna().all(), (code, tissue, j[j.val_accuracy_real.isna()].fold.tolist())
        for _, r in j.iterrows():
            rows.append(dict(output_name=r.output_name, fold=r.fold, attempt=int(r.attempt), seed=int(r.seed),
                             val_loss=round(r.val_loss, 6), val_accuracy_shuffled=round(r.val_accuracy, 6),
                             val_auROC_shuffled=round(r.val_auROC, 6), val_auPR=round(r.val_auPR, 6),
                             accepted=bool(r.accepted), code=code, tissue=tissue,
                             val_accuracy_real=round(r.val_accuracy_real, 6),
                             val_auROC_real=round(r.val_auROC_real, 6)))
    d = pd.DataFrame(rows)
    d["_s"] = d.code.map(ORDER.index)
    d = d.sort_values(["tissue", "_s", "fold"]).drop(columns="_s")
    d.to_csv(os.path.join(out_dir, "ssc_control_per_fold.tsv"), sep="\t", index=False)

    g = d.groupby(["tissue", "code"]).agg(val_accuracy_shuffled=("val_accuracy_shuffled", "mean"),
                                          val_auROC_shuffled=("val_auROC_shuffled", "mean"),
                                          val_accuracy_real=("val_accuracy_real", "mean"),
                                          val_auROC_real=("val_auROC_real", "mean"),
                                          n_folds=("fold", "size")).round(3).reset_index()
    g["_s"] = g.code.map(ORDER.index)
    g = g.sort_values(["tissue", "_s"]).drop(columns="_s")
    g.to_csv(os.path.join(out_dir, "ssc_control_per_species.tsv"), sep="\t", index=False)

    missing = sorted(set(ORDER) - set(d.code), key=ORDER.index)
    print(f"{len(d)} folds, {d.code.nunique()} species" + (f" | MISSING {missing}" if missing else ""))
    for t in ("leaf", "root"):
        s = d[d.tissue == t]
        print(f"  {t}: shuffled accuracy {s.val_accuracy_shuffled.mean():.3f} auROC {s.val_auROC_shuffled.mean():.3f}"
              f" | real accuracy {s.val_accuracy_real.mean():.3f} auROC {s.val_auROC_real.mean():.3f}")
    print(g[g.code == "Hvul"].to_string(index=False))
    print("written:", out_dir)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
