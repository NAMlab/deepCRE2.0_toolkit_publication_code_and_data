"""
Purpose:
    Resolve the final set of trained single-species fold `.h5` models per species x tissue x
    chromosome by merging base training runs, rescue restarts, and `restart_log_*.tsv` status
    records against the expected chromosome fold list.

Manuscript (Peleke et al., 2026):
    Produces `code/02_training/final_models.tsv` (`231` converged folds across 13 species x 2 tissues),
    superseding interrupted partial runs (e.g. *C. sativa* OOM rerun) and tracking `converged` /
    `converged_rescue` / `not_converged` flags.

Usage & I/O:
    python3 code/02_training/s2v8_select_models.1.py \
        <model_dir> <restart_log_dir> <expected_folds.tsv> <out.tsv>
    - Reads `<model_dir>/*.h5` + `<restart_log_dir>/restart_log_*.tsv` + `<expected_folds.tsv>`
    - Writes resolved model manifest `<out.tsv>`

How to change or extend:
    Update `NAME_RE` if checkpoint filenames use a different validation tag or label code (`dPm`).
"""
import glob
import os
import re
import sys

import pandas as pd

NAME_RE = re.compile(r"^([A-Za-z]{4})_(S0X|S1X)0\.75(dC|dPm)_v8val_seed2026(_rescue\d*)?_(N[CW]_\d+\.\d+)_ssr_train_models_(\d{6}_\d{6})\.h5$")


def main(model_dir, log_dir, expected_tsv, out_tsv):
    files = []
    for f in glob.glob(os.path.join(model_dir, "*.h5")):
        m = NAME_RE.match(os.path.basename(f))
        if m:
            files.append(dict(code=m.group(1), tissue="leaf" if m.group(2) == "S0X" else "root", rescue=m.group(4) or "",
                              chrom=m.group(5), ts=m.group(6), file=f,
                              output_name=f"{m.group(1)}_{m.group(2)}0.75{m.group(3)}_v8val_seed2026{m.group(4) or ''}"))
    files = pd.DataFrame(files)
    logs = pd.concat([pd.read_csv(f, sep="\t", dtype=str) for f in glob.glob(os.path.join(log_dir, "restart_log_*.tsv"))])
    logs = logs[logs.output_name != "output_name"]
    logs["accepted"] = logs.accepted == "True"
    acc = logs[logs.accepted].drop_duplicates(["output_name", "fold"], keep="last").set_index(["output_name", "fold"])
    tried = logs.groupby(["output_name", "fold"]).size()

    base = files[files.rescue == ""]
    best_ts = base.groupby(["code", "tissue", "ts"]).size().reset_index(name="n").sort_values("n").drop_duplicates(["code", "tissue"], keep="last")
    base = base.merge(best_ts[["code", "tissue", "ts"]], on=["code", "tissue", "ts"])
    exp = pd.read_csv(expected_tsv, sep="\t")[["tissue", "code", "chrom"]].drop_duplicates()
    rows = []
    for _, e in exp.iterrows():
        b = base[(base.code == e.code) & (base.tissue == e.tissue) & (base.chrom == e.chrom)]
        r = dict(tissue=e.tissue, code=e.code, chrom=e.chrom, status="missing", file="", attempts_base=0, attempts_rescue=0,
                 train_val_loss=None, train_val_accuracy=None, train_val_auROC=None)
        if len(b):
            b = b.iloc[0]
            key = (b.output_name, e.chrom)
            r["attempts_base"] = int(tried.get(key, 0))
            if key in acc.index:
                r.update(status="converged", file=b.file)
                a = acc.loc[key]
            else:
                r.update(status="not_converged", file=b.file)
                a = None
                resc = files[(files.code == e.code) & (files.tissue == e.tissue) & (files.chrom == e.chrom) & (files.rescue != "")]
                for _, rr in resc.sort_values("rescue").iterrows():
                    k2 = (rr.output_name, e.chrom)
                    r["attempts_rescue"] += int(tried.get(k2, 0))
                    if k2 in acc.index:
                        r.update(status="converged_rescue", file=rr.file)
                        a = acc.loc[k2]
            if a is not None:
                r.update(train_val_loss=float(a.val_loss), train_val_accuracy=float(a.val_accuracy), train_val_auROC=float(a.val_auROC))
        rows.append(r)
    out = pd.DataFrame(rows)
    out.to_csv(out_tsv, sep="\t", index=False)
    print(out.groupby(["tissue", "status"]).size().unstack(fill_value=0).to_string())
    print(out[out.status.isin(["not_converged", "converged_rescue", "missing"])][["tissue", "code", "chrom", "status", "attempts_base", "attempts_rescue"]].to_string(index=False))


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
