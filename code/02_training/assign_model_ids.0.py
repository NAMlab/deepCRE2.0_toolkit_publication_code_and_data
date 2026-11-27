"""
Purpose:
    Rename/hard-link trained `.h5` fold and MSR checkpoints into the standardized publication
    ID scheme (`S0X`/`S1X` single-species, `M0X`/`M1X` MSR-all, `M0M`/`M1M` MSR-clade, `0.75dPm`)
    and emit the master index table `model_id_mapping.tsv`.

Manuscript (Peleke et al., 2026):
    Generates `models/model_id_mapping.tsv` indexing all 501 archived 2026 model weights
    (`232` SSR, `217` SSC, `52` MSR) and sets the `deployed` flag for the 26 SSR web-tool models.

Usage & I/O:
    python3 code/02_training/assign_model_ids.0.py \
        <ssr_model_dir> <msr_model_dir> <final_models.tsv> <S2v8_per_fold.tsv> <out_dir>
    - Writes hard links `<out_dir>/<new_id>` and `<out_dir>/model_id_mapping.tsv`

How to change or extend:
    Extend `SSR_RE` / `MSR_RE` regexes or naming rules inside `main()` if new label tags (`dPm`)
    or training regimes are introduced.
"""
import glob
import os
import re
import sys

import pandas as pd

SSR_RE = re.compile(r"^([A-Za-z]{4})_(S0X|S1X)0\.75dPm_v8val_seed2026(_rescue\d*)?(_ctrlSSC)?_(N[CW]_\d+\.\d+)_(ssr|ssc)_train_models_(\d{6}_\d{6})\.h5$")
MSR_RE = re.compile(r"^MSR_dPm_(all|dicots|monocots)_(leaf|root)_seed2026_([A-Z][a-z]{3})_dPm_(leaf|root)_msr_train_models_(\d{6}_\d{6})\.h5$")


def main(ssr_dir, msr_dir, final_tsv, per_fold_tsv, out):
    os.makedirs(out, exist_ok=True)
    final = pd.read_csv(final_tsv, sep="\t")
    used = {os.path.basename(f): s for f, s in zip(final.file, final.status)}
    pf = pd.read_csv(per_fold_tsv, sep="\t")
    deployed = set(pf.loc[pf.deployed == True, "model"])
    rows = []
    for f in sorted(glob.glob(os.path.join(ssr_dir, "*.h5"))):
        b = os.path.basename(f)
        m = SSR_RE.match(b)
        if not m:
            continue
        code, st, resc, ctrl, chrom, case, ts = m.groups()
        new = f"{code}_{st}0.75dPm_v8{ctrl or ''}_{chrom}_{case}_train_models_{ts}.h5"
        rows.append(dict(kind="SSC control" if ctrl else "single-species", code=code, tissue="leaf" if st == "S0X" else "root",
                         fold_or_heldout=chrom, old_file=b, new_id=new, in_final_set=b in used, status=used.get(b, "control" if ctrl else "superseded/partial"),
                         deployed=b in deployed, source=f))
    for f in sorted(glob.glob(os.path.join(msr_dir, "*.h5"))):
        b = os.path.basename(f)
        m = MSR_RE.match(b)
        if not m:
            continue
        sset, tis, held, _, ts = m.groups()
        t = "0" if tis == "leaf" else "1"
        code = f"M{t}X" if sset == "all" else f"M{t}M"
        new = f"MSR_{sset}_{code}0.75dPm_{held}_msr_train_models_{ts}.h5"
        rows.append(dict(kind=f"multi-species ({sset})", code=held, tissue=tis, fold_or_heldout=held, old_file=b, new_id=new,
                         in_final_set=True, status="converged", deployed=False, source=f))
    r = pd.DataFrame(rows)
    assert not r.new_id.duplicated().any(), r[r.new_id.duplicated(keep=False)]
    for _, x in r.iterrows():
        dst = os.path.join(out, x.new_id)
        if not os.path.exists(dst):
            os.link(x.source, dst)
    r.drop(columns="source").to_csv(os.path.join(out, "model_id_mapping.tsv"), sep="\t", index=False)
    print(r.groupby(["kind", "tissue"]).size().to_string())
    print("deployed:", int(r.deployed.sum()), "| links written to", out)


if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    main(*sys.argv[1:])
