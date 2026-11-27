"""
Purpose:
    Compare the Peleke et al. (2024) expression threshold definition (`low <= P25, high >= P75`)
    against the 2026 deepPHY rule (`high > P75`, `low` matched count with `<= 50%` zero-TPM)
    computed directly from `log10(mean TPM + 1)` label files and TPM equivalents (`10**x - 1`).

Manuscript (Peleke et al., 2026):
    Produces the intermediate per-species/tissue threshold summary consumed by
    `label_threshold_summary.1.py` for Table 2 and Methods.

Usage & I/O (run from repo root):
    python3 code/01_data_preparation/label_threshold_summary.0.py <labels_root> <out_tsv>
    - Reads `<labels_root>/<sp>/<sp>_<tissue>_log10_mean_tpm_Q0.75_labeled_dPHYlogic.csv`
    - Writes `<out_tsv>` with quantile cut-offs, counts, zero-TPM shares, and TPM equivalents.

How to change or extend:
    Pass a custom `<labels_root>` holding unpacked species label CSVs, or edit `SPECIES` /
    `TISSUES` if new organisms or tissues are added.
"""
import glob
import os
import sys

import numpy as np
import pandas as pd


def tpm(x):
    return 10 ** x - 1


root, out = sys.argv[1], sys.argv[2]
rows = []
for f in sorted(glob.glob(os.path.join(root, "labels", "*", "*_labeled_dPHYlogic.csv"))):
    sp = f.split(os.sep)[-2]
    tissue = os.path.basename(f)[len(sp) + 1:].split("_log10_")[0]
    d = pd.read_csv(f)
    v = d.log10MeanTPM.to_numpy()
    p25, p75 = np.percentile(v, 25), np.percentile(v, 75)
    n_zero = int((v == 0).sum())
    # old rule (Peleke et al. 2024, DeepCRE/model/utils.py): low <= P25, high >= P75
    old_low, old_high = int((v <= p25).sum()), int((v >= p75).sum())
    old_low_zero = int(((v <= p25) & (v == 0)).sum())
    # new rule as stored in the deployed label files
    hi, lo = d[d.target == 1].log10MeanTPM, d[d.target == 0].log10MeanTPM
    lo_nz = lo[lo > 0]
    # sanity: high class must equal "> P75" (pandas quantile == numpy percentile, linear)
    assert int((v > p75).sum()) == len(hi), f"{f}: high class != >P75"
    rows.append(dict(
        species=sp, tissue=tissue, n_genes=len(v), n_zeroTPM=n_zero, pct_zeroTPM=round(100 * n_zero / len(v), 1),
        P25_log=round(p25, 4), P25_TPM=round(tpm(p25), 3), P75_log=round(p75, 4), P75_TPM=round(tpm(p75), 2),
        old_n_high=old_high, old_n_low=old_low, old_low_pct_zero=round(100 * old_low_zero / old_low, 1) if old_low else np.nan,
        old_low_max_TPM=round(tpm(p25), 3),
        new_n_high=len(hi), new_n_low=len(lo), new_low_n_zero=int((lo == 0).sum()), new_low_n_nonzero=len(lo_nz),
        new_low_pct_zero=round(100 * (lo == 0).mean(), 1),
        new_low_max_log=round(lo.max(), 4), new_low_max_TPM=round(tpm(lo.max()), 3),
        new_low_max_percentile=round(100 * (v <= lo.max()).mean(), 1),
        new_n_excluded=int((d.target == 2).sum()),
        old_low_genes_dropped_by_new=int(((v <= p25) & (d.target != 0)).sum()),
    ))
s = pd.DataFrame(rows)
s.to_csv(out, sep="\t", index=False)
pd.set_option("display.width", 250)
print(s.to_string(index=False))
num = s.select_dtypes("number")
print("\nAcross 26 species x tissue (median [min-max]):")
for c in ["pct_zeroTPM", "P25_TPM", "P75_TPM", "old_n_low", "old_low_pct_zero", "new_n_low", "new_low_pct_zero",
          "new_low_max_TPM", "new_low_max_percentile"]:
    print(f"  {c:24s} {num[c].median():9.3f} [{num[c].min():.3f}-{num[c].max():.3f}]")
print("  P25 == 0 in", int((s.P25_log == 0).sum()), "of", len(s))
