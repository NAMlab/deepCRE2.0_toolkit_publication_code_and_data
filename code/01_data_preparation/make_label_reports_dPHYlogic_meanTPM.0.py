"""
Purpose:
    Verify stored binary `target` labels against `classify_genes()` (allowing boundary ties)
    and generate diagnostic threshold text reports and PNG histograms for all species and tissues.

Manuscript (Peleke et al., 2026):
    Produces every report and histogram under `data/<species>/label_reports/{reports,plots}/`.

Usage & I/O (run from repo root):
    python3 code/01_data_preparation/make_label_reports_dPHYlogic_meanTPM.0.py data [--q 0.75] [--seed 2026]
    - Reads `data/<species>/<species>_<tissue>_log10_mean_tpm_Q0.75_labeled_dPHYlogic.csv.gz`
    - Writes `*_classification_report.txt` and `*_classification_histogram.png` under `label_reports/`

How to change or extend:
    Pass a custom data directory as `data_root`, override `--q` / `--seed`, or extend `TISSUES`.
"""
import argparse
import glob
import importlib.util
import os

import pandas as pd

BIN = "code/01_data_preparation/02_classify_libaries_dPHY.0.3.py"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root")
    ap.add_argument("--q", type=float, default=0.75)
    ap.add_argument("--seed", type=int, default=2026)
    a = ap.parse_args()
    spec = importlib.util.spec_from_file_location("dphy_classify", BIN)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for f in sorted(glob.glob(os.path.join(a.data_root, "*", f"*_log10_mean_tpm_Q{a.q}_labeled_dPHYlogic.csv.gz"))):
        species = os.path.basename(os.path.dirname(f))
        tissue = os.path.basename(f)[len(species) + 1:].split("_log10_")[0]
        stored = pd.read_csv(f)
        d = stored[["gene_id", "log10MeanTPM"]]
        lab, a_idx, b_idx = mod.classify_genes(d, "log10MeanTPM", a.q, a.seed)
        diff = lab[mod.TARGET_COL].values != stored["target"].values
        low = stored[stored["target"] == 0]
        boundary = low.loc[low["log10MeanTPM"] > 0, "log10MeanTPM"].max()
        bad = diff & (stored["log10MeanTPM"].values != boundary)
        if bad.any() or sorted(lab[mod.TARGET_COL].value_counts().items()) != sorted(stored["target"].value_counts().items()):
            raise SystemExit(f"{f}: {int(diff.sum())} genes differ from the stored labels, {int(bad.sum())} not at a tie")
        st = stored[["gene_id", "log10MeanTPM", "target"]].rename(columns={"target": mod.TARGET_COL})
        a_st = low.index[low["log10MeanTPM"] == 0].tolist()
        b_st = low.index[low["log10MeanTPM"] > 0].tolist()
        mod.generate_report_and_plot(st, "log10MeanTPM", species, f"{tissue}_log10_mean_tpm", a.q, a_st, b_st,
                                     os.path.join(a.data_root, species, "label_reports"), True)
        note = f", {int(diff.sum())} genes swapped within a tie at {boundary:.6f}" if diff.any() else ""
        print(f"{species} {tissue}: classes reproduced ({len(stored)} genes{note})")


if __name__ == "__main__":
    main()
