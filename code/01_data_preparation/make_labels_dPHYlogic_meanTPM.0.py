"""
Purpose:
    Apply `classify_genes()` of `02_classify_libaries_dPHY.0.3.py` across all normalized
    species expression tables (`log10(mean TPM + 1)`) using deepPHY low-class balancing.

Manuscript (Peleke et al., 2026):
    Generates `data/<species>/<species>_<tissue>_log10_mean_tpm_Q0.75_labeled_dPHYlogic.csv.gz`
    and `data/label_summary.tsv` for all 13 species x leaf/root.

Usage & I/O (run from repo root):
    python3 code/01_data_preparation/make_labels_dPHYlogic_meanTPM.0.py <normalized_root> <out_root> [--q 0.75] [--seed 2026]
    - Reads `<normalized_root>/<species>/normalized/<species>_<tissue>_log10_mean_tpm.csv`
    - Writes labeled CSVs (`target` in {0, 1, 2}) and `<out_root>/label_summary.tsv`

How to change or extend:
    Adjust `--q` (default `0.75`) or `--seed` (default `2026`), or add tissue filenames inside `main()`.
"""
import argparse
import glob
import importlib.util
import os

import pandas as pd

BIN = "code/01_data_preparation/02_classify_libaries_dPHY.0.3.py"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("normalized_root")
    ap.add_argument("out_root")
    ap.add_argument("--q", type=float, default=0.75)
    ap.add_argument("--seed", type=int, default=2026)
    a = ap.parse_args()
    spec = importlib.util.spec_from_file_location("dphy_classify", BIN)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    rows = []
    for f in sorted(glob.glob(os.path.join(a.normalized_root, "*", "normalized", "*_log10_mean_tpm.csv"))):
        species = f.split(os.sep)[-3]
        tissue = os.path.basename(f)[len(species) + 1:].split("_log10_")[0]
        d = pd.read_csv(f)[["gene_id", "log10MeanTPM"]]
        lab, a_idx, b_idx = mod.classify_genes(d, "log10MeanTPM", a.q, a.seed)
        os.makedirs(os.path.join(a.out_root, species), exist_ok=True)
        out = os.path.join(a.out_root, species, f"{species}_{tissue}_log10_mean_tpm_Q{a.q}_labeled_dPHYlogic.csv")
        lab[["gene_id", "log10MeanTPM", "target"]].to_csv(out, index=False)
        low = lab[lab.target == 0]
        rows.append(dict(species=species, tissue=tissue, genes=len(lab), high=int((lab.target == 1).sum()), low=len(low),
                         low_zeros=int((low.log10MeanTPM == 0).sum()), low_nonzero=int((low.log10MeanTPM > 0).sum()),
                         low_max_value=round(float(low.log10MeanTPM.max()), 4), q_value=round(float(lab.log10MeanTPM.quantile(a.q)), 4),
                         zeros_total=int((lab.log10MeanTPM == 0).sum())))
    s = pd.DataFrame(rows)
    s.to_csv(os.path.join(a.out_root, "label_summary.tsv"), sep="\t", index=False)
    print(s.to_string(index=False))


if __name__ == "__main__":
    main()
