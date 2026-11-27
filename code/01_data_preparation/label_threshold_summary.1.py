"""
Purpose:
    Format the summary table produced by `label_threshold_summary.0.py` into a side-by-side
    comparison of old (2024) vs. new (2026) expression class boundaries per species and tissue.

Manuscript (Peleke et al., 2026):
    Generates `results/label_thresholds_20261126/label_threshold_old_vs_new_20261126.tsv` (Table 2).

Usage & I/O (run from repo root):
    python3 code/01_data_preparation/label_threshold_summary.1.py <summary_tsv_from_.0> <out_tsv>
    - Reads output TSV from `label_threshold_summary.0.py`
    - Writes formatted comparison table `<out_tsv>`

How to change or extend:
    Adjust column formatting or rounding inside `main()` to add additional summary columns
    from `.0.py`.
"""
import sys
import pandas as pd

s = pd.read_csv(sys.argv[1], sep="\t")
ab = {"arabidopsis_thaliana": "A. thaliana", "beta_vulgaris": "B. vulgaris", "brassica_oleracea": "B. oleracea",
      "camelina_sativa": "Ca. sativa", "chenopodium_quinoa": "C. quinoa", "cucumis_sativus": "Cu. sativus",
      "daucus_carota": "D. carota", "glycine_max": "G. max", "hordeum_vulgare": "H. vulgare", "oryza_sativa": "O. sativa",
      "solanum_lycopersicum": "S. lycopersicum", "sorghum_bicolor": "S. bicolor", "zea_mays": "Z. mays"}
o = pd.DataFrame({
    "species": s.species.map(ab), "tissue": s.tissue, "genes": s.n_genes, "zeroTPM_%": s.pct_zeroTPM,
    "high_thr_TPM": s.P75_TPM,
    "high_n_old": s.old_n_high, "high_n_new": s.new_n_high,
    "low_max_TPM_old": s.P25_TPM, "low_max_TPM_new": s.new_low_max_TPM,
    "low_max_pctile_old": 25.0, "low_max_pctile_new": s.new_low_max_percentile,
    "low_n_old": s.old_n_low, "low_n_new": s.new_n_low,
    "low_zero_%_old": s.old_low_pct_zero, "low_zero_%_new": s.new_low_pct_zero,
    "low/high_old": (s.old_n_low / s.old_n_high).round(2), "low/high_new": (s.new_n_low / s.new_n_high).round(2),
})
o = o.sort_values(["tissue", "species"])
o.to_csv(sys.argv[2], sep="\t", index=False)
print(o.to_markdown(index=False))
