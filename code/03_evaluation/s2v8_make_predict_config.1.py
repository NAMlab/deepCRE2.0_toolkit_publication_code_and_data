"""
Purpose:
    Generate the `deepcre_crosspredict.py` JSON configuration predicting every converged/final
    single-species fold model in `final_models.tsv` on its own held-out chromosome using the
    deepPHY mean-TPM label targets.

Manuscript (Peleke et al., 2026):
    Creates `code/03_evaluation/own_chromosome_configs/own_chromosome_predict.json` executed
    by `run_s2v8_predict.0.sh` before running `s2v8_validation_metrics.2.py`.

Usage & I/O (run from repo root):
    python3 code/03_evaluation/s2v8_make_predict_config.1.py \
        <final_models.tsv> <out.json> <pred_out_dir>
    - Expects unpacked species inputs at `EXTERNAL/data_root/` (`EXTERNAL_INPUTS.md`)

How to change or extend:
    Point `E` / `MSR_LEAF` to local `data/` + `code/02_training/` paths or update `SP_INFO`
    when adding new species.
"""
import json
import os
import sys

import pandas as pd

E = "EXTERNAL/data_root/MSR_eval_inputs_20261008"  # EXTERNAL (data_root): see EXTERNAL_INPUTS.md
G = "data/"
SP = {"Atha": "arabidopsis_thaliana", "Bvul": "beta_vulgaris", "Bole": "brassica_oleracea", "Csat": "camelina_sativa",
      "Cqui": "chenopodium_quinoa", "Cucu": "cucumis_sativus", "Dcar": "daucus_carota", "Gmax": "glycine_max",
      "Osat": "oryza_sativa", "Slyc": "solanum_lycopersicum", "Sbic": "sorghum_bicolor", "Zmay": "zea_mays",
      "Hvul": "hordeum_vulgare"}
LAB = "data/"
MSR_LEAF = "EXTERNAL/predictions/msr_crosspredict_leaf_20261008.json"  # EXTERNAL (predictions): see EXTERNAL_INPUTS.md


def inputs():
    out = {}
    for e in json.load(open(MSR_LEAF)):
        sp = e["genome"].split("/")[-2]
        out.setdefault(sp, (e["genome"], e["annotation"]))
    out["hordeum_vulgare"] = (f"{E}/genome/hlx-Hordeum_vulgare_supsp_vulgare-GCF_904849725.1-112509.fna",
                              f"{E}/annotation/hordeum_vulgare/hlx-Hordeum_vulgare_supsp_vulgare-GCF_904849725.1-112509.agat.gtf")
    return out


def main(final_tsv, out_json, pred_dir):
    fm = pd.read_csv(final_tsv, sep="\t")
    io = inputs()
    cfg = []
    for _, r in fm.iterrows():
        if r.status == "missing":
            print("missing, skipped:", r.tissue, r.code, r.chrom)
            continue
        sp = SP[r.code]
        genome, annot = io[sp]
        targets = f"{LAB}/{sp}/{sp}_{r.tissue}_log10_mean_tpm_Q0.75_labeled_dPHYlogic.csv"
        for f in (genome, annot, targets, r.file):
            assert os.path.exists(f), f
        base = os.path.basename(r.file)[:-3]
        cfg.append(dict(genome=genome, annotation=annot, targets=targets, prediction_models=[r.file], chromosomes=[r.chrom],
                        ignore_small_genes=True, output_path=os.path.join(pred_dir, f"{base}_deepcre_crosspredict_{r.code}_S2v8.csv")))
    json.dump(cfg, open(out_json, "w"), indent=2)
    print(len(cfg), "fold models to predict")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
