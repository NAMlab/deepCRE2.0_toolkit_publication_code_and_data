"""
Purpose:
    Package final SSR, deployed MSR, per-fold count, MSR evaluation, and shuffled-sequence
    control evaluation TSVs from `<S2_dir>` into a formatted multi-sheet Excel workbook.
    Copies stored metrics verbatim without recomputing predictions.

Manuscript (Peleke et al., 2026):
    Produces Supplementary Data 1 (`deployed_models`, `deployed_msr_models`, `per_fold_metrics`,
    `fold_sizes`, `msr_evaluation`, `msr_summary`, `shuffled_control`; includes the two deployed
    *B. oleracea* held-out all-species MSR models on `deployed_msr_models`).

Usage & I/O (run from repo root; works offline on `results/evaluation`):
    python3 code/03_evaluation/build_supplementary_data1.1.py \
        <S2_dir> <ssc_per_fold.tsv|NONE> <out.xlsx>
    Example:
    python3 code/03_evaluation/build_supplementary_data1.1.py \
        results/evaluation results/evaluation/ssc_control/ssc_control_per_fold.tsv out.xlsx

How to change or extend:
    Update `DEPLOY_MSR`, `SP`, or `ORDER` to add species or switch which MSR hold-out models
    appear on `deployed_msr_models`. Pass `NONE` as arg 2 to omit `shuffled_control`.
"""
import os
import sys

import pandas as pd

SP = {"Atha": "Arabidopsis thaliana", "Bvul": "Beta vulgaris", "Bole": "Brassica oleracea",
      "Csat": "Camelina sativa", "Cqui": "Chenopodium quinoa", "Cucu": "Cucumis sativus",
      "Dcar": "Daucus carota", "Gmax": "Glycine max", "Hvul": "Hordeum vulgare",
      "Osat": "Oryza sativa", "Slyc": "Solanum lycopersicum", "Sbic": "Sorghum bicolor",
      "Zmay": "Zea mays"}
ORDER = ["Atha", "Csat", "Bole", "Cucu", "Gmax", "Dcar", "Slyc", "Cqui", "Bvul", "Osat", "Hvul", "Sbic", "Zmay"]
SET_NAME = {"all": "all 13 species", "dicots": "eudicots only (9 species)", "monocots": "grasses only (4 species)"}
N_TRAIN_SPECIES = {"all": 12, "dicots": 8, "monocots": 3}
DEPLOYED_SETS = ("all",)
DEPLOY_HELD_OUT = {"all": "Bole"}

# Published as Supplementary Data 1 in the final manuscript (formerly numbered Supplementary Data 6).
LEGEND = (
    "Supplementary Data 6. Evaluation of the deepCRE models of this study (labels: deepPHY logic on log10 mean TPM; "
    "validation genes: paralog-filtered list v8; seed 2026). Sheets: deployed_models (the 26 single-species models "
    "provided on the web tool), deployed_msr_models (the 2 multi-species models provided on the web tool), "
    "per_fold_metrics (all converged single-species folds), fold_sizes (training and validation genes per fold, "
    "Supplementary Table 1), msr_evaluation (every multi-species and deployed single-species model on every target "
    "species; status held_out / seen / outside, own_species / other_species), msr_summary (per target species), "
    "shuffled_control (models trained on shuffled training sequences vs the real model of the same fold; "
    "training-time validation metrics)."
)
H_DEPLOYED = ("Single-species models provided on the deepCRE web tool (one per species and tissue: the chromosome-fold "
              "model with the highest validation auROC). Metrics on the paralog-filtered validation genes of its held-out "
              "chromosome, classes balanced by random downsampling (mean of 100 draws, threshold 0.5). Means over all "
              "folds: Figure 2b, c and sheet per_fold_metrics.")
H_DEPLOYED_MSR = ("Multi-species models provided on the deepCRE web tool: the all-species model of each tissue. Each was "
                  "trained in the leave-one-species-out design on twelve of the thirteen species, with Brassica oleracea "
                  "held out to serve as the validation set for early stopping. The metrics are those of the model on that "
                  "held-out species, i.e. on a species absent from its training set (all low- and high-expression genes, "
                  "balanced, 100 draws), and are therefore an estimate of transfer to an unseen genome; for the twelve "
                  "species contained in the training set, performance is expected to be higher (sheet msr_evaluation, "
                  "status = seen). The eudicot-only and grass-only models of the same experiment are not served by the web "
                  "tool but are reported in full in msr_evaluation and msr_summary.")
H_PER_FOLD = ("Figure 2b, c. Accuracy, auROC and MCC per chromosome fold (one model per fold, trained on all other "
              "chromosomes), evaluated on the paralog-filtered validation genes of the held-out chromosome, low and high "
              "classes balanced by random downsampling (mean of 100 draws, threshold 0.5). Mean over folds and s.d. are "
              "shown in Figure 2b, c; deployed_on_web_tool = fold with the highest validation auROC.")
H_FOLD_SIZES = "Training and validation genes per fold after reduction to low/high classes and balancing."
H_MSR_EVAL = "Accuracy, auROC, MCC on all low/high genes of the target species, balanced (100 draws)."
H_MSR_SUM = "Median-type summaries per target species (MCC)."
H_SSC = "Validation accuracy and auROC at training time; real = deployed-setup model of the same fold."


def order_key(d, col):
    return d[col].map({c: i for i, c in enumerate(ORDER)})


def build(s2, ssc_path, out):
    pf = pd.read_csv(os.path.join(s2, "final", "S2v8_per_fold.tsv"), sep="\t")
    pf = pf[pf.status.isin(["converged", "converged_rescue"])].copy()
    pf["_s"] = order_key(pf, "code")
    pf = pf.sort_values(["tissue", "_s", "chrom"])

    per_fold = pd.DataFrame({
        "tissue": pf.tissue, "species": pf.code.map(SP), "fold_validation_chromosome": pf.chrom,
        "model_id": pf.model, "accuracy": pf.accuracy_v8.round(4), "auROC": pf.auROC_v8.round(4),
        "MCC": pf.MCC_v8.round(4), "n_validation_low": pf.n_low_v8.astype(int),
        "n_validation_high": pf.n_high_v8.astype(int), "n_validation_balanced": pf.n_v8_balanced.astype(int),
        "deployed_on_web_tool": pf.deployed.astype(bool)})
    deployed = per_fold[per_fold.deployed_on_web_tool].drop(
        columns=["n_validation_low", "n_validation_high", "deployed_on_web_tool"])

    fs = pd.read_csv(os.path.join(s2, "final", "SuppTable1_per_fold_counts.tsv"), sep="\t")
    fs["_s"] = order_key(fs, "code")
    fs = fs.sort_values(["tissue", "_s", "fold_chromosome"])
    fold_sizes = fs[["tissue", "species", "code", "fold_chromosome", "train_low", "train_high", "train_balanced",
                     "val_low", "val_high", "val_balanced", "deployed", "accuracy", "auROC", "MCC"]]

    msr, summ = [], []
    for t in ("leaf", "root"):
        d = pd.read_csv(os.path.join(s2, "msr", f"msr_eval_{t}_full.tsv"), sep="\t")
        d["_h"] = order_key(d, "held_out")
        d["_t"] = order_key(d, "target")
        d = d.sort_values(["kind", "set", "_h", "_t"])
        msr.append(pd.DataFrame({
            "tissue": t, "kind": d.kind, "set": d["set"], "held_out_species": d.held_out.map(SP),
            "target_species": d.target.map(SP), "status": d.status, "model_id": d.model,
            "n_balanced": d.n_balanced.astype(int), "accuracy": d.accuracy.round(4),
            "auROC": d.auROC.round(4), "MCC": d.MCC.round(4)}))
        s = pd.read_csv(os.path.join(s2, "msr", f"msr_eval_{t}_summary.tsv"), sep="\t")
        s["_t"] = order_key(s, "target")
        s = s.sort_values("_t").drop(columns="_t")
        s["target"] = s.target.map(SP)
        summ.append(s)
    msr = pd.concat(msr, ignore_index=True)
    msr_summary = pd.concat(summ, ignore_index=True)

    rows = []
    for t in ("leaf", "root"):
        for st in DEPLOYED_SETS:
            ho = DEPLOY_HELD_OUT[st]
            r = msr[(msr.tissue == t) & (msr["set"] == st) & (msr.status == "held_out") &
                    (msr.held_out_species == SP[ho])]
            assert len(r) == 1, (t, st, len(r))
            r = r.iloc[0]
            seen = msr[(msr.tissue == t) & (msr["set"] == st) & (msr.status == "seen") &
                       (msr.held_out_species == SP[ho])]
            assert len(seen) == N_TRAIN_SPECIES[st], (t, st, len(seen))
            rows.append(dict(tissue=t, training_set=SET_NAME[st], n_training_species=N_TRAIN_SPECIES[st],
                             held_out_species=SP[ho], model_id=r.model_id, accuracy=r.accuracy, auROC=r.auROC,
                             MCC=r.MCC, n_balanced=r.n_balanced,
                             mean_accuracy_training_species=round(seen.accuracy.mean(), 4),
                             mean_auROC_training_species=round(seen.auROC.mean(), 4),
                             mean_MCC_training_species=round(seen.MCC.mean(), 4)))
    deployed_msr = pd.DataFrame(rows)

    ssc = None
    if ssc_path and ssc_path.upper() != "NONE":
        c = pd.read_csv(ssc_path, sep="\t")
        c["_s"] = order_key(c, "code")
        c = c.sort_values(["tissue", "_s", "fold"])
        ssc = pd.DataFrame({"tissue": c.tissue, "species": c.code.map(SP), "fold": c.fold,
                            "val_loss": c.val_loss.round(4), "val_accuracy_shuffled": c.val_accuracy_shuffled.round(4),
                            "val_auROC_shuffled": c.val_auROC_shuffled.round(4),
                            "val_accuracy_real": c.val_accuracy_real.round(4),
                            "val_auROC_real": c.val_auROC_real.round(4)})
        miss = sorted(set(ORDER) - set(c.code), key=ORDER.index)
        if miss:
            print(f"  WARNING shuffled_control is missing: {', '.join(miss)}")
            globals()["H_SSC"] += (" Controls were trained for the " + str(len(ORDER) - len(miss)) +
                                   " species listed here; no shuffled-sequence control was trained for " +
                                   ", ".join(SP[m] for m in miss) + ".")

    sheets = [("deployed_models", H_DEPLOYED, deployed), ("deployed_msr_models", H_DEPLOYED_MSR, deployed_msr),
              ("per_fold_metrics", H_PER_FOLD, per_fold), ("fold_sizes", H_FOLD_SIZES, fold_sizes),
              ("msr_evaluation", H_MSR_EVAL, msr), ("msr_summary", H_MSR_SUM, msr_summary)]
    if ssc is not None:
        sheets.append(("shuffled_control", H_SSC, ssc))

    with pd.ExcelWriter(out, engine="openpyxl") as w:
        pd.DataFrame([[LEGEND]]).to_excel(w, sheet_name="legend", index=False, header=False)
        for name, head, df in sheets:
            pd.DataFrame([[head]]).to_excel(w, sheet_name=name, index=False, header=False)
            df.to_excel(w, sheet_name=name, index=False, startrow=2)
            print(f"  {name}: {len(df)} rows x {len(df.columns)} cols")
    print("written:", out)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    build(sys.argv[1], sys.argv[2], sys.argv[3])
