"""
Purpose:
    List all models served on the web tool (`models/deployed_models.tsv`): the 26 deployed
    single-species fold models and 2 all-species MSR models (*B. oleracea* held out) with their
    validation metrics, plus the 8 Peleke et al. (2024) models (no metrics: not evaluated in this
    study), each with its path under `models/` and SHA-256 checksum.
    Version 2 of `build_web_model_manifest.1.py`: models in `models/{ssr,msr,peleke2024}/`,
    adds the 2024 models, drops the comparison with a previous online manifest.

Manuscript (Peleke et al., 2026):
    Defines the exact 28 2026 models served at https://deepcre.ipk-gatersleben.de and
    referenced via `EXTERNAL/web_model_manifest/` by `code/05_supplementary_figures/` and
    `code/06_annotation_frame/`.

Usage & I/O (run from repo root):
    python3 code/03_evaluation/build_web_model_manifest.2.py \
        models results/evaluation/final/S2v8_per_fold.tsv results/evaluation/msr models/deployed_models.tsv
    Metrics: single-species = held-out chromosome, paralog-filtered validation genes (`*_v8` columns);
    MSR = held-out species (`msr_eval_<tissue>_full.tsv`).

How to change or extend:
    Update `DEPLOY_HELD_OUT` (default `{"all": "Bole"}`), `NICE`, or `ORDER` to deploy a
    different MSR hold-out species or register newly trained species/tissues.
"""
import hashlib
import os
import sys

import pandas as pd

SP = {"Atha": "Arabidopsis thaliana", "Bvul": "Beta vulgaris", "Bole": "Brassica oleracea",
      "Csat": "Camelina sativa", "Cqui": "Chenopodium quinoa", "Cucu": "Cucumis sativus",
      "Dcar": "Daucus carota", "Gmax": "Glycine max", "Hvul": "Hordeum vulgare",
      "Osat": "Oryza sativa", "Slyc": "Solanum lycopersicum", "Sbic": "Sorghum bicolor",
      "Zmay": "Zea mays"}
ORDER = ["Atha", "Csat", "Bole", "Cucu", "Gmax", "Dcar", "Slyc", "Cqui", "Bvul", "Osat", "Hvul", "Sbic", "Zmay"]
MSR_WEB = {"all": "MSR_all_species"}
MSR_HELD_OUT = {"all": "Bole"}
MSR_SET_LABEL = {"all": "all 13 species"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(id_dir, per_fold, msr_dir, out):
    mapping = pd.read_csv(os.path.join(id_dir, "model_id_mapping.tsv"), sep="\t")
    pf = pd.read_csv(per_fold, sep="\t")
    pf = pf[pf.deployed == True]

    rows = []
    for _, r in pf.iterrows():
        m = mapping[(mapping.old_file == r.model)]
        assert len(m) == 1, (r.model, len(m))
        new_id = m.new_id.iloc[0]
        rows.append(dict(web_file=f"{SP[r.code].replace(' ', '_')}_{r.tissue}.h5", kind="single-species",
                         species_or_set=SP[r.code], tissue=r.tissue, fold_or_held_out=r.chrom, model_id=new_id,
                         path=os.path.join(id_dir, "ssr", new_id), accuracy=round(r.accuracy_v8, 4),
                         auROC=round(r.auROC_v8, 4), MCC=round(r.MCC_v8, 4),
                         _sort=(0, ORDER.index(r.code), r.tissue)))

    for tis in ("leaf", "root"):
        ev = pd.read_csv(os.path.join(msr_dir, f"msr_eval_{tis}_full.tsv"), sep="\t")
        for i, (st, web) in enumerate(MSR_WEB.items()):
            ho = MSR_HELD_OUT[st]
            e = ev[(ev["set"] == st) & (ev.status == "held_out") & (ev.held_out == ho)]
            assert len(e) == 1, (tis, st, len(e))
            e = e.iloc[0]
            m = mapping[mapping.old_file == e.model]
            assert len(m) == 1, (e.model, len(m))
            new_id = m.new_id.iloc[0]
            rows.append(dict(web_file=f"{web}_{tis}.h5", kind=f"multi-species ({st})",
                             species_or_set=MSR_SET_LABEL[st],
                             tissue=tis, fold_or_held_out=f"held out: {SP[ho]}", model_id=new_id,
                             path=os.path.join(id_dir, "msr", new_id), accuracy=round(e.accuracy, 4),
                             auROC=round(e.auROC, 4), MCC=round(e.MCC, 4), _sort=(1, i, tis)))

    p24 = pd.read_csv(os.path.join(id_dir, "peleke2024_models.tsv"), sep="\t")
    for i, r in p24[p24["set"] == "web_tool_2024"].reset_index(drop=True).iterrows():
        rows.append(dict(web_file=r.file, kind="Peleke et al. (2024)", species_or_set=r.species, tissue=r.tissue,
                         fold_or_held_out="", model_id=r.file, path=os.path.join(id_dir, "peleke2024", r.file),
                         accuracy=None, auROC=None, MCC=None, _sort=(2, i, r.tissue), _sha=r.sha256))

    d = pd.DataFrame(rows).sort_values("_sort").drop(columns="_sort").reset_index(drop=True)
    for p in d.path:
        assert os.path.exists(p), p
    d["sha256"] = [sha256(p) for p in d.path]
    known = d._sha.notna()
    assert (d.sha256[known] == d._sha[known]).all(), "2024 model checksum differs from peleke2024_models.tsv"
    d = d.drop(columns=["_sha", "model_id"])
    d["path"] = [os.path.relpath(p, id_dir) for p in d.path]
    d = d.rename(columns={"path": "file_in_models_dir"})
    d.to_csv(out, sep="\t", index=False)
    print(d.groupby("kind", sort=False).size().to_string())
    print("written:", out)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
