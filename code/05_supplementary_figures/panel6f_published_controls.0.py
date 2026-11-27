"""
Purpose:
    Score the 9 *RAP2.12* `IP-Cat-0` contributor SNPs against the exact 5 published control
    substitution sets (`ALTCtrMUT-1..5` columns of Supplementary Data 5 exported to TSV) across all
    deployed 2026 and 2024 models.

Manuscript (Peleke et al., 2026):
    Generates `panel6f_published_20261109.tsv` for Supplementary Figure 6f/6g comparisons.

Usage & I/O (TF 2.10 env; run from directory holding `lib/`, `genomes/`, `tutorial/`, `data/`, `models_old/`):
    python code/05_supplementary_figures/panel6f_published_controls.0.py <suppdata5_tsv> <out_tsv>
    - `<suppdata5_tsv>`: Supplementary Data 5 exported to TSV (`POS, REF, ALT, HAPLOTYPE, CONTRIBUTOR, ALTCtrMUT-1..5`)

How to change or extend:
    Export any variant table with `ALTCtrMUT-*` columns to TSV and pass as `<suppdata5_tsv>`.
"""
import sys

import numpy as np
import pandas as pd
import tensorflow as tf

tf.compat.v1.disable_eager_execution()
tf.compat.v1.disable_v2_behavior()
tf.config.set_visible_devices([], "GPU")

from lib.utils import one_hot_encode, one_hot_to_dna, prepare_dataset, make_predictions  # noqa: E402

import os

GENOME = "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
ANNOT = "Arabidopsis_thaliana.TAIR10.55.gff3"
RAP212 = "AT1G53910"
MANIFEST = ("EXTERNAL/web_model_manifest/web_model_manifest_20261106.tsv")  # EXTERNAL (web_model_manifest): see EXTERNAL_INPUTS.md


def load_models():
    man = pd.read_csv(MANIFEST, sep="\t")
    models = {}
    for _, r in man.iterrows():
        sp = r.species_or_set.replace("all 13 species", "MSR all")
        models[f"r3 {sp} ({r.tissue})"] = r.path
    for f in sorted(os.listdir("models_old")):
        if f.endswith(".h5"):
            sp, tis = f[:-3].rsplit("_", 1)
            models[f"2024 {sp.replace('_', ' ')} ({tis})"] = f"models_old/{f}"
    return models


def main(d5_path, out_tsv):
    d5 = pd.read_csv(d5_path, sep="\t")
    d5 = d5[(d5.HAPLOTYPE.astype(str).str.startswith("I-Cat"))
            & (d5.CONTRIBUTOR.astype(str) == "yes")].copy()
    ctr_cols = [c for c in d5.columns if c.startswith("ALTCtrMUT")]
    print(f"{len(d5)} contributor SNPs, control columns: {ctr_cols}")
    assert len(d5) == 9 and len(ctr_cols) == 5

    x, gene_ids, chroms, starts, ends, sizes, gc, strands = prepare_dataset(
        genome=GENOME, annot=ANNOT, gene_list=[RAP212])
    base = one_hot_to_dna(x[0])[0]
    strand, start, end = strands[0], starts[0], ends[0]
    u = min(500, abs(end - start) // 2)
    prom_start, term_start = start - 1000, end - u
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}

    def apply(snps):
        s, hit = base, 0
        for pos, alt in snps:
            mapped = (pos - prom_start - 1 if strand == "+"
                      else (1000 + u) - (pos - term_start - 1) - 1)
            if not 0 <= mapped < len(s):
                continue
            s = s[:mapped] + (alt if strand == "+" else comp.get(alt, "N")) + s[mapped + 1:]
            hit += 1
        return s, hit

    seqs = {"reference": base}
    seqs["nine contributors (real ALT)"], h = apply(
        list(zip(d5.POS.astype(int), d5.ALT.astype(str).str.upper())))
    print(f"  nine contributors: {h} mapped in-window")
    for c in ctr_cols:
        seqs[f"control {c[-1]}"], h = apply(
            list(zip(d5.POS.astype(int), d5[c].astype(str).str.upper())))
        print(f"  {c}: {h} mapped in-window")

    models = load_models()
    keys = list(seqs)
    xx = np.array([one_hot_encode(seqs[k]) for k in keys])
    rows = []
    for i, (name, path) in enumerate(models.items(), 1):
        print(f"  [{i}/{len(models)}] {name}", flush=True)
        p = dict(zip(keys, np.asarray(make_predictions(model=path, x=xx)).ravel()))
        for k in keys:
            rows.append(dict(analysis="panel6f", target=RAP212, variant=k, model=name,
                             generation="round-3" if name.startswith("r3") else "2024",
                             P=round(float(p[k]), 4),
                             delta=np.nan if k == "reference"
                             else round(float(p[k] - p["reference"]), 4)))
    pd.DataFrame(rows).to_csv(out_tsv, sep="\t", index=False)
    print(f"written {out_tsv} ({len(rows)} rows)")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
