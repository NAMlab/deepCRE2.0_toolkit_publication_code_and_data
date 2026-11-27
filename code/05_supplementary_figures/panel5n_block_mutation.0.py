"""
Purpose:
    Transvert every base of window block `1185–1230` (`A<->C, G<->T`) inside the *OsACT1* `gTUR`
    swap across all 28 deployed 2026 models and all 26 Peleke et al. (2024) models on the 9
    Supplementary file 1 reporter genes.

Manuscript (Peleke et al., 2026):
    Generates `panel5n_block_20261109.tsv` supporting Supplementary Figure 5n–o whole-block checks.

Usage & I/O (TF 2.10 env; run from directory holding `lib/`, `genomes/`, `tutorial/`, `models_old/`):
    python code/05_supplementary_figures/panel5n_block_mutation.0.py <out_tsv>
    - Reads `EXTERNAL/web_model_manifest/web_model_manifest_20261106.tsv`

How to change or extend:
    Modify block interval boundaries `1185`–`1230` or test additional transversion rules.
"""
import os
import sys

import numpy as np
import pandas as pd
import tensorflow as tf

tf.compat.v1.disable_eager_execution()
tf.compat.v1.disable_v2_behavior()
tf.config.set_visible_devices([], "GPU")

from lib.utils import one_hot_encode, one_hot_to_dna, prepare_dataset, make_predictions  # noqa: E402

GENOME = "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
ANNOT = "Arabidopsis_thaliana.TAIR10.55.gff3"
MANIFEST = ("EXTERNAL/web_model_manifest/web_model_manifest_20261106.tsv")  # EXTERNAL (web_model_manifest): see EXTERNAL_INPUTS.md
FASTA = ("tutorial/File_2_deepCRE_tutorial_promoter_sequences_"
         "Supplementary-file-2_dCRE-Pelekeetal2025.txt")
OS_gTUR = "gTUR_Osativa_OsACT1_KP100426-PIG2_5UTR500BP"
GENES = ["AT1G67090", "AT2G22200", "AT2G37620", "AT2G46800", "AT3G02150",
         "AT4G30120", "AT5G04050", "AT5G52630", "ATCG00270"]
CONS_LO, CONS_HI = 1185, 1230
TRANSVERT = {"A": "C", "C": "A", "G": "T", "T": "G", "N": "N"}


def read_fasta(path):
    d, n = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            n = line[1:]
            d[n] = ""
        elif n:
            d[n] += line.upper()
    return d


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


def main(out_tsv):
    fa = read_fasta(FASTA)
    don = fa[OS_gTUR]
    x, gene_ids, chroms, starts, ends, sizes, gc, strands = prepare_dataset(
        genome=GENOME, annot=ANNOT, gene_list=GENES)
    seqs = {g: one_hot_to_dna(x[i])[0] for i, g in enumerate(gene_ids)}
    idx = {g: i for i, g in enumerate(gene_ids)}
    a, b = CONS_LO - 1001, CONS_HI - 1000          # donor-local slice

    sm = {}
    for g in gene_ids:
        r = seqs[g]
        u = min(500, abs(ends[idx[g]] - starts[idx[g]]) // 2)
        sm[f"{g}|reference"] = r
        base_don = don[:u]
        sm[f"{g}|OsACT1 gTUR"] = r[:1000] + base_don + r[1000 + u:]
        mut = list(base_don)
        n = 0
        for j in range(a, min(b, u)):
            mut[j] = TRANSVERT[mut[j]]
            n += 1
        sm[f"{g}|OsACT1 gTUR, block {CONS_LO}-{CONS_HI} mutated"] = (
            r[:1000] + "".join(mut) + r[1000 + u:])
    print(f"{n} of {b - a} block positions mutated per construct")

    models = load_models()
    keys = list(sm)
    xx = np.array([one_hot_encode(sm[k]) for k in keys])
    rows = []
    for i, (name, path) in enumerate(models.items(), 1):
        print(f"  [{i}/{len(models)}] {name}", flush=True)
        p = dict(zip(keys, np.asarray(make_predictions(model=path, x=xx)).ravel()))
        for k in keys:
            g, var = k.split("|")
            rows.append(dict(analysis="panel5n_block", target=g, variant=var, model=name,
                             generation="round-3" if name.startswith("r3") else "2024",
                             P=round(float(p[k]), 4),
                             delta=np.nan if var == "reference"
                             else round(float(p[k] - p[f"{g}|reference"]), 4)))
    pd.DataFrame(rows).to_csv(out_tsv, sep="\t", index=False)
    print(f"written {out_tsv} ({len(rows)} rows)")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
