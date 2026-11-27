"""
Purpose:
    Align *OsACT1*, *AtACT1*, *PpACT1*, and *ZmUBI1* `gTUR` donor sequences pairwise to
    *OsACT1* (`Bio.Align.PairwiseAligner`), identify columns conserved across *OsACT1*/
    *AtACT1*/*PpACT1* inside window `[1185, 1230]`, and re-run the mutagenesis predictions.

Manuscript (Peleke et al., 2026):
    Stage 1 computation script for Supplementary Figure 5n–o (`aligned_20261109_*` tables).

Usage & I/O (TF 2.10 env; run from directory holding `lib/`, `genomes/`, `tutorial/`, `models_old/`):
    python code/05_supplementary_figures/align_donors_panel5n.0.py <zmubi1_donor_fasta> <out_prefix>
    - Reads `EXTERNAL/web_model_manifest/web_model_manifest_20261106.tsv` (`EXTERNAL_INPUTS.md`)
    - Writes `<out_prefix>_alignment.tsv`, `<out_prefix>_conserved.tsv`, `<out_prefix>_predictions.tsv`

How to change or extend:
    Pass a new donor FASTA as arg 1, or adjust window interval `[1185, 1230]` / model sets at
    the top of the script.
"""
import os
import sys

import numpy as np
import pandas as pd
import tensorflow as tf

tf.compat.v1.disable_eager_execution()
tf.compat.v1.disable_v2_behavior()
tf.config.set_visible_devices([], "GPU")

from Bio import Align  # noqa: E402
from lib.utils import one_hot_encode, one_hot_to_dna, prepare_dataset, make_predictions  # noqa: E402

GENOME = "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
ANNOT = "Arabidopsis_thaliana.TAIR10.55.gff3"
MANIFEST = ("EXTERNAL/web_model_manifest/web_model_manifest_20261106.tsv")  # EXTERNAL (web_model_manifest): see EXTERNAL_INPUTS.md
FASTA = ("tutorial/File_2_deepCRE_tutorial_promoter_sequences_"
         "Supplementary-file-2_dCRE-Pelekeetal2025.txt")
HDR = {"OsACT1": "gTUR_Osativa_OsACT1_KP100426-PIG2_5UTR500BP",
       "AtACT1": "gTUR_AthalianaTARI10_ACT1-AT2G37620_5UTR500BP",
       "PpACT1": "gTUR_Ppatens_PpACT1_Physcomitrium-patens-v6.1_5UTR500BP"}
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


def project(ref, other, aligner):
    """other projected onto ref coordinates; '-' where ref has no aligned base."""
    a = aligner.align(ref, other)[0]
    out = ["-"] * len(ref)
    for (r0, r1), (q0, q1) in zip(a.aligned[0], a.aligned[1]):
        for i in range(r1 - r0):
            out[r0 + i] = other[q0 + i]
    return "".join(out), float(a.score)


def main(zm_fasta, prefix):
    fa = read_fasta(FASTA)
    zm = read_fasta(zm_fasta)
    seqs = {k: fa[h] for k, h in HDR.items()}
    seqs["ZmUBI1"] = zm["gTUR_Zmays_ZmUBI1_NC_050100.1_5UTR500BP"]
    ref = seqs["OsACT1"]

    aligner = Align.PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                                    open_gap_score=-5, extend_gap_score=-0.5)
    proj = {"OsACT1": ref}
    for k in ("AtACT1", "PpACT1", "ZmUBI1"):
        proj[k], sc = project(ref, seqs[k], aligner)
        print(f"  {k}: alignment score {sc:.0f}, "
              f"{sum(c != '-' for c in proj[k])}/500 columns aligned")

    pd.DataFrame({"window_position": np.arange(1001, 1501),
                  **{k: list(proj[k]) for k in proj}}).to_csv(
        f"{prefix}_alignment.tsv", sep="\t", index=False)

    lo, hi = CONS_LO - 1001, CONS_HI - 1000
    cons3 = [i for i in range(500)
             if ref[i] == proj["AtACT1"][i] == proj["PpACT1"][i]]
    in_block = [i for i in cons3 if lo <= i < hi]
    rate = len(cons3) / 500
    print(f"\nwhole gTUR: {len(cons3)}/500 columns identical in OsACT1/AtACT1/PpACT1 "
          f"({100 * rate:.1f} %)")
    print(f"block {CONS_LO}-{CONS_HI} ({hi - lo} bp): {len(in_block)} conserved, "
          f"{(hi - lo) * rate:.1f} expected at the whole-gTUR rate")
    print("  conserved window positions:", [1001 + i for i in in_block])

    pd.DataFrame(dict(window_position=[1001 + i for i in cons3],
                      in_block=[lo <= i < hi for i in cons3],
                      base=[ref[i] for i in cons3])).to_csv(
        f"{prefix}_conserved.tsv", sep="\t", index=False)

    # ── rerun panel 5n on the aligned-conserved columns ─────────────────────
    x, gene_ids, chroms, starts, ends, sizes, gc, strands = prepare_dataset(
        genome=GENOME, annot=ANNOT, gene_list=GENES)
    gseq = {g: one_hot_to_dna(x[i])[0] for i, g in enumerate(gene_ids)}
    idx = {g: i for i, g in enumerate(gene_ids)}

    sm = {}
    for g in gene_ids:
        r = gseq[g]
        u = min(500, abs(ends[idx[g]] - starts[idx[g]]) // 2)
        sm[f"{g}|reference"] = r
        don = list(ref[:u])
        sm[f"{g}|OsACT1 gTUR"] = r[:1000] + "".join(don) + r[1000 + u:]
        for i in in_block:
            if i < u:
                don[i] = TRANSVERT[don[i]]
        sm[f"{g}|OsACT1 gTUR, {len(in_block)} aligned-conserved sites mutated"] = (
            r[:1000] + "".join(don) + r[1000 + u:])

    models = load_models()
    keys = list(sm)
    xx = np.array([one_hot_encode(sm[k]) for k in keys])
    rows = []
    for i, (name, path) in enumerate(models.items(), 1):
        print(f"  [{i}/{len(models)}] {name}", flush=True)
        p = dict(zip(keys, np.asarray(make_predictions(model=path, x=xx)).ravel()))
        for k in keys:
            g, var = k.split("|")
            rows.append(dict(analysis="panel5n_aligned", target=g, variant=var, model=name,
                             generation="round-3" if name.startswith("r3") else "2024",
                             P=round(float(p[k]), 4),
                             delta=np.nan if var == "reference"
                             else round(float(p[k] - p[f"{g}|reference"]), 4)))
    pd.DataFrame(rows).to_csv(f"{prefix}_predictions.tsv", sep="\t", index=False)
    print(f"\nwritten {prefix}_alignment.tsv, {prefix}_conserved.tsv, "
          f"{prefix}_predictions.tsv")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
