"""
Purpose:
    Insert only the 46 bp *OsACT1* `gTUR` fragment (window positions `1185–1230`) into each of
    the 9 native test genes from Supplementary file 1 — keeping all flanking bases native — and
    compare intact vs. conserved-site transversions vs. size-matched non-conserved transversions.

Manuscript (Peleke et al., 2026):
    Generates `fragment_20261109_{predictions,conserved}.tsv` plotted in Supplementary Figure 5o.

Usage & I/O (TF 2.10 env; run from directory holding `lib/`, `genomes/`, `tutorial/`, `models_old/`):
    python code/05_supplementary_figures/fragment_swap_panel5n.0.py <donor_alignment_fasta> <out_prefix>
    - Reads `EXTERNAL/web_model_manifest/web_model_manifest_20261106.tsv`
    - Writes `<out_prefix>_predictions.tsv` and `<out_prefix>_conserved.tsv`

How to change or extend:
    Adjust fragment coordinates `FRAG_LO`, `FRAG_HI` or random transversion draw counts to test
    other candidate cis-regulatory blocks.
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
HDR = {"OsACT1": "gTUR_Osativa_OsACT1_KP100426-PIG2_5UTR500BP",
       "AtACT1": "gTUR_AthalianaTARI10_ACT1-AT2G37620_5UTR500BP",
       "PpACT1": "gTUR_Ppatens_PpACT1_Physcomitrium-patens-v6.1_5UTR500BP"}
GENES = ["AT1G67090", "AT2G22200", "AT2G37620", "AT2G46800", "AT3G02150",
         "AT4G30120", "AT5G04050", "AT5G52630", "ATCG00270"]
FRAG_LO, FRAG_HI = 1185, 1230
ALIGN_OFFSET = 150                      # each donor's aligned 250 bp starts here in its gTUR
TRANSVERT = {"A": "C", "C": "A", "G": "T", "T": "G", "N": "N"}
N_CONTROL = 20


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


def conserved_positions(align_path, donor_gtur):
    """OsACT1 window positions whose alignment column is identical in all three donors."""
    raw = read_fasta(align_path)
    rows = {k.replace("250", ""): v for k, v in raw.items()}
    for k in ("OsACT1", "AtACT1", "PpACT1"):
        assert donor_gtur[k].find(rows[k].replace("-", "")) == ALIGN_OFFSET, \
            f"{k} aligned block is not at gTUR offset {ALIGN_OFFSET}"
    pos, j, out = [], 0, []
    for i, c in enumerate(rows["OsACT1"]):
        if c == "-":
            continue
        wp = 1001 + ALIGN_OFFSET + j
        j += 1
        if rows["OsACT1"][i] == rows["AtACT1"][i] == rows["PpACT1"][i]:
            out.append(wp)
    return out


def main(align_path, prefix):
    fa = read_fasta(FASTA)
    donor_gtur = {k: fa[h] for k, h in HDR.items()}
    cons_all = conserved_positions(align_path, donor_gtur)
    cons = [p for p in cons_all if FRAG_LO <= p <= FRAG_HI]
    frag_pos = list(range(FRAG_LO, FRAG_HI + 1))
    noncons = [p for p in frag_pos if p not in cons]
    print(f"{len(cons_all)} conserved columns in the alignment, "
          f"{len(cons)} inside the fragment {FRAG_LO}-{FRAG_HI} of {len(frag_pos)} bp")
    print("  conserved:", cons)

    os_g = donor_gtur["OsACT1"]
    frag = os_g[FRAG_LO - 1001:FRAG_HI - 1000]
    assert len(frag) == len(frag_pos)

    deg = list(frag)
    for p in cons:
        deg[p - FRAG_LO] = TRANSVERT[deg[p - FRAG_LO]]
    deg = "".join(deg)

    rng = np.random.default_rng(2026)
    ctrls = []
    for _ in range(N_CONTROL):
        pick = rng.choice(len(noncons), size=len(cons), replace=False)
        s = list(frag)
        for k in pick:
            p = noncons[k]
            s[p - FRAG_LO] = TRANSVERT[s[p - FRAG_LO]]
        ctrls.append("".join(s))
    print(f"  degenerated {len(cons)} conserved sites; {N_CONTROL} control sets "
          f"degenerate {len(cons)} non-conserved sites of the same fragment")

    x, gene_ids, chroms, starts, ends, sizes, gc, strands = prepare_dataset(
        genome=GENOME, annot=ANNOT, gene_list=GENES)
    gseq = {g: one_hot_to_dna(x[i])[0] for i, g in enumerate(gene_ids)}

    a, b = FRAG_LO - 1, FRAG_HI        # 0-based slice of the 3020 bp window
    sm = {}
    for g in gene_ids:
        r = gseq[g]
        sm[f"{g}|native"] = r
        sm[f"{g}|fragment {FRAG_LO}-{FRAG_HI} intact"] = r[:a] + frag + r[b:]
        sm[f"{g}|conserved sites degenerated"] = r[:a] + deg + r[b:]
        for i, c in enumerate(ctrls, 1):
            sm[f"{g}|control {i}"] = r[:a] + c + r[b:]
        for k, v in sm.items():
            assert len(v) == len(r), k

    models = load_models()
    keys = list(sm)
    xx = np.array([one_hot_encode(sm[k]) for k in keys])
    print(f"{len(keys)} sequences x {len(models)} models")
    rows = []
    for i, (name, path) in enumerate(models.items(), 1):
        print(f"  [{i}/{len(models)}] {name}", flush=True)
        p = dict(zip(keys, np.asarray(make_predictions(model=path, x=xx)).ravel()))
        for k in keys:
            g, var = k.split("|")
            rows.append(dict(analysis="panel5n_fragment", target=g, variant=var, model=name,
                             generation="round-3" if name.startswith("r3") else "2024",
                             n_sites_changed=0 if var == "native" else
                             (0 if var.endswith("intact") else len(cons)),
                             P=round(float(p[k]), 4),
                             delta=np.nan if var == "native"
                             else round(float(p[k] - p[f"{g}|native"]), 4)))
    pd.DataFrame(rows).to_csv(f"{prefix}_predictions.tsv", sep="\t", index=False)
    pd.DataFrame(dict(window_position=cons_all,
                      in_fragment=[FRAG_LO <= p <= FRAG_HI for p in cons_all],
                      base=[os_g[p - 1001] for p in cons_all])).to_csv(
        f"{prefix}_conserved.tsv", sep="\t", index=False)
    print(f"\nwritten {prefix}_predictions.tsv, {prefix}_conserved.tsv")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
