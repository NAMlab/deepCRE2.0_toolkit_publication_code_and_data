"""
Purpose:
    Locate each model's own 46 bp peak of summed `|DeepSHAP contribution|` inside *OsACT1* `gTUR`
    on *MEF1* and compare whole-block (`1185–1230`) mutation vs. model-matched peak mutation across
    all 9 reporter genes.

Manuscript (Peleke et al., 2026):
    Generates `peak_20261109_{peaks,predictions}.tsv` checking whether model-specific saliency peaks
    outperform the fixed `1185–1230` block across species/models.

Usage & I/O (TF 2.10 env; run from directory holding `lib/`, `genomes/`, `tutorial/`, `models_old/`):
    python code/05_supplementary_figures/peak_matched_panel5n.0.py \
        <saliency_npz> <zmubi1_donor_fasta> <out_prefix>

How to change or extend:
    Pass a custom `<saliency_npz>` (e.g. from `supp_figures_extra_data.0.py`) or vary the 46 bp
    sliding window width.
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
PUB_LO, PUB_HI = 1185, 1230
WIN = PUB_HI - PUB_LO + 1
TRANSVERT = {"A": "C", "C": "A", "G": "T", "T": "G", "N": "N"}
ALN = "aligned_20261109_alignment.tsv"


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


def conserved_in(aln, lo, hi):
    b = aln[(aln.window_position >= lo) & (aln.window_position <= hi)]
    return int(((b.OsACT1 == b.AtACT1) & (b.OsACT1 == b.PpACT1)).sum())


def mutate(don, lo, hi, u):
    """Transvert window positions lo..hi inside a donor gTUR of length u."""
    s = list(don[:u])
    n = 0
    for p in range(lo, hi + 1):
        j = p - 1001
        if 0 <= j < u:
            s[j] = TRANSVERT[s[j]]
            n += 1
    return "".join(s), n


def main(npz_path, zm_fasta, prefix):
    sal = np.load(npz_path, allow_pickle=True)
    aln = pd.read_csv(ALN, sep="\t")
    don = read_fasta(FASTA)[OS_gTUR]

    peaks = {}
    rows_pk = []
    for k in sal.files:
        if "||MEF1 Os_ACT1 gTUR" not in k:
            continue
        m = k.split("||")[0]
        g = sal[k].sum(axis=1)[1000:1500]
        cs = np.convolve(g, np.ones(WIN), "valid")
        i = int(np.argmax(cs))
        lo, hi = 1001 + i, 1000 + i + WIN
        peaks[m] = (lo, hi)
        rows_pk.append(dict(model=m, peak_lo=lo, peak_hi=hi,
                            peak_sum=round(float(cs[i]), 4),
                            offset_from_published=lo - PUB_LO,
                            conserved_in_peak=conserved_in(aln, lo, hi),
                            conserved_in_published=conserved_in(aln, PUB_LO, PUB_HI),
                            score_at_published=round(float(cs[PUB_LO - 1001]), 4)))
    pk = pd.DataFrame(rows_pk).sort_values("model")
    pk.to_csv(f"{prefix}_peaks.tsv", sep="\t", index=False)
    print(pk.to_string(index=False))

    x, gene_ids, chroms, starts, ends, sizes, gc, strands = prepare_dataset(
        genome=GENOME, annot=ANNOT, gene_list=GENES)
    gseq = {g: one_hot_to_dna(x[i])[0] for i, g in enumerate(gene_ids)}
    idx = {g: i for i, g in enumerate(gene_ids)}
    models = load_models()

    rows = []
    for m, (lo, hi) in peaks.items():
        sm = {}
        for g in gene_ids:
            r = gseq[g]
            u = min(500, abs(ends[idx[g]] - starts[idx[g]]) // 2)
            sm[f"{g}|reference"] = r
            sm[f"{g}|OsACT1 gTUR"] = r[:1000] + don[:u] + r[1000 + u:]
            pub, _ = mutate(don, PUB_LO, PUB_HI, u)
            sm[f"{g}|published block {PUB_LO}-{PUB_HI} mutated"] = (
                r[:1000] + pub + r[1000 + u:])
            own, _ = mutate(don, lo, hi, u)
            sm[f"{g}|own peak {lo}-{hi} mutated"] = r[:1000] + own + r[1000 + u:]
        keys = list(sm)
        xx = np.array([one_hot_encode(sm[k]) for k in keys])
        print(f"  scoring {m} (peak {lo}-{hi})", flush=True)
        p = dict(zip(keys, np.asarray(make_predictions(model=models[m], x=xx)).ravel()))
        for k in keys:
            g, var = k.split("|")
            rows.append(dict(analysis="panel5n_peak", target=g, variant=var, model=m,
                             generation="round-3" if m.startswith("r3") else "2024",
                             peak_lo=lo, peak_hi=hi,
                             P=round(float(p[k]), 4),
                             delta=np.nan if var == "reference"
                             else round(float(p[k] - p[f"{g}|reference"]), 4)))
    pd.DataFrame(rows).to_csv(f"{prefix}_predictions.tsv", sep="\t", index=False)
    print(f"\nwritten {prefix}_peaks.tsv, {prefix}_predictions.tsv")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
