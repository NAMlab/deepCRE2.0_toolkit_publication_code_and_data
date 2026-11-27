"""
Purpose:
    Score and compute DeepSHAP contribution tracks (plus 5 random non-ref/non-alt controls)
    for the optimal 10-SNP *IP-Cat-0* *RAP2.12* subset discovered by `rap212_subset_search.0.py`.

Manuscript (Peleke et al., 2026):
    Produces `best10_20261109_{predictions.tsv,saliency.npz,positions.tsv}` for Supplementary Figure 6
    follow-up checks under the 2026 *A. thaliana* leaf model.

Usage & I/O (TF 2.10 env; run from directory holding `lib/`, `genomes/`, `data/`, `models_old/`):
    python code/05_supplementary_figures/rap212_best10.0.py <out_prefix>

How to change or extend:
    Replace `BEST10_WINDOW` at the top with any subset output from `rap212_subset_search.0.py`.
"""
import os
import sys

import numpy as np
import pandas as pd
import tensorflow as tf

tf.compat.v1.disable_eager_execution()
tf.compat.v1.disable_v2_behavior()
tf.config.set_visible_devices([], "GPU")

from lib.utils import (one_hot_encode, one_hot_to_dna, prepare_dataset,  # noqa: E402
                       make_predictions, compute_actual_hypothetical_scores)

GENOME = "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
ANNOT = "Arabidopsis_thaliana.TAIR10.55.gff3"
MANIFEST = ("EXTERNAL/web_model_manifest/web_model_manifest_20261106.tsv")  # EXTERNAL (web_model_manifest): see EXTERNAL_INPUTS.md
GENE = "AT1G53910"
# window positions of the deepCRE v3 optimum, from subsetsearch_20261109_best_by_size.tsv
BEST10_WINDOW = [1170, 1179, 1214, 1219, 1230, 1234, 1239, 1319, 1369, 1393]
N_CONTROL = 5


def main(prefix):
    man = pd.read_csv(MANIFEST, sep="\t")
    r3 = man[(man.species_or_set == "Arabidopsis thaliana") & (man.tissue == "leaf")]
    models = {"r3 Arabidopsis thaliana (leaf)": r3.path.iloc[0],
              "2024 Arabidopsis thaliana (leaf)": "models_old/Arabidopsis_thaliana_leaf.h5"}

    x, gene_ids, chroms, starts, ends, sizes, gc, strands = prepare_dataset(
        genome=GENOME, annot=ANNOT, gene_list=[GENE])
    base = one_hot_to_dna(x[0])[0]
    prom_start = starts[0] - 1000

    icat = pd.read_csv("tutorial/Tutorial_table3_icat.csv")
    alt_of, ref_of = {}, {}
    for _, r in icat.iterrows():
        w = int(r.POS) - prom_start - 1 + 1            # 1-based window position
        alt_of[w], ref_of[w] = str(r.ALT).upper(), str(r.REF).upper()
    missing = [w for w in BEST10_WINDOW if w not in alt_of]
    assert not missing, f"no IP-Cat-0 allele for window {missing}"
    print("deepCRE v3 optimum:", ", ".join(
        f"{w} {ref_of[w]}>{alt_of[w]} (chr1:{w + prom_start})" for w in BEST10_WINDOW))

    def apply(pairs):
        s = base
        for w, b in pairs:
            s = s[:w - 1] + b + s[w:]
        return s

    seqs = {"reference": base,
            "best ten (v3)": apply([(w, alt_of[w]) for w in BEST10_WINDOW])}
    rng = np.random.default_rng(2026)
    for k in range(1, N_CONTROL + 1):
        pairs = []
        for w in BEST10_WINDOW:
            choices = [c for c in "ACGT" if c not in (ref_of[w], alt_of[w])]
            pairs.append((w, str(rng.choice(choices))))
        seqs[f"control {k}"] = apply(pairs)
    for k, v in seqs.items():
        assert len(v) == len(base), k

    keys = list(seqs)
    xx = np.array([one_hot_encode(seqs[k]) for k in keys])
    rows = []
    for name, path in models.items():
        p = dict(zip(keys, np.asarray(make_predictions(model=path, x=xx)).ravel()))
        for k in keys:
            rows.append(dict(analysis="best10", target=GENE, variant=k, model=name,
                             generation="round-3" if name.startswith("r3") else "2024",
                             P=round(float(p[k]), 4),
                             delta=np.nan if k == "reference"
                             else round(float(p[k] - p["reference"]), 4)))
        print(f"  {name}: reference {p['reference']:.4f}, "
              f"best ten {p['best ten (v3)']:.4f} "
              f"({100 * (p['best ten (v3)'] - p['reference']):+.2f} %)")
    pd.DataFrame(rows).to_csv(f"{prefix}_predictions.tsv", sep="\t", index=False)

    sal_seqs = {"RAP2.12 best ten (v3)": seqs["best ten (v3)"]}
    xs = np.array([one_hot_encode(v) for v in sal_seqs.values()])
    store = {}
    for name, path in models.items():
        print(f"  saliency {name}", flush=True)
        sc = compute_actual_hypothetical_scores(x=xs, model=path)
        for k, arr in zip(sal_seqs, np.asarray(sc)):
            store[f"{name}||{k}"] = arr.astype(np.float32)
    np.savez_compressed(f"{prefix}_saliency.npz", **store)

    pd.DataFrame(dict(window_position=BEST10_WINDOW,
                      genomic_position=[w + prom_start for w in BEST10_WINDOW],
                      ref=[ref_of[w] for w in BEST10_WINDOW],
                      alt=[alt_of[w] for w in BEST10_WINDOW])).to_csv(
        f"{prefix}_positions.tsv", sep="\t", index=False)
    print(f"\nwritten {prefix}_predictions.tsv, {prefix}_saliency.npz, "
          f"{prefix}_positions.tsv")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
