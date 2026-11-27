"""
Purpose:
    Compute DeepSHAP contribution arrays for the full 17-SNP *IP-Cat-0* *RAP2.12* promoter
    haplotype under both the 2024 and 2026 *A. thaliana* leaf models.

Manuscript (Peleke et al., 2026):
    Produces `saliency_icat_20261109.npz` plotted alongside Supplementary Figure 6b–d.

Usage & I/O (TF 2.10 env; run from directory holding `lib/`, `genomes/`, `tutorial/`, `data/`, `models_old/`):
    python code/05_supplementary_figures/rap212_saliency_icat.0.py <out_npz>
    - Writes compressed arrays keyed as `<model>||<variant>` to `<out_npz>`

How to change or extend:
    Add variant keys or additional model names from `EXTERNAL/web_model_manifest/`.
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
                       compute_actual_hypothetical_scores)

GENOME = "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
ANNOT = "Arabidopsis_thaliana.TAIR10.55.gff3"
MANIFEST = ("EXTERNAL/web_model_manifest/web_model_manifest_20261106.tsv")  # EXTERNAL (web_model_manifest): see EXTERNAL_INPUTS.md
GENE = "AT1G53910"
MODELS = {"2024 Arabidopsis thaliana (leaf)": "models_old/Arabidopsis_thaliana_leaf.h5"}


def main(out_npz):
    man = pd.read_csv(MANIFEST, sep="\t")
    r3 = man[(man.species_or_set == "Arabidopsis thaliana") & (man.tissue == "leaf")]
    MODELS["r3 Arabidopsis thaliana (leaf)"] = r3.path.iloc[0]

    x, gene_ids, chroms, starts, ends, sizes, gc, strands = prepare_dataset(
        genome=GENOME, annot=ANNOT, gene_list=[GENE])
    base = one_hot_to_dna(x[0])[0]
    start, end, strand = starts[0], ends[0], strands[0]
    u = min(500, abs(end - start) // 2)
    prom_start = start - 1000

    icat = pd.read_csv("tutorial/Tutorial_table3_icat.csv")
    snps = list(zip(icat.POS.astype(int), icat.ALT.astype(str).str.upper()))
    s, hit = base, 0
    for pos, alt in snps:
        mapped = pos - prom_start - 1
        if 0 <= mapped < len(s):
            s = s[:mapped] + alt + s[mapped + 1:]
            hit += 1
    print(f"IP-Cat-0: {len(snps)} SNPs listed, {hit} mapped, "
          f"{sum(a != b for a, b in zip(base, s))} positions changed")

    seqs = {"RAP2.12 IP-Cat-0 17 SNPs": s}
    xx = np.array([one_hot_encode(v) for v in seqs.values()])
    store = {}
    for name, path in MODELS.items():
        print(f"  {name}", flush=True)
        sc = compute_actual_hypothetical_scores(x=xx, model=path)
        for k, arr in zip(seqs, np.asarray(sc)):
            store[f"{name}||{k}"] = arr.astype(np.float32)
    np.savez_compressed(out_npz, **store)
    print(f"written {out_npz} ({len(store)} arrays)")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
