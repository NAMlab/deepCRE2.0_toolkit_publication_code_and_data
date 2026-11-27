"""
Purpose:
    Exhaustively enumerate and batch-score all `2^17 = 131,072` subsets of the 17 *IP-Cat-0*
    *RAP2.12* SNPs directly in one-hot encoding under 2024 and 2026 *A. thaliana* leaf models.

Manuscript (Peleke et al., 2026):
    Generates `subsetsearch_20261109_best_by_size.tsv` benchmarking the published 9-SNP set
    against all subsets of size `1..17`.

Usage & I/O (TF 2.10 env; run from directory holding `lib/`, `genomes/`, `data/`, `models_old/`):
    python code/05_supplementary_figures/rap212_subset_search.0.py <out_prefix>
    - Writes `<out_prefix>_best_by_size.tsv` and `<out_prefix>_subset_index.tsv`

How to change or extend:
    Adapt the SNP list for another promoter haplotype (feasible for `k <= ~20` SNPs in batch).
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd
import tensorflow as tf

tf.compat.v1.disable_eager_execution()
tf.compat.v1.disable_v2_behavior()
tf.config.set_visible_devices([], "GPU")

from lib.utils import one_hot_encode, one_hot_to_dna, prepare_dataset  # noqa: E402
from tensorflow.keras.models import load_model  # noqa: E402

GENOME = "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
ANNOT = "Arabidopsis_thaliana.TAIR10.55.gff3"
MANIFEST = ("EXTERNAL/web_model_manifest/web_model_manifest_20261106.tsv")  # EXTERNAL (web_model_manifest): see EXTERNAL_INPUTS.md
GENE = "AT1G53910"
BATCH = 1024
BASES = "ACGT"


def main(prefix):
    man = pd.read_csv(MANIFEST, sep="\t")
    r3 = man[(man.species_or_set == "Arabidopsis thaliana") & (man.tissue == "leaf")]
    models = {"deepCRE v3 Arabidopsis thaliana (leaf)": r3.path.iloc[0],
              "2024 Arabidopsis thaliana (leaf)": "models_old/Arabidopsis_thaliana_leaf.h5"}

    x, gene_ids, chroms, starts, ends, sizes, gc, strands = prepare_dataset(
        genome=GENOME, annot=ANNOT, gene_list=[GENE])
    ref_oh = np.asarray(x[0], dtype=np.float32)
    base = one_hot_to_dna(x[0])[0]
    prom_start = starts[0] - 1000

    icat = pd.read_csv("tutorial/Tutorial_table3_icat.csv")
    snps = [(int(r.POS), str(r.REF).upper(), str(r.ALT).upper(), r.CONTRIBUTOR == "yes")
            for _, r in icat.iterrows()]
    idx = [p - prom_start - 1 for p, _, _, _ in snps]
    for (p, ref_b, alt_b, _), i in zip(snps, idx):
        assert base[i] == ref_b, f"reference base at chr1:{p} is {base[i]}, expected {ref_b}"
    n = len(snps)
    published = tuple(k for k, (_, _, _, c) in enumerate(snps) if c)
    print(f"{n} IP-Cat-0 SNPs, window positions "
          f"{[i + 1 for i in idx]}; published key set = positions "
          f"{[idx[k] + 1 for k in published]}")

    alt_rows = np.zeros((n, 4), dtype=np.float32)
    for k, (_, _, alt_b, _) in enumerate(snps):
        alt_rows[k, BASES.index(alt_b)] = 1.0

    subsets = list(itertools.chain.from_iterable(
        itertools.combinations(range(n), k) for k in range(n + 1)))
    print(f"{len(subsets)} subsets to score")

    rows = []
    for name, path in models.items():
        print(f"\n=== {name}")
        model = load_model(path)
        preds = np.empty(len(subsets), dtype=np.float32)
        buf = np.repeat(ref_oh[None], BATCH, axis=0)
        for s0 in range(0, len(subsets), BATCH):
            chunk = subsets[s0:s0 + BATCH]
            b = buf[:len(chunk)]
            b[:] = ref_oh                       # reset the 17 mutable rows and the rest
            for r, sub in enumerate(chunk):
                for k in sub:
                    b[r, idx[k]] = alt_rows[k]
            preds[s0:s0 + len(chunk)] = model.predict(b, verbose=0).ravel()
            if (s0 // BATCH) % 20 == 0:
                print(f"  {s0 + len(chunk)}/{len(subsets)}", flush=True)

        p_ref = float(preds[0])                 # the empty subset
        delta = (preds - p_ref) * 100
        order = np.argsort(delta)               # most negative first
        pub_i = subsets.index(published)
        pub_d = float(delta[pub_i])
        rank = int(np.where(order == pub_i)[0][0]) + 1
        best_i = int(order[0])
        print(f"  reference P = {p_ref:.4f}")
        print(f"  published nine: {pub_d:+.2f} %, rank {rank} of {len(subsets)}")
        print(f"  best overall:   {delta[best_i]:+.2f} % with {len(subsets[best_i])} SNPs "
              f"at window {[idx[k] + 1 for k in subsets[best_i]]}")

        by_size = {}
        for i, sub in enumerate(subsets):
            k = len(sub)
            if k not in by_size or delta[i] < delta[by_size[k]]:
                by_size[k] = i
        smallest = min((len(subsets[i]) for i in range(len(subsets))
                        if delta[i] <= pub_d), default=None)
        print(f"  smallest subset reaching the published effect: {smallest} SNPs")

        for k in sorted(by_size):
            i = by_size[k]
            rows.append(dict(model=name, n_snps=k,
                             best_delta_pct=round(float(delta[i]), 3),
                             window_positions=",".join(str(idx[j] + 1) for j in subsets[i]),
                             genomic_positions=",".join(str(snps[j][0]) for j in subsets[i]),
                             is_published_set=subsets[i] == published))
        rows.append(dict(model=name, n_snps=len(published),
                         best_delta_pct=round(pub_d, 3),
                         window_positions=",".join(str(idx[j] + 1) for j in published),
                         genomic_positions=",".join(str(snps[j][0]) for j in published),
                         is_published_set=True))
        np.save(f"{prefix}_delta_{name.split()[0]}.npy", delta)

    pd.DataFrame(rows).to_csv(f"{prefix}_best_by_size.tsv", sep="\t", index=False)
    pd.DataFrame(dict(subset_index=range(len(subsets)),
                      n_snps=[len(s) for s in subsets])).to_csv(
        f"{prefix}_subset_index.tsv", sep="\t", index=False)
    print(f"\nwritten {prefix}_best_by_size.tsv")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
