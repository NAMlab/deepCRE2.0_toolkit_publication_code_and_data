"""
Purpose:
    Compute 1,000-draw permutation tests for the 9 *IP-Cat-0* key SNPs (allele-null randomizing
    substitutions at fixed positions; position-null randomizing both position and allele in the
    aligned block) in either `"tair10"` or `"helixer"` coordinate frame.

Manuscript (Peleke et al., 2026):
    Produces `perm_20261110_{tair10,helixer}_summary.tsv` and `.npz` null arrays plotted in
    Supplementary Figure 6g (`make_suppfig6_classic.14.py`).

Usage & I/O (TF 2.10 env; run from directory holding `lib/`, `genomes/`, `data/`):
    python code/05_supplementary_figures/rap212_permutation.0.py <frame> <out_prefix>
    - `<frame>` must be `tair10` or `helixer`
    - Writes `<out_prefix>_<frame>_summary.tsv` and `.npz` permutation distributions

How to change or extend:
    Adjust permutation count `N = 1000` or random seed constants to run higher-precision nulls.
"""
import os
import re
import sys

import numpy as np
import pandas as pd
import tensorflow as tf

tf.compat.v1.disable_eager_execution()
tf.compat.v1.disable_v2_behavior()
tf.config.set_visible_devices([], "GPU")

from lib.utils import one_hot_encode, one_hot_to_dna, prepare_dataset, make_predictions  # noqa: E402

MANIFEST = ("EXTERNAL/web_model_manifest/web_model_manifest_20261106.tsv")  # EXTERNAL (web_model_manifest): see EXTERNAL_INPUTS.md
FRAMES = {
    "tair10": dict(genome="Arabidopsis_thaliana.TAIR10.dna.toplevel.fa",
                   annot="Arabidopsis_thaliana.TAIR10.55.gff3", gene="AT1G53910"),
    "helixer": dict(genome="arabidopsis_thaliana_genome_mod.fa",
                    annot="Arabidopsis_thaliana-GCF_000001735.4-genome.helixer_modified.gtf",
                    gene="Arabidopsis_thaliana_NC_003070.9_002266"),
}
# tested sets, as chromosome 1 coordinates
SETS = {
    "published nine": [20135013, 20135057, 20135073, 20135077, 20135082,
                       20135096, 20135143, 20135162, 20135212],
    "best ten (v3, TAIR10 frame)": [20135013, 20135022, 20135057, 20135062, 20135073,
                                    20135077, 20135082, 20135162, 20135212, 20135236],
    "best eight (v3, Helixer frame, leaf)": [20135013, 20135022, 20135057, 20135060,
                                             20135073, 20135162, 20135212, 20135236],
    "best eight (v3, Helixer frame, root)": [20135013, 20135022, 20135060, 20135073,
                                             20135077, 20135212, 20135233, 20135236],
}
N_PERM = 1000
SEED = 2026
BATCH = 256


def main(frame, prefix):
    cfg = FRAMES[frame]
    man = pd.read_csv(MANIFEST, sep="\t")
    models = {}
    if frame == "tair10":
        models["2024 Arabidopsis thaliana (leaf)"] = "models_old/Arabidopsis_thaliana_leaf.h5"
        models["2024 Arabidopsis thaliana (root)"] = "models_old/Arabidopsis_thaliana_root.h5"
    for t in ("leaf", "root"):
        r = man[(man.species_or_set == "Arabidopsis thaliana") & (man.tissue == t)]
        models[f"r3 Arabidopsis thaliana ({t})"] = r.path.iloc[0]

    x, gene_ids, chroms, starts, ends, sizes, gc, strands = prepare_dataset(
        genome=cfg["genome"], annot=cfg["annot"], gene_list=[cfg["gene"]])
    base = one_hot_to_dna(x[0])[0]
    prom_start = starts[0] - 1000
    def win(pos):
        return pos - prom_start                     # 1-based window position

    icat = pd.read_csv("tutorial/Tutorial_table3_icat.csv")
    alt_of = {int(r.POS): str(r.ALT).upper() for _, r in icat.iterrows()}

    # the aligned 250 bp block of Lou et al., located in this frame
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
    aln = read_fasta("lou2022_data_s1.fas")
    col = aln["Col-0"]
    blk = "".join(c for c in col if c != "-")
    hits = [i for i in range(len(base) - len(blk) + 1) if base[i:i + len(blk)] == blk]
    assert len(hits) == 1
    block = list(range(hits[0] + 1, hits[0] + len(blk) + 1))     # 1-based window positions
    print(f"[{frame}] window anchored at {starts[0]}; aligned block at "
          f"{block[0]}-{block[-1]}")

    def build(pairs):
        s = base
        for w, b in pairs:
            s = s[:w - 1] + b + s[w:]
        return s

    rng = np.random.default_rng(SEED)
    rows = []
    for set_name, positions in SETS.items():
        wpos = [win(p) for p in positions]
        if not all(1 <= w <= len(base) for w in wpos):
            print(f"  {set_name}: outside the window in this frame, skipped")
            continue
        real = build([(win(p), alt_of[p]) for p in positions])
        k = len(positions)

        allele_null, position_null = [], []
        for _ in range(N_PERM):
            pr = []
            for p in positions:
                w = win(p)
                choices = [c for c in "ACGT" if c not in (base[w - 1], alt_of[p])]
                pr.append((w, str(rng.choice(choices))))
            allele_null.append(build(pr))
            ws = rng.choice(block, size=k, replace=False)
            pr = []
            for w in ws:
                choices = [c for c in "ACGT" if c != base[int(w) - 1]]
                pr.append((int(w), str(rng.choice(choices))))
            position_null.append(build(pr))

        seqs = [base, real] + allele_null + position_null
        X = np.array([one_hot_encode(s) for s in seqs])
        for mname, path in models.items():
            preds = np.concatenate([
                np.asarray(make_predictions(model=path, x=X[i:i + BATCH])).ravel()
                for i in range(0, len(X), BATCH)])
            p_ref, p_real = float(preds[0]), float(preds[1])
            a_null = (preds[2:2 + N_PERM] - p_ref) * 100
            p_null = (preds[2 + N_PERM:] - p_ref) * 100
            obs = (p_real - p_ref) * 100
            for nm, null in (("allele null", a_null), ("position null", p_null)):
                pv = (int((null <= obs).sum()) + 1) / (N_PERM + 1)
                rows.append(dict(frame=frame, model=mname, snp_set=set_name, n_snps=k,
                                 null=nm, reference_P=round(p_ref, 4),
                                 observed_delta_pct=round(obs, 3),
                                 null_mean_pct=round(float(null.mean()), 3),
                                 null_sd_pct=round(float(null.std(ddof=1)), 3),
                                 null_min_pct=round(float(null.min()), 3),
                                 n_perm=N_PERM, p_empirical=pv))
                print(f"  {mname:34s} {set_name:38s} {nm:14s} "
                      f"obs {obs:+7.2f} %  null {null.mean():+6.2f} "
                      f"+-{null.std(ddof=1):.2f}  p = {pv:.4f}")
            slug = re.sub(r"[^a-z0-9]+", "-", set_name.lower()).strip("-")
            mslug = re.sub(r"[^a-z0-9]+", "-", mname.lower()).strip("-")
            np.savez_compressed(
                f"{prefix}_{frame}__{mslug}__{slug}.npz",
                allele_null=a_null, position_null=p_null, observed=obs, reference=p_ref)

    pd.DataFrame(rows).to_csv(f"{prefix}_{frame}_summary.tsv", sep="\t", index=False)
    print(f"\nwritten {prefix}_{frame}_summary.tsv")


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] not in FRAMES:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
