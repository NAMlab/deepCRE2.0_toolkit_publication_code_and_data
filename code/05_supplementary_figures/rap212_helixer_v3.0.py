"""
Purpose:
    Re-evaluate the full *RAP2.12* (`AT1G53910`) experiment for 2026 deepCRE models using the
    Helixer gene annotation (`Arabidopsis_thaliana_NC_003070.9_002266`, TSS +250 bp vs. TAIR10)
    across 672 Lou et al. haplotypes, variant sets, control substitutions, and DeepSHAP tracks.

Manuscript (Peleke et al., 2026):
    Generates `helixerv3_20261109_{per_accession.tsv,variants.tsv,saliency.npz,meta.tsv}`
    quantifying the *RAP2.12* annotation-frame effect noted in Supplementary Note 4 / `make_suppfig6_classic.14.py`.

Usage & I/O (TF 2.10 env; run from directory holding `lib/`, `genomes/`, `tutorial/`, `data/`):
    python code/05_supplementary_figures/rap212_helixer_v3.0.py \
        <lou_fasta> <suppdata5_tsv> <out_prefix>

How to change or extend:
    Pass any gene's haplotype alignment FASTA or modify the target Helixer locus coordinates inside `main()`.
"""
import gzip
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

GENOME = "arabidopsis_thaliana_genome_mod.fa"
ANNOT = "Arabidopsis_thaliana-GCF_000001735.4-genome.helixer_modified.gtf"
HELIXER_GENE = "Arabidopsis_thaliana_NC_003070.9_002266"
MANIFEST = ("EXTERNAL/web_model_manifest/web_model_manifest_20261106.tsv")  # EXTERNAL (web_model_manifest): see EXTERNAL_INPUTS.md
VCF = ("tutorial/File_3_deepCRE_tutorial_rap12-2-variant_"
       "Supplementary-file-3_dCRE_Peleketal2025.vcf.gz")


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


def main(lou_fa, d5_path, prefix):
    man = pd.read_csv(MANIFEST, sep="\t")
    models = {t: man[(man.species_or_set == "Arabidopsis thaliana")
                     & (man.tissue == t)].path.iloc[0] for t in ("leaf", "root")}

    x, gene_ids, chroms, starts, ends, sizes, gc, strands = prepare_dataset(
        genome=GENOME, annot=ANNOT, gene_list=[HELIXER_GENE])
    base = one_hot_to_dna(x[0])[0]
    start, end, strand = starts[0], ends[0], strands[0]
    prom_start = start - 1000
    def win(pos):
        return pos - prom_start            # 1-based window position
    print(f"{HELIXER_GENE}: start {start}, end {end}, strand {strand}; "
          f"window position of a genomic coordinate = coordinate - {prom_start}")

    # ── the 672 Lou haplotypes ──────────────────────────────────────────────
    aln = read_fasta(lou_fa)
    col = aln["Col-0"]
    keep = [i for i, c in enumerate(col) if c != "-"]
    col_block = "".join(col[i] for i in keep)
    hits = [i for i in range(len(base) - len(col_block) + 1)
            if base[i:i + len(col_block)] == col_block]
    assert len(hits) == 1, f"Col-0 block maps to {len(hits)} positions"
    off = hits[0]
    print(f"Lou alignment block at window {off + 1}-{off + len(col_block)} "
          f"({'gUR' if off + len(col_block) <= 1000 else 'spans gUR/gTUR'})")

    seqs = {}
    for acc, a in aln.items():
        s = "".join(a[i] if a[i] != "-" else "N" for i in keep)
        seqs[acc] = base[:off] + s + base[off + len(col_block):]
        assert len(seqs[acc]) == len(base)
    assert seqs["Col-0"] == base

    d5 = pd.read_csv(d5_path, sep="\t")
    icat = pd.read_csv("tutorial/Tutorial_table3_icat.csv")
    nine = d5[(d5.HAPLOTYPE.astype(str).str.startswith("I-Cat"))
              & (d5.CONTRIBUTOR.astype(str) == "yes")]
    vcf = []
    with gzip.open(VCF, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                f = line.split("\t")
                vcf.append((int(f[1]), f[4].strip().upper()))

    def apply(pairs):
        s, hit = base, 0
        for pos, alt in pairs:
            w = win(pos)
            if 1 <= w <= len(s):
                s = s[:w - 1] + alt + s[w:]
                hit += 1
        return s, hit

    variants = {"reference": (base, 0)}
    variants["IP-Cat-0 (17 SNPs)"] = apply(
        list(zip(icat.POS.astype(int), icat.ALT.astype(str).str.upper())))
    variants["nine contributors"] = apply(
        list(zip(nine.POS.astype(int), nine.ALT.astype(str).str.upper())))
    variants[f"all {len(vcf)} SNPs"] = apply(vcf)
    for c in [c for c in d5.columns if c.startswith("ALTCtrMUT")]:
        variants[f"control {c[-1]}"] = apply(
            list(zip(nine.POS.astype(int), nine[c].astype(str).str.upper())))
    for k, (s, hit) in variants.items():
        print(f"  {k}: {hit} SNPs placed in-window")

    # ── predictions ─────────────────────────────────────────────────────────
    accs = list(seqs)
    Xh = np.array([one_hot_encode(seqs[a]) for a in accs])
    vkeys = list(variants)
    Xv = np.array([one_hot_encode(variants[k][0]) for k in vkeys])
    grp = pd.read_csv("../lou2022_suppl/table_s2_groups.tsv", sep="\t").dropna(
        subset=["Accession"])
    group = dict(zip(grp.Accession, grp.group))

    rows_a, rows_v = [], []
    for tis, path in models.items():
        print(f"  predicting {tis}", flush=True)
        p = np.asarray(make_predictions(model=path, x=Xh)).ravel()
        for a, v in zip(accs, p):
            rows_a.append(dict(model=f"r3 Arabidopsis thaliana ({tis})", accession=a,
                               group=group[a], P=round(float(v), 4)))
        q = dict(zip(vkeys, np.asarray(make_predictions(model=path, x=Xv)).ravel()))
        for k in vkeys:
            rows_v.append(dict(model=f"r3 Arabidopsis thaliana ({tis})", variant=k,
                               P=round(float(q[k]), 4),
                               delta=np.nan if k == "reference"
                               else round(float(q[k] - q["reference"]), 4)))
        print(f"    reference P = {q['reference']:.4f}; "
              f"IP-Cat-0 {q['IP-Cat-0 (17 SNPs)']:.4f}; "
              f"nine {q['nine contributors']:.4f}")
    pd.DataFrame(rows_a).to_csv(f"{prefix}_per_accession.tsv", sep="\t", index=False)
    pd.DataFrame(rows_v).to_csv(f"{prefix}_variants.tsv", sep="\t", index=False)

    # ── saliency, leaf model ────────────────────────────────────────────────
    want = ["reference", f"all {len(vcf)} SNPs", "nine contributors",
            "IP-Cat-0 (17 SNPs)"]
    Xs = np.array([one_hot_encode(variants[k][0]) for k in want])
    store = {}
    for tis in ("leaf",):
        print(f"  saliency {tis}", flush=True)
        sc = compute_actual_hypothetical_scores(x=Xs, model=models[tis])
        for k, arr in zip(want, np.asarray(sc)):
            store[f"r3 Arabidopsis thaliana ({tis})||RAP2.12 {k}"] = arr.astype(np.float32)
    np.savez_compressed(f"{prefix}_saliency.npz", **store)

    meta = pd.DataFrame([
        dict(key="helixer_gene", value=HELIXER_GENE),
        dict(key="gene_start", value=start),
        dict(key="prom_start", value=prom_start),
        dict(key="block_window", value=f"{off + 1}-{off + len(col_block)}"),
        dict(key="nine_window_positions",
             value=",".join(str(win(p)) for p in nine.POS.astype(int))),
        dict(key="icat_window_positions",
             value=",".join(str(win(p)) for p in icat.POS.astype(int))),
        dict(key="vcf_window_positions",
             value=",".join(str(win(p)) for p, _ in vcf))])
    meta.to_csv(f"{prefix}_meta.tsv", sep="\t", index=False)
    print(f"\nwritten {prefix}_per_accession.tsv, {prefix}_variants.tsv, "
          f"{prefix}_saliency.npz, {prefix}_meta.tsv")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
