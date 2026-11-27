"""
Purpose:
    Compute the four intermediate datasets required to rebuild Supplementary Figures 5 & 6 across
    all deployed 2026 (`28`) and 2024 (`26`) models: (1) *ZmUBI1* `gUR`/`gTUR`/`gUR+gTUR` swaps,
    (2) *MEF1* and *RAP2.12* DeepSHAP saliency tracks, (3) *OsACT1* conserved-site mutations, and
    (4) *IP-Cat-0* 5-control random substitutions.

Manuscript (Peleke et al., 2026):
    Produces `suppextra_20261109_{predictions.tsv,saliency.npz,meta.tsv}` read by
    `make_supp_figures_full.15.py`.

Usage & I/O (TF 2.10 env; run from directory holding `lib/`, `genomes/`, `tutorial/`, `data/`, `models_old/`):
    python code/05_supplementary_figures/supp_figures_extra_data.0.py <zmubi1_donor_fasta> <out_prefix>

How to change or extend:
    Pass any donor FASTA matching *ZmUBI1* section headers (`gUR`/`gTUR`/`gUR+gTUR`) or edit
    `EXTERNAL/web_model_manifest/` to run on newly trained models.
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

GENOME = "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
ANNOT = "Arabidopsis_thaliana.TAIR10.55.gff3"
MANIFEST = ("EXTERNAL/web_model_manifest/web_model_manifest_20261106.tsv")  # EXTERNAL (web_model_manifest): see EXTERNAL_INPUTS.md
FASTA = ("tutorial/File_2_deepCRE_tutorial_promoter_sequences_"
         "Supplementary-file-2_dCRE-Pelekeetal2025.txt")
VCF = ("tutorial/File_3_deepCRE_tutorial_rap12-2-variant_"
       "Supplementary-file-3_dCRE_Peleketal2025.vcf.gz")

SWAP_GENES = ["AT1G67090", "AT2G22200", "AT2G37620", "AT2G46800", "AT3G02150",
              "AT4G30120", "AT5G04050", "AT5G52630", "ATCG00270"]
MEF1 = "AT5G52630"
RAP212 = "AT1G53910"
CONS = (1185, 1230)                      # window positions of the conserved gTUR block
TRANSVERT = {"A": "C", "C": "A", "G": "T", "T": "G", "N": "N"}

DONOR_gTUR = {
    "Os_ACT1": "gTUR_Osativa_OsACT1_KP100426-PIG2_5UTR500BP",
    "At_ACT1": "gTUR_AthalianaTARI10_ACT1-AT2G37620_5UTR500BP",
    "Pp_ACT1": "gTUR_Ppatens_PpACT1_Physcomitrium-patens-v6.1_5UTR500BP",
}
SAL_MODELS_SWAP = ["2024 Arabidopsis thaliana (leaf)", "r3 Arabidopsis thaliana (leaf)",
                   "2024 Zea mays (leaf)", "r3 Zea mays (leaf)"]
SAL_MODELS_VCF = ["2024 Arabidopsis thaliana (leaf)", "r3 Arabidopsis thaliana (leaf)"]


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


def predict_all(seq_map, models):
    keys = list(seq_map)
    x = np.array([one_hot_encode(seq_map[k]) for k in keys])
    out = {}
    for i, (name, path) in enumerate(models.items(), 1):
        print(f"    [{i}/{len(models)}] {name}", flush=True)
        out[name] = np.asarray(make_predictions(model=path, x=x)).ravel()
    return pd.DataFrame(out, index=keys)


def main(zm_fasta, prefix):
    fa = read_fasta(FASTA)
    zm = read_fasta(zm_fasta)
    zm_gUR = zm["gUR_Zmays_ZmUBI1_NC_050100.1_Promoter1000BP"]
    zm_gTUR = zm["gTUR_Zmays_ZmUBI1_NC_050100.1_5UTR500BP"]
    assert len(zm_gUR) == 1000 and len(zm_gTUR) == 500

    x, gene_ids, chroms, starts, ends, sizes, gc, strands = prepare_dataset(
        genome=GENOME, annot=ANNOT, gene_list=SWAP_GENES + [RAP212])
    seqs = {g: one_hot_to_dna(x[i])[0] for i, g in enumerate(gene_ids)}
    idx = {g: i for i, g in enumerate(gene_ids)}
    utr = {g: min(500, abs(ends[idx[g]] - starts[idx[g]]) // 2) for g in gene_ids}
    print("utr_len:", {g: utr[g] for g in gene_ids})

    models = load_models()
    print(f"{len(models)} models")
    rows = []

    # ── 1a. ZmUBI1 swaps ────────────────────────────────────────────────────
    print("\n[1a] ZmUBI1 swaps")
    sm = {}
    for g in SWAP_GENES:
        r, u = seqs[g], utr[g]
        sm[f"{g}|reference"] = r
        sm[f"{g}|ZmUBI1 gUR"] = zm_gUR + r[1000:]
        sm[f"{g}|ZmUBI1 gTUR"] = r[:1000] + zm_gTUR[:u] + r[1000 + u:]
        sm[f"{g}|ZmUBI1 gUR+gTUR"] = zm_gUR + zm_gTUR[:u] + r[1000 + u:]
    pz = predict_all(sm, models)
    for k in pz.index:
        g, var = k.split("|")
        for m in models:
            base = pz.loc[f"{g}|reference", m]
            rows.append(dict(analysis="swap", target=g, variant=var, model=m,
                             generation="round-3" if m.startswith("r3") else "2024",
                             P=round(float(pz.loc[k, m]), 4),
                             delta=np.nan if var == "reference"
                             else round(float(pz.loc[k, m] - base), 4)))

    # ── 1b. panel 5n: conserved-site mutations in the OsACT1 gTUR swap ──────
    print("\n[1b] panel 5n conserved-site mutations")
    os_g, at_g, pp_g = (fa[DONOR_gTUR[d]] for d in ("Os_ACT1", "At_ACT1", "Pp_ACT1"))
    a, b = CONS[0] - 1001, CONS[1] - 1000          # donor-local slice of the block
    cons_cols = [j for j in range(a, b)
                 if os_g[j] == at_g[j] == pp_g[j]]
    print(f"  window {CONS[0]}-{CONS[1]} = donor positions {a + 1}-{b}; "
          f"{len(cons_cols)} of {b - a} columns identical in OsACT1/AtACT1/PpACT1")
    print("  conserved window positions:", [1001 + j for j in cons_cols])
    sn = {}
    for g in SWAP_GENES:
        r, u = seqs[g], utr[g]
        don = list(os_g[:u])
        sn[f"{g}|OsACT1 gTUR"] = r[:1000] + "".join(don) + r[1000 + u:]
        for j in cons_cols:
            if j < u:
                don[j] = TRANSVERT[don[j]]
        sn[f"{g}|OsACT1 gTUR, conserved sites mutated"] = (
            r[:1000] + "".join(don) + r[1000 + u:])
    pn = predict_all(sn, models)
    for k in pn.index:
        g, var = k.split("|")
        for m in models:
            rows.append(dict(analysis="panel5n", target=g, variant=var, model=m,
                             generation="round-3" if m.startswith("r3") else "2024",
                             P=round(float(pn.loc[k, m]), 4),
                             delta=round(float(pn.loc[k, m] - pz.loc[f"{g}|reference", m]), 4)))

    # ── 1c. panel 6f: nine contributors vs five random substitutions ────────
    print("\n[1c] panel 6f random-substitution controls")
    i = idx[RAP212]
    start, end, strand, u = starts[i], ends[i], strands[i], utr[RAP212]
    prom_start, term_start = start - 1000, end - u
    base_seq = seqs[RAP212]
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}

    def apply(snps):
        s, hit = base_seq, 0
        for pos, alt in snps:
            mapped = (pos - prom_start - 1 if strand == "+"
                      else (1000 + u) - (pos - term_start - 1) - 1)
            if not 0 <= mapped < len(s):
                continue
            aa = alt.upper() if strand == "+" else comp.get(alt.upper(), "N")
            s = s[:mapped] + aa + s[mapped + 1:]
            hit += 1
        return s, hit

    icat = pd.read_csv("tutorial/Tutorial_table3_icat.csv")
    icat_9 = list(zip(icat[icat.CONTRIBUTOR == "yes"].POS.astype(int),
                      icat[icat.CONTRIBUTOR == "yes"].ALT.astype(str)))
    vcf41 = []
    with gzip.open(VCF, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                f = line.split("\t")
                vcf41.append((int(f[1]), f[4]))
    print(f"  {len(icat_9)} contributor SNPs, {len(vcf41)} VCF records")

    sv = {"reference": base_seq}
    sv["nine contributors"], _ = apply(icat_9)
    sv[f"all {len(vcf41)} SNPs"], _ = apply(vcf41)
    rng = np.random.default_rng(2026)
    for k in range(1, 6):
        rnd = []
        for pos, alt in icat_9:
            mapped = (pos - prom_start - 1 if strand == "+"
                      else (1000 + u) - (pos - term_start - 1) - 1)
            ref_b = base_seq[mapped] if 0 <= mapped < len(base_seq) else "N"
            choices = [c for c in "ACGT" if c not in (ref_b, alt.upper())]
            rnd.append((pos, str(rng.choice(choices))))
        sv[f"random substitution set {k}"], _ = apply(rnd)
    pv = predict_all(sv, models)
    for k in pv.index:
        for m in models:
            rows.append(dict(analysis="panel6f", target=RAP212, variant=k, model=m,
                             generation="round-3" if m.startswith("r3") else "2024",
                             P=round(float(pv.loc[k, m]), 4),
                             delta=np.nan if k == "reference"
                             else round(float(pv.loc[k, m] - pv.loc["reference", m]), 4)))

    pd.DataFrame(rows).to_csv(f"{prefix}_predictions.tsv", sep="\t", index=False)
    print(f"\nwritten {prefix}_predictions.tsv ({len(rows)} rows)")

    # ── 2. saliency ─────────────────────────────────────────────────────────
    print("\n[2] saliency")
    sal_seqs = {"MEF1 native": seqs[MEF1]}
    r, u = seqs[MEF1], utr[MEF1]
    for d, hdr in DONOR_gTUR.items():
        sal_seqs[f"MEF1 {d} gTUR"] = r[:1000] + fa[hdr][:u] + r[1000 + u:]
    sal_seqs["MEF1 ZmUBI1 gTUR"] = r[:1000] + zm_gTUR[:u] + r[1000 + u:]
    vcf_seqs = {"RAP2.12 reference": base_seq,
                f"RAP2.12 all {len(vcf41)} SNPs": sv[f"all {len(vcf41)} SNPs"],
                "RAP2.12 nine contributors": sv["nine contributors"]}

    store = {}
    for group, smap, mlist in (("swap", sal_seqs, SAL_MODELS_SWAP),
                               ("vcf", vcf_seqs, SAL_MODELS_VCF)):
        keys = list(smap)
        xx = np.array([one_hot_encode(smap[k]) for k in keys])
        for m in mlist:
            print(f"  {group}: {m}", flush=True)
            sc = compute_actual_hypothetical_scores(x=xx, model=models[m])
            for k, s in zip(keys, np.asarray(sc)):
                store[f"{m}||{k}"] = s.astype(np.float32)
        store[f"__seq__{group}"] = np.array(keys, dtype=object)
    np.savez_compressed(f"{prefix}_saliency.npz", **store)
    print(f"written {prefix}_saliency.npz ({len(store)} arrays)")

    meta = pd.DataFrame(dict(
        key=["conserved_window", "conserved_positions", "n_conserved",
             "zmubi1_source", "zmubi1_tss", "panel5n_scheme"],
        value=[f"{CONS[0]}-{CONS[1]}",
               ",".join(str(1001 + j) for j in cons_cols), len(cons_cols),
               os.path.abspath(zm_fasta),
               "NC_050100.1 minus strand, Zea mays GCF_902167145.1",
               "transversion A<->C, G<->T at every conserved column"]))
    meta.to_csv(f"{prefix}_meta.tsv", sep="\t", index=False)
    print(f"written {prefix}_meta.tsv")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
