"""
Purpose:
    Build paralog-filtered, chromosome-stratified validation gene sets (`{species_code: [genes]}`)
    from Helixer peptide FASTAs and DIAMOND self-blastp outputs.
    Requires `pident >= --min_pident` (50%), reciprocal alignment coverage `>= --min_cov` (0.50)
    on both query and subject, and `evalue <= --max_evalue` (1e-10), selecting genes without
    inter-chromosomal paralogs strictly without tier-3 fallback (seeded, equalized across species).

Manuscript (Peleke et al., 2026):
    Produces `data/validation_set_dCRE2026.v8.p50c50e10.pickle` (Methods, training/validation splits,
    Figure 2, Supplementary Figure 1, Supplementary Table 1).

Usage & I/O (run from repo root):
    python3 code/01_data_preparation/prepare_non_homologous_validation_sets_dPHY.8.py \
        --species_config code/01_data_preparation/validation_v8/species_config_v8.json \
        --output data/validation_set_dCRE2026.v8.p50c50e10.pickle \
        [--min_pident 50.0] [--min_cov 0.5] [--max_evalue 1e-10] \
        [--n_per_species 7000] [--equalize] [--allow_tier3] [--seed 2026]

How to change or extend:
    Add new species entries to `validation_v8/species_config_v8.json` (with paths to unpacked
    Helixer `.faa`/`.fasta` and `EXTERNAL/diamond/` outputs) or tune stringency flags on the CLI.
"""
import argparse
import json
import pickle
import re

import numpy as np
import pandas as pd

CHROM_RE = re.compile(r"_(N[CW]_\d+\.\d+)_\d+$")
COLS = ["q", "s", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]


def protein_lengths(fasta):
    lengths, name, n = {}, None, 0
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip().rstrip("*"))
    if name is not None:
        lengths[name] = n
    return lengths


def gene_of(pid):
    return re.sub(r"\.\d+$", "", pid)


def eligible_genes(sp, args):
    plen = protein_lengths(sp["proteome"])
    chroms = set(sp["chromosomes"])
    genes = {}
    for pid in plen:
        g = gene_of(pid)
        m = CHROM_RE.search(g)
        if m and m.group(1) in chroms:
            genes[g] = m.group(1)
    hits = pd.read_csv(sp["diamond"], sep="\t", header=None, names=COLS, usecols=["q", "s", "pident", "length", "evalue"])
    hits["qg"], hits["sg"] = hits.q.map(gene_of), hits.s.map(gene_of)
    hits = hits[(hits.qg != hits.sg) & hits.qg.isin(genes) & hits.sg.isin(genes)]
    qcov = hits.length / hits.q.map(plen)
    scov = hits.length / hits.s.map(plen)
    par = hits[(hits.pident >= args.min_pident) & (qcov >= args.min_cov) & (scov >= args.min_cov) & (hits.evalue <= args.max_evalue)]
    par = par.assign(qc=par.qg.map(genes), sc=par.sg.map(genes))
    # symmetric: a pair is a paralog pair if either direction passes
    pairs = pd.concat([par[["qg", "sg", "qc", "sc"]],
                       par[["sg", "qg", "sc", "qc"]].set_axis(["qg", "sg", "qc", "sc"], axis=1)])
    with_par = set(pairs.qg)
    other = set(pairs.loc[pairs.qc != pairs.sc, "qg"])
    tier1 = set(genes) - with_par
    tier2 = with_par - other
    tier3 = other
    return genes, tier1, tier2, tier3


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--species_config", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--min_pident", type=float, default=50.0)
    p.add_argument("--min_cov", type=float, default=0.5)
    p.add_argument("--max_evalue", type=float, default=1e-10)
    p.add_argument("--n_per_species", type=int, default=7000)
    p.add_argument("--equalize", action="store_true", help="use the smallest eligible pool among species if below n_per_species")
    p.add_argument("--allow_tier3", action="store_true", help="v7 behaviour: fill with tier-3 genes (not recommended)")
    p.add_argument("--seed", type=int, default=2026)
    args = p.parse_args()

    species = json.load(open(args.species_config))
    pools, report = {}, []
    for sp in species:
        genes, t1, t2, t3 = eligible_genes(sp, args)
        pools[sp["id"]] = (genes, t1, t2, t3)
    n_target = args.n_per_species
    if args.equalize:
        n_target = min(n_target, min(len(t1 | t2) for _, t1, t2, _ in pools.values()))
    rng = np.random.default_rng(args.seed)
    out = {}
    for sid, (genes, t1, t2, t3) in pools.items():
        elig = t1 | t2
        df = pd.DataFrame({"gene": sorted(elig)})
        df["chrom"] = df.gene.map(genes)
        n = min(n_target, len(df))
        # proportional allocation per chromosome (largest remainder), then seeded sampling within chromosome
        share = df.chrom.value_counts() / len(df) * n
        alloc = np.floor(share).astype(int)
        rest = n - alloc.sum()
        alloc[(share - alloc).sort_values(ascending=False).index[:rest]] += 1
        chosen = []
        for chrom, k in alloc.items():
            pool = df.loc[df.chrom == chrom, "gene"].to_numpy()
            chosen += list(rng.choice(pool, size=k, replace=False))
        if args.allow_tier3 and len(chosen) < args.n_per_species:
            extra = sorted(t3 - set(chosen))
            chosen += list(rng.choice(extra, size=min(len(extra), args.n_per_species - len(chosen)), replace=False))
        out[sid] = sorted(chosen)
        cand = pd.Series(genes).value_counts()
        for chrom in sorted(cand.index):
            report.append(dict(species=sid, chrom=chrom, candidates=int(cand[chrom]),
                               tier1=sum(genes[g] == chrom for g in t1), tier2=sum(genes[g] == chrom for g in t2),
                               eligible=int((df.chrom == chrom).sum()), selected=sum(genes[g] == chrom for g in out[sid])))
    pickle.dump(out, open(args.output, "wb"))
    r = pd.DataFrame(report)
    r.to_csv(args.output + ".report.tsv", sep="\t", index=False)
    s = r.groupby("species")[["candidates", "tier1", "tier2", "eligible", "selected"]].sum()
    print(f"thresholds: pident>={args.min_pident} cov>={args.min_cov} evalue<={args.max_evalue}; target per species {n_target}")
    print(s.to_string())


if __name__ == "__main__":
    main()
