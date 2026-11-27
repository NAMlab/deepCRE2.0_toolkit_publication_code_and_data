# Arabidopsis TAIR10 Annotation-Frame Exception (`code/02_training/tair10_exception/`)

Training configuration, translated expression labels, and cluster submission script for training *A. thaliana* single-species models directly on curated TAIR10/RefSeq gene windows rather than Helixer gene models (Supplementary Note 4 & **Run F** in `code/06_annotation_frame/`).

---

## Purpose & Finding in Peleke et al. (2026)

All 13 species' primary deployed models use Helixer gene models (`*.helixer.agat.gtf`) for a consistent TSS/TES frame across genomes. As quantified in [`code/06_annotation_frame/`](../../06_annotation_frame/README.md) (`annotation_frame_summary_20261110.tsv`), evaluating Helixer-trained *A. thaliana* models on curated TAIR10 windows drops auROC from `0.799` to `0.701` (`-0.098`, leaf) and `0.817` to `0.735` (`-0.083`, root) because Helixer under-calls 5'UTRs in compact, UTR-rich *A. thaliana* (median TSS `30 bp` downstream in `64%` of genes), whereas *O. sativa*, *S. lycopersicum*, and *Z. mays* show differences between `-0.026` and `+0.024`.

These configs train companion TAIR10-window *A. thaliana* models using the same hyperparameters, random seed (`2026`), and fold chromosomes (`NC_003070.9`–`NC_003076.8`). Because labels are transferred 1:1 from Helixer genes (`24,600 / 27,456` TAIR10 genes = `89.6%`; `12,981 / 14,940` validation genes under key `atha_tair10`), only the input sequence windows change and the primary Helixer-frame benchmark (Figure 2, Supplementary Figure 1, Supplementary Data 1) remains unchanged.

---

## Files in This Folder

| File | Description |
|---|---|
| `arabidopsis_thaliana_leaf_log10_mean_tpm_Q0.75_labeled_dPHYlogic_TAIR10.csv` | Leaf deepPHY binary labels translated to `AT*G*` TAIR10 locus identifiers. |
| `arabidopsis_thaliana_root_log10_mean_tpm_Q0.75_labeled_dPHYlogic_TAIR10.csv` | Root deepPHY binary labels translated to `AT*G*` TAIR10 locus identifiers. |
| `01_Atha_TAIR10_leaf.json` | Training config (`output_name`: `Atha_TAIR10_S0X0.75dPm_v8_leaf_seed2026`, `pickle_key`: `atha_tair10`). |
| `02_Atha_TAIR10_root.json` | Training config (`output_name`: `Atha_TAIR10_S0X0.75dPm_v8_root_seed2026`, `pickle_key`: `atha_tair10`). |
| `job_train_atha_tair10.0.sh` | Original SGE 2-task array job script referencing `$REPO/EXTERNAL/environment/` and `$REPO/EXTERNAL/cluster_scratch/`. |

Additional inputs referenced by the JSON configs:
* `data/arabidopsis_thaliana/genome/GCF_000001735.4.fna` (fetched by `bash data/download_genomes.sh`)
* `data/validation_set_dCRE2026.v8.p50c50e10.tair10.pickle` (key `atha_tair10`)
* Protein-coding RefSeq GTF (`Arabidopsis_thaliana_TAIR10_RefSeq.gene.gtf`, derived from `EXTERNAL/refseq_gff/GCF_000001735.4_TAIR10.1_genomic.gff.gz` and staged at `code/06_benchmarking/workflows/deepcre/src/deepCRE/inputs/tair10_exception_20261110/`)

---

## Usage & How to Extend

```bash
# Stage train_models_seed_restart.2.py alongside train_models.py in code/06_benchmarking/workflows/deepcre/src/deepCRE/,
# stage Arabidopsis_thaliana_TAIR10_RefSeq.gene.gtf at the path in the JSON (or edit "annotation"), then run:
python train_models_seed_restart.2.py \
  --input code/02_training/tair10_exception/01_Atha_TAIR10_leaf.json \
  --seed 2026 --max_attempts 5 --plateau_loss 0.68
```

Evaluation of the resulting models against Helixer-frame and deployed models is performed in `code/06_annotation_frame/annotation_frame_runF.0.py`. To repeat this pattern for another curated annotation (e.g., Ensembl/RefSeq in another species), prepare a gene-level GTF, map expression labels and validation pickle keys to the new gene IDs, and copy one of the JSON configs above.
