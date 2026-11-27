# Supplementary Figures 5 & 6 (`code/05_supplementary_figures/`)

Inference/DeepSHAP computation scripts, pre-computed intermediate tables/arrays, and publication figure builders reproducing **Supplementary Figure 5** (cross-species promoter swaps in *A. thaliana* reporter genes + *ZmUBI1* donor + conserved-block mutagenesis) and **Supplementary Figure 6** (*RAP2.12* / `AT1G53910` 672-accession promoter haplotypes, DeepSHAP saliency, and permutation controls across 2024 and 2026 deepCRE models).

---

## Two Stages of Execution

1. **TensorFlow 2.10 Inference & DeepSHAP Stage (11 scripts):** Import `from lib.utils import ...` from https://github.com/NAMlab/streamlit_deepcre_app and read `EXTERNAL/web_model_manifest/web_model_manifest_20261106.tsv`. Run from a web-app working directory holding `lib/`, `genomes/`, `tutorial/`, `data/`, `models/`, and (for 2024 baseline comparisons) `models_old/`.
2. **Offline Figure Plotting Stage (`make_supp_figures_full.15.py`, `make_suppfig6_classic.14.py`):** Run directly with Python 3.11 + `requirements_analysis.txt` on the pre-computed `.tsv`/`.npz`/`.txt` intermediates committed in this folder + `tutorial/`.

---

## 1. Inference, Mutagenesis, and DeepSHAP Scripts (TF 2.10)

| Script | Purpose & Manuscript Role | Key CLI Arguments |
|---|---|---|
| `supp_figures_extra_data.0.py` | Runs *ZmUBI1* (`gUR`/`gTUR`/`gUR+gTUR`) swaps across all deployed 2026 & 2024 models (`suppextra_*`), *MEF1*/*RAP2.12* DeepSHAP saliency (`suppextra_*_saliency.npz`), *OsACT1* conserved-site mutations, and *IP-Cat-0* 5-control substitutions. | `<zmubi1_donor_fasta> <out_prefix>` |
| `align_donors_panel5n.0.py` | Aligns *OsACT1*, *AtACT1*, *PpACT1*, and *ZmUBI1* donors with `Bio.Align.PairwiseAligner` and identifies conserved columns in `[1185, 1230]`. | `<zmubi1_donor_fasta> <out_prefix>` |
| `fragment_swap_panel5n.0.py` | Inserts only the 46 bp *OsACT1* `gTUR` fragment (`1185–1230`) into the 9 Supplementary file 1 test genes (intact vs. conserved-site or non-conserved-site transversion controls) for **Supplementary Figure 5o**. | `<donor_alignment_fasta> <out_prefix>` |
| `panel5n_block_mutation.0.py` | Transverts all 46 bases of window positions `1185–1230` in *OsACT1* `gTUR` swaps across 28 deployed 2026 and 26 vanilla 2024 models (`panel5n_block_20261109.tsv`). | `<out_tsv>` |
| `peak_matched_panel5n.0.py` | Finds each model's own 46 bp maximum-`|saliency|` window inside *OsACT1* `gTUR` and compares whole-block mutation vs. model-matched peak mutation. | `<saliency_npz> <zmubi1_donor_fasta> <out_prefix>` |
| `panel6f_published_controls.0.py` | Evaluates the 9 *IP-Cat-0* contributor SNPs against the exact 5 published control substitution sets (`ALTCtrMUT-1..5`) in Supplementary Data 5. | `<suppdata5_tsv> <out_tsv>` |
| `rap212_saliency_icat.0.py` | Computes DeepSHAP saliency for the full 17-SNP *IP-Cat-0* *RAP2.12* haplotype (`saliency_icat_20261109.npz`) matching published **Supp. Fig. 6b**. | `<out_npz>` |
| `rap212_subset_search.0.py` | Exhaustively scores all `2^17 = 131,072` subsets of the 17 *IP-Cat-0* SNPs under 2026 and 2024 *A. thaliana* leaf models (`subsetsearch_20261109_best_by_size.tsv`). | `<out_prefix>` |
| `rap212_best10.0.py` | Scores & computes DeepSHAP saliency + 5 matched controls for the optimal 10-SNP subset identified by `rap212_subset_search.0.py`. | `<out_prefix>` |
| `rap212_permutation.0.py` | Runs 1,000-draw allele-null and position+allele-null permutation tests for the 9 contributor SNPs in `"tair10"` or `"helixer"` frame (**Supp. Fig. 6g**, `perm_20261110_*`). | `<frame> <out_prefix>` |
| `rap212_helixer_v3.0.py` | Re-runs the full 672-accession *RAP2.12* analysis in Helixer frame (`Arabidopsis_thaliana_NC_003070.9_002266`, TSS `+250 bp` vs. TAIR10) for 2026 models (`helixerv3_20261109_*`). | `<lou_fasta> <suppdata5_tsv> <out_prefix>` |

---

## 2. Figure Plotting Scripts & Pre-Computed Data

| File | Description |
|---|---|
| `make_supp_figures_full.0.py` → `make_supp_figures_full.15.py` | **Supplementary Figure 5** builder (`.15.py` is final: `5a–o`, using published alignment `donor_alignment_Supplementary-file-3.txt` in `5n` and fragment-only box plots in `5o`, outputting both `2024models` and 2026 `r3` variants). Usage: `python3 code/05_supplementary_figures/make_supp_figures_full.15.py <results_dir> <scratch_dir> <out_dir>`. |
| `make_suppfig6_classic.0.py` → `make_suppfig6_classic.14.py` | **Supplementary Figure 6** builder (`.14.py` is final: `6a–g` in classic published layout on TAIR10 windows with `g` showing the 1,000-permutation allele-null histogram). *(Note: `TFBS` intervals for Panel 6f are transcribed from PlantPAN 4.0 annotations in the published figure.)* Usage: `python3 code/05_supplementary_figures/make_suppfig6_classic.14.py <results_dir> <scratch_dir> <out_dir>`. |
| `aligned_20261109_*`, `best10_20261109_*`, `fragment_20261109_*`, `helixerv3_20261109_*`, `panel5n_block_20261109.tsv`, `panel6f_published_20261109.tsv`, `peak_20261109_*`, `perm_20261110_*`, `saliency_icat_20261109.npz`, `subsetsearch_20261109_best_by_size.tsv`, `suppextra_20261109_*`, `helixer_gene_map_20261109.tsv`, `donor_provenance_rice_blast_20261109.tsv` | Pre-computed intermediate TSV/NPZ outputs saved from Stage 1 runs above. |
| `donor_alignment_Supplementary-file-3.txt` | Published 339-column *OsACT1*/*AtACT1*/*PpACT1* `gTUR` alignment drawn in Supplementary Figure 5n. |
| `SuppFig5_full_2024models_legend.txt`, `SuppFig6_classic_2024models_legend.txt` | Generated caption/statistics text companions accompanying Supplementary Figures 5 and 6. |

---

## How to Change or Extend

* To test additional donor promoters or SNPs, run Stage 1 scripts with a custom FASTA/TSV and point `<results_dir>`/`<scratch_dir>` of `make_supp_figures_full.15.py` or `make_suppfig6_classic.14.py` to a folder containing the matching `.tsv`/`.npz` tables and `tutorial/` files.
