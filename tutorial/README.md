# Web-Tool Tutorial Files (`tutorial/`)

Exact copies of the tutorial inputs (`Supplementary files 1–3` of Peleke et al., 2026) and exported tutorial output tables from `https://github.com/NAMlab/streamlit_deepcre_app`.

## Contents

| File | Manuscript Role & Contents |
|---|---|
| `File_1_deepCRE_tutorial_gene_IDs_Supplementary-file-1_dCRE-Pelekeetal2025.txt` | **Supplementary file 1:** 9 *A. thaliana* gene IDs used as reporter targets in the web-tool tutorial and in all Supplementary Figure 5 promoter-swap experiments. |
| `File_2_deepCRE_tutorial_promoter_sequences_Supplementary-file-2_dCRE-Pelekeetal2025.txt` | **Supplementary file 2:** Donor regulatory sequences (*AtACT1*, *OsACT1*, *PpACT1*, *CaMV 35S*) tested in the Supplementary Figure 5 `gUR` / `gTUR` / `gUR+gTUR` swap assays. |
| `File_3_deepCRE_tutorial_rap12-2-variant_Supplementary-file-3_dCRE_Peleketal2025.vcf.gz` | **Supplementary file 3:** VCF of *RAP2.12* (`AT1G53910`) promoter variants in accession `IP-Cat-0` analyzed in Supplementary Figure 6. |
| `Tutorial_table1_Atleaf.csv` | Exported web-tool predictions for *A. thaliana* leaf tutorial run. |
| `Tutorial_table2_Atroot.csv` | Exported web-tool predictions for *A. thaliana* root tutorial run. |
| `Tutorial_table3_icat.csv` | Position/allele/contributor annotations (`POS`, `CONTRIBUTOR`) for the *IP-Cat-0* *RAP2.12* SNPs read by `code/05_supplementary_figures/` scripts (`make_suppfig6_classic.*.py`, `rap212_*.py`). |

## How to Modify or Extend

* Scripts in `code/05_supplementary_figures/` read `File_1`, `File_2`, `File_3`, and `Tutorial_table3_icat.csv` directly. To test another set of target genes, donor promoters, or variants, replace or point CLI arguments to new `.txt`/`.vcf.gz` files formatted identically to `File_1`–`File_3`.
