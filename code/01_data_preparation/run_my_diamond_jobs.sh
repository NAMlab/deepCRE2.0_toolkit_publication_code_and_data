#!/bin/bash
# Purpose: Create DIAMOND databases and run sensitive all-vs-all self-`blastp` per species proteome.
# Manuscript (Peleke et al., 2026): Generated `EXTERNAL/diamond/BLASTdP_<sp>_to_<sp>` inputs for
#   paralog-filtered validation set selection (`prepare_non_homologous_validation_sets_dPHY.8.py`).
# Usage: Decompress proteomes in `data/<species>/` (or update filenames via `generate_diamond_commands.py`) and run from repo root: `bash code/01_data_preparation/run_my_diamond_jobs.sh`.
# How to change or extend: Regenerate or add species using `generate_diamond_commands.py`.
set -e # Exit immediately if a command exits with a non-zero status.
set -u # Treat unset variables as an error when substituting.
set -o pipefail # Causes a pipeline to return the exit status of the last command in the pipe that failed.

NUM_THREADS=4
MAX_TARGET_SEQS=25
DIAMOND_EXEC="${DIAMOND_EXEC:-diamond}"
EVALUE_THRESHOLD="1e-3"
mkdir -p EXTERNAL/diamond/databases EXTERNAL/diamond

# --- Processing atha (arabidopsis_thaliana/GCF_000001735.4-genome_peptides_V3.dPHY.faa) ---
echo "Processing atha: Creating DIAMOND database from data/arabidopsis_thaliana/GCF_000001735.4-genome_peptides_V3.dPHY.faa..."
$DIAMOND_EXEC makedb --in "data/arabidopsis_thaliana/GCF_000001735.4-genome_peptides_V3.dPHY.faa" -d "EXTERNAL/diamond/databases/atha.dmnd" --threads $NUM_THREADS
echo "Processing atha: Running DIAMOND blastp..."
$DIAMOND_EXEC blastp -d "EXTERNAL/diamond/databases/atha.dmnd" -q "data/arabidopsis_thaliana/GCF_000001735.4-genome_peptides_V3.dPHY.faa" -o "EXTERNAL/diamond/BLASTdP_atha_to_atha" --outfmt 6 --threads $NUM_THREADS --ultra-sensitive --max-hsps 1 --max-target-seqs $MAX_TARGET_SEQS --evalue $EVALUE_THRESHOLD
echo "Done."

# --- Processing ahal (data/arabidopsis_halleri/GCA_900078215.1-63677.proteome.dPHY.faa) ---
echo "Processing ahal: Creating DIAMOND database from data/arabidopsis_halleri/GCA_900078215.1-63677.proteome.dPHY.faa..."
$DIAMOND_EXEC makedb --in "data/arabidopsis_halleri/GCA_900078215.1-63677.proteome.dPHY.faa" -d "EXTERNAL/diamond/databases/ahal.dmnd" --threads $NUM_THREADS
echo "Processing ahal: Running DIAMOND blastp..."
$DIAMOND_EXEC blastp -d "EXTERNAL/diamond/databases/ahal.dmnd" -q "data/arabidopsis_halleri/GCA_900078215.1-63677.proteome.dPHY.faa" -o "EXTERNAL/diamond/BLASTdP_ahal_to_ahal" --outfmt 6 --threads $NUM_THREADS --ultra-sensitive --max-hsps 1 --max-target-seqs $MAX_TARGET_SEQS --evalue $EVALUE_THRESHOLD
echo "Done."

# --- Processing bvul (data/beta_vulgaris/GCF_002917755.1-genome_peptides_V3.dPHY.faa) ---
echo "Processing bvul: Creating DIAMOND database from data/beta_vulgaris/GCF_002917755.1-genome_peptides_V3.dPHY.faa..."
$DIAMOND_EXEC makedb --in "data/beta_vulgaris/GCF_002917755.1-genome_peptides_V3.dPHY.faa" -d "EXTERNAL/diamond/databases/bvul.dmnd" --threads $NUM_THREADS
echo "Processing bvul: Running DIAMOND blastp..."
$DIAMOND_EXEC blastp -d "EXTERNAL/diamond/databases/bvul.dmnd" -q "data/beta_vulgaris/GCF_002917755.1-genome_peptides_V3.dPHY.faa" -o "EXTERNAL/diamond/BLASTdP_bvul_to_bvul" --outfmt 6 --threads $NUM_THREADS --ultra-sensitive --max-hsps 1 --max-target-seqs $MAX_TARGET_SEQS --evalue $EVALUE_THRESHOLD
echo "Done."

# --- Processing bnap (data/brassica_napus/Brassica_napus-GCF_020379485.1-proteome.helixer.dPHY.faa) ---
echo "Processing bnap: Creating DIAMOND database from data/brassica_napus/Brassica_napus-GCF_020379485.1-proteome.helixer.dPHY.faa..."
$DIAMOND_EXEC makedb --in "data/brassica_napus/Brassica_napus-GCF_020379485.1-proteome.helixer.dPHY.faa" -d "EXTERNAL/diamond/databases/bnap.dmnd" --threads $NUM_THREADS
echo "Processing bnap: Running DIAMOND blastp..."
$DIAMOND_EXEC blastp -d "EXTERNAL/diamond/databases/bnap.dmnd" -q "data/brassica_napus/Brassica_napus-GCF_020379485.1-proteome.helixer.dPHY.faa" -o "EXTERNAL/diamond/BLASTdP_bnap_to_bnap" --outfmt 6 --threads $NUM_THREADS --ultra-sensitive --max-hsps 1 --max-target-seqs $MAX_TARGET_SEQS --evalue $EVALUE_THRESHOLD
echo "Done."

# --- Processing bole (data/brassica_oleracea/GCF_000695525.1-genome_peptides_V3.dPHY.faa) ---
echo "Processing bole: Creating DIAMOND database from data/brassica_oleracea/GCF_000695525.1-genome_peptides_V3.dPHY.faa..."
$DIAMOND_EXEC makedb --in "data/brassica_oleracea/GCF_000695525.1-genome_peptides_V3.dPHY.faa" -d "EXTERNAL/diamond/databases/bole.dmnd" --threads $NUM_THREADS
echo "Processing bole: Running DIAMOND blastp..."
$DIAMOND_EXEC blastp -d "EXTERNAL/diamond/databases/bole.dmnd" -q "data/brassica_oleracea/GCF_000695525.1-genome_peptides_V3.dPHY.faa" -o "EXTERNAL/diamond/BLASTdP_bole_to_bole" --outfmt 6 --threads $NUM_THREADS --ultra-sensitive --max-hsps 1 --max-target-seqs $MAX_TARGET_SEQS --evalue $EVALUE_THRESHOLD
echo "Done."

# --- Processing csat (data/camelina_sativa/GCF_000633955.1-genome_peptides_V3.dPHY.faa) ---
echo "Processing csat: Creating DIAMOND database from data/camelina_sativa/GCF_000633955.1-genome_peptides_V3.dPHY.faa..."
$DIAMOND_EXEC makedb --in "data/camelina_sativa/GCF_000633955.1-genome_peptides_V3.dPHY.faa" -d "EXTERNAL/diamond/databases/csat.dmnd" --threads $NUM_THREADS
echo "Processing csat: Running DIAMOND blastp..."
$DIAMOND_EXEC blastp -d "EXTERNAL/diamond/databases/csat.dmnd" -q "data/camelina_sativa/GCF_000633955.1-genome_peptides_V3.dPHY.faa" -o "EXTERNAL/diamond/BLASTdP_csat_to_csat" --outfmt 6 --threads $NUM_THREADS --ultra-sensitive --max-hsps 1 --max-target-seqs $MAX_TARGET_SEQS --evalue $EVALUE_THRESHOLD
echo "Done."

# --- Processing cqui (data/chenopodium_quinoa/chenopodium_quinoa_peptides_090724_V2.dPHY.faa) ---
echo "Processing cqui: Creating DIAMOND database from data/chenopodium_quinoa/chenopodium_quinoa_peptides_090724_V2.dPHY.faa..."
$DIAMOND_EXEC makedb --in "data/chenopodium_quinoa/chenopodium_quinoa_peptides_090724_V2.dPHY.faa" -d "EXTERNAL/diamond/databases/cqui.dmnd" --threads $NUM_THREADS
echo "Processing cqui: Running DIAMOND blastp..."
$DIAMOND_EXEC blastp -d "EXTERNAL/diamond/databases/cqui.dmnd" -q "data/chenopodium_quinoa/chenopodium_quinoa_peptides_090724_V2.dPHY.faa" -o "EXTERNAL/diamond/BLASTdP_cqui_to_cqui" --outfmt 6 --threads $NUM_THREADS --ultra-sensitive --max-hsps 1 --max-target-seqs $MAX_TARGET_SEQS --evalue $EVALUE_THRESHOLD
echo "Done."

# --- Processing cucu (data/cucumis_sativus/GCF_000004075.3-genome_peptides_V3.dPHY.faa) ---
echo "Processing cucu: Creating DIAMOND database from data/cucumis_sativus/GCF_000004075.3-genome_peptides_V3.dPHY.faa..."
$DIAMOND_EXEC makedb --in "data/cucumis_sativus/GCF_000004075.3-genome_peptides_V3.dPHY.faa" -d "EXTERNAL/diamond/databases/cucu.dmnd" --threads $NUM_THREADS
echo "Processing cucu: Running DIAMOND blastp..."
$DIAMOND_EXEC blastp -d "EXTERNAL/diamond/databases/cucu.dmnd" -q "data/cucumis_sativus/GCF_000004075.3-genome_peptides_V3.dPHY.faa" -o "EXTERNAL/diamond/BLASTdP_cucu_to_cucu" --outfmt 6 --threads $NUM_THREADS --ultra-sensitive --max-hsps 1 --max-target-seqs $MAX_TARGET_SEQS --evalue $EVALUE_THRESHOLD
echo "Done."

# --- Processing dcar (data/daucus_carota/GCF_001625215.1-genome_peptides_V3.dPHY.faa) ---
echo "Processing dcar: Creating DIAMOND database from data/daucus_carota/GCF_001625215.1-genome_peptides_V3.dPHY.faa..."
$DIAMOND_EXEC makedb --in "data/daucus_carota/GCF_001625215.1-genome_peptides_V3.dPHY.faa" -d "EXTERNAL/diamond/databases/dcar.dmnd" --threads $NUM_THREADS
echo "Processing dcar: Running DIAMOND blastp..."
$DIAMOND_EXEC blastp -d "EXTERNAL/diamond/databases/dcar.dmnd" -q "data/daucus_carota/GCF_001625215.1-genome_peptides_V3.dPHY.faa" -o "EXTERNAL/diamond/BLASTdP_dcar_to_dcar" --outfmt 6 --threads $NUM_THREADS --ultra-sensitive --max-hsps 1 --max-target-seqs $MAX_TARGET_SEQS --evalue $EVALUE_THRESHOLD
echo "Done."

# --- Processing gmax (data/glycine_max/GCF_000004515.6-genome_peptides_V3.dPHY.faa) ---
echo "Processing gmax: Creating DIAMOND database from data/glycine_max/GCF_000004515.6-genome_peptides_V3.dPHY.faa..."
$DIAMOND_EXEC makedb --in "data/glycine_max/GCF_000004515.6-genome_peptides_V3.dPHY.faa" -d "EXTERNAL/diamond/databases/gmax.dmnd" --threads $NUM_THREADS
echo "Processing gmax: Running DIAMOND blastp..."
$DIAMOND_EXEC blastp -d "EXTERNAL/diamond/databases/gmax.dmnd" -q "data/glycine_max/GCF_000004515.6-genome_peptides_V3.dPHY.faa" -o "EXTERNAL/diamond/BLASTdP_gmax_to_gmax" --outfmt 6 --threads $NUM_THREADS --ultra-sensitive --max-hsps 1 --max-target-seqs $MAX_TARGET_SEQS --evalue $EVALUE_THRESHOLD
echo "Done."

# --- Processing ghir (data/gossypium_hirsutum/GCF_007990345.1-genome_peptides_V3.dPHY.faa) ---
echo "Processing ghir: Creating DIAMOND database from data/gossypium_hirsutum/GCF_007990345.1-genome_peptides_V3.dPHY.faa..."
$DIAMOND_EXEC makedb --in "data/gossypium_hirsutum/GCF_007990345.1-genome_peptides_V3.dPHY.faa" -d "EXTERNAL/diamond/databases/ghir.dmnd" --threads $NUM_THREADS
echo "Processing ghir: Running DIAMOND blastp..."
$DIAMOND_EXEC blastp -d "EXTERNAL/diamond/databases/ghir.dmnd" -q "data/gossypium_hirsutum/GCF_007990345.1-genome_peptides_V3.dPHY.faa" -o "EXTERNAL/diamond/BLASTdP_ghir_to_ghir" --outfmt 6 --threads $NUM_THREADS --ultra-sensitive --max-hsps 1 --max-target-seqs $MAX_TARGET_SEQS --evalue $EVALUE_THRESHOLD
echo "Done."

# --- Processing hvul (data/hordeum_vulgare/hlx-Hordeum_vulgare_supsp_vulgare-GCF_904849725.1-112509_peptides_V3.dPHY.faa) ---
echo "Processing hvul: Creating DIAMOND database from data/hordeum_vulgare/hlx-Hordeum_vulgare_supsp_vulgare-GCF_904849725.1-112509_peptides_V3.dPHY.faa..."
$DIAMOND_EXEC makedb --in "data/hordeum_vulgare/hlx-Hordeum_vulgare_supsp_vulgare-GCF_904849725.1-112509_peptides_V3.dPHY.faa" -d "EXTERNAL/diamond/databases/hvul.dmnd" --threads $NUM_THREADS
echo "Processing hvul: Running DIAMOND blastp..."
$DIAMOND_EXEC blastp -d "EXTERNAL/diamond/databases/hvul.dmnd" -q "data/hordeum_vulgare/hlx-Hordeum_vulgare_supsp_vulgare-GCF_904849725.1-112509_peptides_V3.dPHY.faa" -o "EXTERNAL/diamond/BLASTdP_hvul_to_hvul" --outfmt 6 --threads $NUM_THREADS --ultra-sensitive --max-hsps 1 --max-target-seqs $MAX_TARGET_SEQS --evalue $EVALUE_THRESHOLD
echo "Done."

# --- Processing ntab (data/nicotiana_tabacum/nicotiana_tabacum_peptides_090724.dPHY.faa) ---
echo "Processing ntab: Creating DIAMOND database from data/nicotiana_tabacum/nicotiana_tabacum_peptides_090724.dPHY.faa..."
$DIAMOND_EXEC makedb --in "data/nicotiana_tabacum/nicotiana_tabacum_peptides_090724.dPHY.faa" -d "EXTERNAL/diamond/databases/ntab.dmnd" --threads $NUM_THREADS
echo "Processing ntab: Running DIAMOND blastp..."
$DIAMOND_EXEC blastp -d "EXTERNAL/diamond/databases/ntab.dmnd" -q "data/nicotiana_tabacum/nicotiana_tabacum_peptides_090724.dPHY.faa" -o "EXTERNAL/diamond/BLASTdP_ntab_to_ntab" --outfmt 6 --threads $NUM_THREADS --ultra-sensitive --max-hsps 1 --max-target-seqs $MAX_TARGET_SEQS --evalue $EVALUE_THRESHOLD
echo "Done."

# --- Processing osat (data/oryza_sativa/Oryza_sativa_Japonica_Group-GCF_001433935.1-proteome.helixer.dPHY.faa) ---
echo "Processing osat: Creating DIAMOND database from data/oryza_sativa/Oryza_sativa_Japonica_Group-GCF_001433935.1-proteome.helixer.dPHY.faa..."
$DIAMOND_EXEC makedb --in "data/oryza_sativa/Oryza_sativa_Japonica_Group-GCF_001433935.1-proteome.helixer.dPHY.faa" -d "EXTERNAL/diamond/databases/osat.dmnd" --threads $NUM_THREADS
echo "Processing osat: Running DIAMOND blastp..."
$DIAMOND_EXEC blastp -d "EXTERNAL/diamond/databases/osat.dmnd" -q "data/oryza_sativa/Oryza_sativa_Japonica_Group-GCF_001433935.1-proteome.helixer.dPHY.faa" -o "EXTERNAL/diamond/BLASTdP_osat_to_osat" --outfmt 6 --threads $NUM_THREADS --ultra-sensitive --max-hsps 1 --max-target-seqs $MAX_TARGET_SEQS --evalue $EVALUE_THRESHOLD
echo "Done."

# --- Processing slyc (data/solanum_lycopersicum/Solanum_lycopersicum-GCF_000188115.5-proteome.helixer.dPHY.faa) ---
echo "Processing slyc: Creating DIAMOND database from data/solanum_lycopersicum/Solanum_lycopersicum-GCF_000188115.5-proteome.helixer.dPHY.faa..."
$DIAMOND_EXEC makedb --in "data/solanum_lycopersicum/Solanum_lycopersicum-GCF_000188115.5-proteome.helixer.dPHY.faa" -d "EXTERNAL/diamond/databases/slyc.dmnd" --threads $NUM_THREADS
echo "Processing slyc: Running DIAMOND blastp..."
$DIAMOND_EXEC blastp -d "EXTERNAL/diamond/databases/slyc.dmnd" -q "data/solanum_lycopersicum/Solanum_lycopersicum-GCF_000188115.5-proteome.helixer.dPHY.faa" -o "EXTERNAL/diamond/BLASTdP_slyc_to_slyc" --outfmt 6 --threads $NUM_THREADS --ultra-sensitive --max-hsps 1 --max-target-seqs $MAX_TARGET_SEQS --evalue $EVALUE_THRESHOLD
echo "Done."

# --- Processing stub (data/solanum_lycopersicum/Solanum_lycopersicum-GCF_000188115.5-proteome.helixer.dPHY.faa) ---
echo "Processing stub: Creating DIAMOND database from data/solanum_lycopersicum/Solanum_lycopersicum-GCF_000188115.5-proteome.helixer.dPHY.faa..."
$DIAMOND_EXEC makedb --in "data/solanum_lycopersicum/Solanum_lycopersicum-GCF_000188115.5-proteome.helixer.dPHY.faa" -d "EXTERNAL/diamond/databases/stub.dmnd" --threads $NUM_THREADS
echo "Processing stub: Running DIAMOND blastp..."
$DIAMOND_EXEC blastp -d "EXTERNAL/diamond/databases/stub.dmnd" -q "data/solanum_lycopersicum/Solanum_lycopersicum-GCF_000188115.5-proteome.helixer.dPHY.faa" -o "EXTERNAL/diamond/BLASTdP_stub_to_stub" --outfmt 6 --threads $NUM_THREADS --ultra-sensitive --max-hsps 1 --max-target-seqs $MAX_TARGET_SEQS --evalue $EVALUE_THRESHOLD
echo "Done."

# --- Processing sbic (data/sorghum_bicolor/Sorghum_bicolor-GCF_000003195.3-proteome.helixer.dPHY.faa) ---
echo "Processing sbic: Creating DIAMOND database from data/sorghum_bicolor/Sorghum_bicolor-GCF_000003195.3-proteome.helixer.dPHY.faa..."
$DIAMOND_EXEC makedb --in "data/sorghum_bicolor/Sorghum_bicolor-GCF_000003195.3-proteome.helixer.dPHY.faa" -d "EXTERNAL/diamond/databases/sbic.dmnd" --threads $NUM_THREADS
echo "Processing sbic: Running DIAMOND blastp..."
$DIAMOND_EXEC blastp -d "EXTERNAL/diamond/databases/sbic.dmnd" -q "data/sorghum_bicolor/Sorghum_bicolor-GCF_000003195.3-proteome.helixer.dPHY.faa" -o "EXTERNAL/diamond/BLASTdP_sbic_to_sbic" --outfmt 6 --threads $NUM_THREADS --ultra-sensitive --max-hsps 1 --max-target-seqs $MAX_TARGET_SEQS --evalue $EVALUE_THRESHOLD
echo "Done."

# --- Processing zmay (data/zea_mays/Zea_mays-GCF_902167145.1-proteome.helixer.dPHY.faa) ---
echo "Processing zmay: Creating DIAMOND database from data/zea_mays/Zea_mays-GCF_902167145.1-proteome.helixer.dPHY.faa..."
$DIAMOND_EXEC makedb --in "data/zea_mays/Zea_mays-GCF_902167145.1-proteome.helixer.dPHY.faa" -d "EXTERNAL/diamond/databases/zmay.dmnd" --threads $NUM_THREADS
echo "Processing zmay: Running DIAMOND blastp..."
$DIAMOND_EXEC blastp -d "EXTERNAL/diamond/databases/zmay.dmnd" -q "data/zea_mays/Zea_mays-GCF_902167145.1-proteome.helixer.dPHY.faa" -o "EXTERNAL/diamond/BLASTdP_zmay_to_zmay" --outfmt 6 --threads $NUM_THREADS --ultra-sensitive --max-hsps 1 --max-target-seqs $MAX_TARGET_SEQS --evalue $EVALUE_THRESHOLD
echo "Done."

echo "All DIAMOND processes completed."
