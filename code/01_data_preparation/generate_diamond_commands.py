"""
Purpose:
    Generate a executable shell script running `diamond makedb` and self-`diamond blastp`
    on each species' Helixer proteome FASTA for paralog detection.

Manuscript (Peleke et al., 2026):
    Generates batch command scripts of the form committed in `run_my_diamond_jobs.sh`
    (DIAMOND `blastp --outfmt 6 --ultra-sensitive --max-hsps 1`).

Usage & I/O:
    python3 code/01_data_preparation/generate_diamond_commands.py \
        --proteome_dir <dir> --diamond_db_dir <db_dir> --diamond_out_dir <out_dir> \
        [--output_script_file run_all_diamond.sh] [--num_threads 4] [--max_target_seqs 25] \
        [--diamond_executable diamond] [--evalue 1e-3]

How to change or extend:
    Edit `species_configs` inside `main()` to add or update `(id_for_db, peptide_fasta_name,
    diamond_output_filename)` entries for additional species.
"""
import os
import argparse
import stat

def main():
    parser = argparse.ArgumentParser(description="Generate a shell script to run DIAMOND makedb and blastp for multiple species.")
    parser.add_argument("--proteome_dir", required=True, help="Directory containing input peptide FASTA files (or parent directory if FASTA names include subpaths).")
    parser.add_argument("--diamond_db_dir", required=True, help="Directory to store generated DIAMOND databases (.dmnd files).")
    parser.add_argument("--diamond_out_dir", required=True, help="Directory to store DIAMOND blastp output files.")
    parser.add_argument("--output_script_file", default="run_all_diamond.sh", help="Name of the shell script to generate (default: run_all_diamond.sh).")
    parser.add_argument("--num_threads", type=int, default=4, help="Number of threads for DIAMOND (makedb and blastp) (default: 4).")
    parser.add_argument("--max_target_seqs", type=int, default=25, help="Value for --max-target-seqs in diamond blastp (default: 25).")
    parser.add_argument("--diamond_executable", default="diamond", help="Path to the DIAMOND executable (default: 'diamond', assumes it's in PATH).")
    parser.add_argument("--evalue", type=str, default="1e-3", help="E-value threshold for diamond blastp (default: 1e-3).")

    args = parser.parse_args()

    # --- Species Configurations ---
    # CRITICAL: User must verify and update 'peptide_fasta_name' for each species.
    # 'peptide_fasta_name' should be the path relative to the '--proteome_dir' argument.
    # For example, if --proteome_dir is data/
    # and the file is in arabidopsis_halleri/peptide.faa, then peptide_fasta_name should be
    # "arabidopsis_halleri/peptide.faa".
    # 'diamond_output_filename' should match what 'produce_non_homologous_sets.py' expects.
    species_configs = [
        {
            "id_for_db": "ahal", 
            "peptide_fasta_name": "arabidopsis_halleri/hlx-Arabidopsis_halleri_subsp_gemmifera-GCA_900078215.1-63677.proteome.faa", # Example, verify
            "diamond_output_filename": "BLAST_ahal_to_ahal"
        },
        {
            "id_for_db": "atha", 
            "peptide_fasta_name": "arabidopsis_thaliana/GCF_000001735.4-genome_peptides_V3.fasta", # Example, verify
            "diamond_output_filename": "BLAST_atha_to_atha"
        },
        {
            "id_for_db": "bvul", 
            "peptide_fasta_name": "beta_vulgaris/GCF_002917755.1-genome_peptides_V3.fasta", # Example, verify
            "diamond_output_filename": "BLAST_bvul_to_bvul"
        },
        {
            "id_for_db": "bnap", 
            "peptide_fasta_name": "brassica_napus/Brassica_napus-GCF_020379485.1-proteome.helixer.faa", # Example, verify
            "diamond_output_filename": "BLAST_bnap_to_bnap"
        },
        {
          "id_for_db": "bole", 
          "peptide_fasta_name": "brassica_oleracea/GCF_000695525.1-genome_peptides_V3.fasta", # Example, verify
          "diamond_output_filename": "BLAST_bole_to_bole"
        },
        {
          "id_for_db": "csat", # Camelina sativa
          "peptide_fasta_name": "camelina_sativa/GCF_000633955.1-genome_peptides_V3.fasta", # Example, verify
          "diamond_output_filename": "BLAST_csat_to_csat" 
        },
        {
          "id_for_db": "cqui", 
          "peptide_fasta_name": "chenopodium_quinoa/chenopodium_quinoa_peptides_090724_V2.fasta", # Example, verify
          "diamond_output_filename": "BLAST_cqui_to_cqui"
        },
        {
          "id_for_db": "cucu", # Cucumis sativus
          "peptide_fasta_name": "cucumis_sativus/GCF_000004075.3-genome_peptides_V3.fasta", # Example, verify
          "diamond_output_filename": "BLAST_cucu_to_cucu"
        },
        {
          "id_for_db": "dcar", # Daucus carota
          "peptide_fasta_name": "daucus_carota/GCF_001625215.1-genome_peptides_V3.fasta", # Example, verify
          "diamond_output_filename": "BLAST_dcar_to_dcar"
        },
        { 
          "id_for_db": "gmax", 
          "peptide_fasta_name": "glycine_max/GCF_000004515.6-genome_peptides_V3.fasta", # Example, verify
          "diamond_output_filename": "BLAST_gmax_to_gmax" 
        },
        { 
          "id_for_db": "ghir", # Gossypium hirsutum
          "peptide_fasta_name": "gossypium_hirsutum/GCF_007990345.1-genome_peptides_V3.fasta", # Example, verify
          "diamond_output_filename": "BLAST_ghir_to_ghir" 
        },
        {
        "id_for_db": "hvul", # Hordeum vulgare
        "peptide_fasta_name": "hordeum_vulgare/hlx-Hordeum_vulgare_supsp_vulgare-GCF_904849725.1-112509_peptides_V3.fasta", # Example, verify
        "diamond_output_filename": "BLAST_hvul_to_hvul" 
        },
        {
        "id_for_db": "ntab", # Nicotiana tabacum
        "peptide_fasta_name": "nicotiana_tabacum/nicotiana_tabacum_peptides_090724.fasta", # Example, verify
        "diamond_output_filename": "BLAST_ntab_to_ntab" 
        },
        {
          "id_for_db": "osat", # Oryza sativa
          "peptide_fasta_name": "oryza_sativa/Oryza_sativa_Japonica_Group-GCF_001433935.1-proteome.helixer.faa", # Example, verify
          "diamond_output_filename": "BLAST_osat_to_osat" 
        },
        {
            "id_for_db": "slyc", 
            "peptide_fasta_name": "solanum_lycopersicum/Solanum_lycopersicum-GCF_000188115.5-proteome.helixer.faa", # Example, verify
            "diamond_output_filename": "BLAST_slyc_to_slyc" 
        },
        {
            "id_for_db": "stub", # Solanum tuberosum
            "peptide_fasta_name": "solanum_tuberosum/hlx-Solanum_tuberosum-GCF_000226075.1-4113.proteome.faa", # PLEASE REPLACE WITH ACTUAL FILENAME
            "diamond_output_filename": "BLAST_stub_to_stub"
        },
        {
            "id_for_db": "sbic", 
            "peptide_fasta_name": "sorghum_bicolor/Sorghum_bicolor-GCF_000003195.3-proteome.helixer.faa", # Example, verify
            "diamond_output_filename": "BLAST_sbic_to_sbic" 
      },
      {
            "id_for_db": "zmay", 
            "peptide_fasta_name": "zea_mays/Zea_mays-GCF_902167145.1-proteome.helixer.faa", # Example, verify
            "diamond_output_filename": "BLAST_zmay_to_zmay" 
        }
    ]
    # Filter out duplicate configurations (based on id_for_db for simplicity here, user should ensure uniqueness)
    unique_configs = []
    seen_ids = set()
    for config_item in species_configs: # Renamed config to config_item to avoid conflict with module
        if config_item["id_for_db"] not in seen_ids:
            unique_configs.append(config_item)
            seen_ids.add(config_item["id_for_db"])
    species_configs = unique_configs


    # Create output directories if they don't exist
    os.makedirs(args.diamond_db_dir, exist_ok=True)
    os.makedirs(args.diamond_out_dir, exist_ok=True)

    script_lines = [
        "#!/bin/bash",
        "# Shell script generated to run DIAMOND makedb and blastp",
        "set -e  # Exit immediately if a command exits with a non-zero status.",
        "set -u  # Treat unset variables as an error when substituting.",
        "set -o pipefail # Causes a pipeline to return the exit status of the last command in the pipe that failed.",
        "",
        f"NUM_THREADS={args.num_threads}",
        f"MAX_TARGET_SEQS={args.max_target_seqs}",
        f"DIAMOND_EXEC=\"{args.diamond_executable}\"", 
        f"EVALUE_THRESHOLD=\"{args.evalue}\"",
        "",
    ]

    for config_item in species_configs: # Renamed config to config_item
        species_id_db = config_item["id_for_db"]
        peptide_filename_cleaned = config_item["peptide_fasta_name"].strip()
        input_fasta_path = os.path.join(args.proteome_dir, peptide_filename_cleaned)
        
        diamond_result_filename = config_item["diamond_output_filename"].strip()

        db_path_diamond = os.path.join(args.diamond_db_dir, f"{species_id_db}.dmnd")
        diamond_output_path = os.path.join(args.diamond_out_dir, diamond_result_filename)

        script_lines.append(f"# --- Processing {species_id_db} ({peptide_filename_cleaned}) ---")
        script_lines.append(f"echo \"Processing {species_id_db}: Creating DIAMOND database from {input_fasta_path}...\"")
        
        makedb_cmd = (
            f"$DIAMOND_EXEC makedb --in \"{input_fasta_path}\" "
            f"-d \"{db_path_diamond}\" --threads $NUM_THREADS"
        )
        script_lines.append(makedb_cmd)
        
        script_lines.append(f"echo \"Processing {species_id_db}: Running DIAMOND blastp...\"")
        blastp_cmd = (
            f"$DIAMOND_EXEC blastp -d \"{db_path_diamond}\" "
            f"-q \"{input_fasta_path}\" "
            f"-o \"{diamond_output_path}\" "
            f"--outfmt 6 --threads $NUM_THREADS --ultra-sensitive --max-hsps 1 "
            f"--max-target-seqs $MAX_TARGET_SEQS --evalue $EVALUE_THRESHOLD"
        )
        script_lines.append(blastp_cmd)
        script_lines.append("echo \"Done.\"")
        script_lines.append("") 

    script_lines.append("echo \"All DIAMOND processes completed.\"")

    try:
        with open(args.output_script_file, "w") as f:
            for line in script_lines:
                f.write(line + "\n")
        
        current_permissions = os.stat(args.output_script_file).st_mode
        os.chmod(
            args.output_script_file,
            current_permissions | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
        )
        print(f"\nSuccessfully generated shell script: {args.output_script_file}")
        print("Please review the script and then run it from your terminal, for example:")
        print(f"  bash ./{args.output_script_file}")
        print(f"\nEnsure that DIAMOND (version supporting these options) is installed and accessible via '{args.diamond_executable}'.")
        print(f"\nCRITICAL: Ensure the '--proteome_dir' argument you provide when running this Python script is set to the parent directory of your species folders (e.g., data/).")
        print(f"AND that 'peptide_fasta_name' entries in 'species_configs' correctly include the species subdirectory and the peptide FASTA filename (e.g., 'arabidopsis_halleri/your_peptide_file.faa').")


    except Exception as e:
        print(f"Error generating shell script: {e}")

if __name__ == "__main__":
    main()
