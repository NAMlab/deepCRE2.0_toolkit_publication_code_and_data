"""
Purpose:
    Earlier validation-gene selector excluding paralogs using all non-self DIAMOND hits
    (`E <= 1e-3`) with tier-3 fallback filling up to the target per-species count.

Manuscript (Peleke et al., 2026):
    Legacy `v7` validation builder preserved for exact comparison against `v8`
    (`prepare_non_homologous_validation_sets_dPHY.8.py`); not used for final 2026 models.

Usage & I/O:
    python code/01_data_preparation/prepare_non_homologous_validation_sets_dPHY.7.py \
        --species_config <config.json> --output_file <out.pickle> [options]

How to change or extend:
    Use `prepare_non_homologous_validation_sets_dPHY.8.py` for final workflows (`--allow_tier3`
    on `.8.py` reproduces the tier-3 fallback behavior when needed).
"""
import os
import pickle
import pandas as pd
from Bio import SeqIO
import argparse
import re
import ast # For safely evaluating string to Python literal
import random # For sampling if needed

# --- Configuration ---
EXCLUDED_CHROMS_GENERAL = ['Pt', 'Mt', 'pt', 'mt', 'PT', 'MT', 'CP',
                           'PLASTID', 'MITOCHONDRION',
                           'chloroplast', 'mitochondrion',
                           'cp'] #inlcude 'scaffold' o.s. if excluded

def parse_fasta_header(protein_id_str, description_line, species_id_for_parsing=None):
    """
    Extracts gene_id and chromosome from FASTA header information.
    Prioritizes 'gene=' and 'seq_id=' tags in the description.
    Then tries to use the protein_id_str itself (minus version) as gene_id.
    """
    gene_id = None
    chromosome = None

    # --- Step 1: Attempt to extract gene_id ---
    # Priority A: Look for 'gene=...' in the description.
    gene_match_desc = re.search(r'gene=([^\s;]+)', description_line)
    if gene_match_desc:
        gene_id = gene_match_desc.group(1)
    
    # Fallback if 'gene=' not found: Use the protein_id_str from the > line,
    # stripping common protein version suffixes.
    if not gene_id and protein_id_str:
        print(f" gene_id not found for '{protein_id_str}'. Attempting to parse from protein ID...")

        temp_gene_id = str(protein_id_str) 
        temp_gene_id = re.sub(r'\.[pt]\d+$', '', temp_gene_id, flags=re.IGNORECASE)
        temp_gene_id = re.sub(r'\.\d+$', '', temp_gene_id) 
        gene_id = temp_gene_id
    
    # --- Step 2: Attempt to extract chromosome ---
    seq_id_match_desc = re.search(r'seq_id=([^\s;]+)', description_line) #this has been tested on the given helixer annotations only
    if seq_id_match_desc:
        chromosome = seq_id_match_desc.group(1)
    
    if not chromosome:
        chrom_match = re.search(r'(?:chromosome|chrom|chr)[:=]\s*([^\s;,]+)', description_line, re.IGNORECASE)
        if chrom_match:
            chromosome = chrom_match.group(1)
        else:
            chrom_match_gb_desc = re.search(r'\[chromosome=([^\]]+)\]', description_line, re.IGNORECASE)
            if chrom_match_gb_desc:
                chromosome = chrom_match_gb_desc.group(1)
            else:
                general_chrom_match_desc = re.search(r'\b(?:Chr|chromosome|scaffold)_?([a-zA-Z0-9_.-]+)\b', description_line, re.IGNORECASE)
                if general_chrom_match_desc:
                    chromosome = general_chrom_match_desc.group(1)

    # --- Step 3: Further fallback for chromosome if still not found, using protein_id_str ---
    if not chromosome and protein_id_str:
        print(f" Chromosome not found for '{protein_id_str}'. Attempting to parse from protein ID...")

        base_id_for_chrom_parse = re.sub(r'\.\d+$', '', str(protein_id_str))
        base_id_for_chrom_parse = re.sub(r'\.[pt]\d+$', '', base_id_for_chrom_parse, flags=re.IGNORECASE)

        acc_match = re.search(r'([A-Z]{2,4}\d{6,}\.\d+)', base_id_for_chrom_parse, re.IGNORECASE) 
        if acc_match:
            chromosome = acc_match.group(1)
        else:
            if '_' in base_id_for_chrom_parse:
                potential_chrom_candidate = base_id_for_chrom_parse.rsplit('_', 1)[0]
                if re.search(r'\.\d+$', potential_chrom_candidate) or \
                   potential_chrom_candidate.lower().startswith(('nc_', 'nw_', 'nt_', 'scaffold_', 'chr', 'chromosome')):
                    chromosome = potential_chrom_candidate

   # Inconsistent labelling and annotation of .fna, .faa, .gtf or .gff makes these fallbacks helpfull (i.e. Arabis cases)                  
    # --- Step 4: Final Cleaning ---
    if gene_id:
        gene_id = str(gene_id).strip()
        gene_id = re.sub(r'\.[pt]\d+$', '', gene_id, flags=re.IGNORECASE)
        gene_id = re.sub(r'\.\d+$', '', gene_id)

    if chromosome:
        chromosome = str(chromosome).strip()
    return gene_id, chromosome

def parse_proteome_fasta(proteome_fasta_path, species_id):
    if not os.path.exists(proteome_fasta_path):
        print(f"Error: Proteome FASTA file not found: {proteome_fasta_path}")
        return None
    parsed_records = []
    print(f"  Parsing proteome: {os.path.basename(proteome_fasta_path)} for species: {species_id}...")
    try:
        for rec in SeqIO.parse(proteome_fasta_path, "fasta"):
            protein_id = rec.id
            gene_id, chrom = parse_fasta_header(protein_id, rec.description, species_id_for_parsing=species_id)
            
            if protein_id and gene_id and chrom:
                parsed_records.append([protein_id, gene_id, chrom])
            else:
                print(f"    Warning: Could not parse full gene/chrom for protein '{protein_id}' from header: '{rec.description[:100]}...' (Parsed ProteinID: {protein_id}, Parsed Gene: {gene_id}, Parsed Chrom: {chrom})")
        
        if not parsed_records:
            print(f"  Error: No records successfully parsed with gene_id and chromosome from {proteome_fasta_path}.")
            return None
            
        info_df = pd.DataFrame(parsed_records, columns=['protein_id', 'gene_id', 'chrom'])
        
        info_df['protein_id'] = info_df['protein_id'].astype(str)
        info_df['gene_id'] = info_df['gene_id'].astype(str)
        info_df['chrom'] = info_df['chrom'].astype(str)

        if info_df['protein_id'].duplicated().any():
            print(f"    Warning: Duplicate protein_ids found in {proteome_fasta_path}. Keeping first occurrence.")
            info_df.drop_duplicates(subset=['protein_id'], keep='first', inplace=True)
        
        info_df.set_index('protein_id', inplace=True)
        print(f"    Successfully parsed {len(info_df)} proteins with gene and chromosome info.")
        return info_df
    except Exception as e:
        print(f"  Error parsing FASTA file {proteome_fasta_path}: {e}")
        return None

def get_info_from_df(protein_id, column_name, df_to_lookup):
    if protein_id is None: return None
    try:
        return df_to_lookup.loc[str(protein_id), column_name]
    except KeyError:
        return None

def process_species_for_validation_set(species_id, peptide_fasta_path, blast_output_path, target_n_genes):
    print(f"\nProcessing species: {species_id} (Target N: {target_n_genes})")

    protein_info_df = parse_proteome_fasta(peptide_fasta_path, species_id)
    if protein_info_df is None or protein_info_df.empty:
        print(f"  Skipping {species_id} due to proteome parsing issues.")
        return []

    excluded_chroms_lower = [ch.lower() for ch in EXCLUDED_CHROMS_GENERAL]
    gene_to_main_chrom_map = {}
    
    protein_info_df_for_grouping = protein_info_df.copy()
    protein_info_df_for_grouping['gene_id'] = protein_info_df_for_grouping['gene_id'].astype(str)
    protein_info_df_for_grouping['chrom'] = protein_info_df_for_grouping['chrom'].astype(str)

    gene_chrom_groups = protein_info_df_for_grouping.groupby('gene_id')['chrom'].apply(lambda x: list(pd.unique(x))).to_dict()

    for gene, chroms in gene_chrom_groups.items():
        valid_chroms = [c for c in chroms if str(c).lower() not in excluded_chroms_lower and pd.notna(c) and str(c).strip() != ""]
        if len(valid_chroms) == 1:
            gene_to_main_chrom_map[str(gene)] = str(valid_chroms[0])
        elif len(valid_chroms) > 1:
            gene_to_main_chrom_map[str(gene)] = str(valid_chroms[0]) 
    
    all_candidate_genes = set(gene_to_main_chrom_map.keys())
    if not all_candidate_genes:
        print(f"  Error: No genes with valid/mapped chromosome information found for {species_id}.")
        return []

    final_validation_set = []
    processed_genes = set() 

    blast_file_usable = False
    non_self_hits_df = pd.DataFrame() 

    if os.path.exists(blast_output_path) and os.path.getsize(blast_output_path) > 0:
        print(f"  Reading DIAMOND/BLAST output: {os.path.basename(blast_output_path)}...")
        try:
            blast_df = pd.read_csv(blast_output_path, sep='\t',
                                   names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                                          'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                                          'evalue', 'bitscore'])
            blast_df['qseqid'] = blast_df['qseqid'].astype(str)
            blast_df['sseqid'] = blast_df['sseqid'].astype(str)

            blast_df['qgene'] = blast_df['qseqid'].map(protein_info_df['gene_id'])
            blast_df['sgene'] = blast_df['sseqid'].map(protein_info_df['gene_id'])            

            #blast_df['qgene'] = blast_df['qseqid'].apply(lambda x: get_info_from_df(x, 'gene_id', protein_info_df))
            #blast_df['sgene'] = blast_df['sseqid'].apply(lambda x: get_info_from_df(x, 'gene_id', protein_info_df))
            
            blast_df.dropna(subset=['qgene', 'sgene'], inplace=True)
            blast_df['qgene'] = blast_df['qgene'].astype(str)
            blast_df['sgene'] = blast_df['sgene'].astype(str)
            blast_df = blast_df[ (blast_df['qgene'] != 'None') & (blast_df['sgene'] != 'None') & \
                                 (blast_df['qgene'] != '') & (blast_df['sgene'] != '')]

            if not blast_df.empty:
                temp_non_self_hits_df = blast_df[blast_df['qgene'] != blast_df['sgene']].copy()
                
                qgenes_in_map = temp_non_self_hits_df['qgene'].isin(gene_to_main_chrom_map.keys())
                sgenes_in_map = temp_non_self_hits_df['sgene'].isin(gene_to_main_chrom_map.keys())

                non_self_hits_df = temp_non_self_hits_df[qgenes_in_map & sgenes_in_map]
                
                if not non_self_hits_df.empty:
                    blast_file_usable = True
                else:
                    print("  No valid non-self hits found after filtering for genes with mapped chromosomes.")
            else:
                print("  BLAST data empty after initial gene mapping or all genes failed mapping.")
        except Exception as e:
            print(f"  Error reading or processing DIAMOND/BLAST output file {blast_output_path}: {e}")
    else:
        print(f"  BLAST output file not found or empty: {blast_output_path}")

    genes_with_valid_non_self_hits = set()
    if blast_file_usable:
        genes_with_valid_non_self_hits = set(non_self_hits_df['qgene'].unique())
    
    tier1_candidates_list = list(all_candidate_genes - genes_with_valid_non_self_hits)
    print(f"  Identified {len(tier1_candidates_list)} Tier 1 candidates with no homologs in genome (single copy genes).")

    tier2_candidates_list = []
    if blast_file_usable:
        for qgene in (genes_with_valid_non_self_hits & all_candidate_genes):
            qgene_chrom = gene_to_main_chrom_map.get(qgene)
            qgene_hits = non_self_hits_df[non_self_hits_df['qgene'] == qgene]
            
            is_strictly_local = True
            if qgene_hits.empty: 
                is_strictly_local = False 
            
            for _, hit_row in qgene_hits.iterrows():
                sgene = str(hit_row['sgene'])
                sgene_chrom = gene_to_main_chrom_map.get(sgene)
                
                if sgene_chrom != qgene_chrom:
                    is_strictly_local = False
                    break
            
            if is_strictly_local:
                tier2_candidates_list.append(qgene)
    print(f"  Identified {len(tier2_candidates_list)} Tier 2 candidates with local homologs on chromosome.")

    combined_tier1_tier2_pool = tier1_candidates_list + tier2_candidates_list
    random.shuffle(combined_tier1_tier2_pool)
    
    print(f"  Combined pool of Tier 1 and Tier 2: {len(combined_tier1_tier2_pool)} unique candidates.")

    for gene in combined_tier1_tier2_pool:
        if len(final_validation_set) < target_n_genes:
            if gene not in processed_genes:
                final_validation_set.append(gene)
                processed_genes.add(gene)
        else:
            break
    print(f"  After selecting from mixed Tier 1/2 pool: {len(final_validation_set)}/{target_n_genes} genes selected.")

    if len(final_validation_set) < target_n_genes and blast_file_usable:
        tier3_candidate_data = []
        potential_tier3_qgenes = (genes_with_valid_non_self_hits - processed_genes) & all_candidate_genes

        for qgene in potential_tier3_qgenes:
            qgene_hits = non_self_hits_df[non_self_hits_df['qgene'] == qgene]
            if qgene_hits.empty: continue 

            max_bitscore_overall = qgene_hits['bitscore'].max()
            tier3_candidate_data.append((qgene, max_bitscore_overall))
        
        print(f"  Tier 3: Evaluating {len(tier3_candidate_data)} raw candidates from remaining genes with lowest pairwise similarities.")
        tier3_candidate_data.sort(key=lambda x: x[1])

        for gene, _ in tier3_candidate_data:
            if len(final_validation_set) < target_n_genes:
                if gene not in processed_genes: 
                    final_validation_set.append(gene)
                    processed_genes.add(gene)
            else:
                break
        print(f"  After Tier 3: {len(final_validation_set)}/{target_n_genes} genes selected.")

#    if len(final_validation_set) < target_n_genes:
#        remaining_candidates_for_random_fill = list(all_candidate_genes - processed_genes)
#        if remaining_candidates_for_random_fill:
#            print(f"  Tier 4: Adding up to {target_n_genes - len(final_validation_set)} remaining genes randomly.")
#            random.shuffle(remaining_candidates_for_random_fill)
#            needed = target_n_genes - len(final_validation_set)
#            
#            added_count = 0
#            for gene_to_add in remaining_candidates_for_random_fill:
#                if added_count >= needed:
#                    break
#                if gene_to_add not in processed_genes:
#                    final_validation_set.append(gene_to_add)
#                    processed_genes.add(gene_to_add) 
#                    added_count +=1
#        print(f"  After Tier 4 (random fill): {len(final_validation_set)}/{target_n_genes} genes selected.")
    
    print(f"  Found {len(final_validation_set)} validation genes for {species_id} (target was {target_n_genes}).")
    return list(set(final_validation_set))

def load_species_configs_from_file(config_filepath):
    if not os.path.exists(config_filepath):
        print(f"Error: Species configuration file not found: {config_filepath}")
        return None
    try:
        with open(config_filepath, 'r') as f:
            content = f.read()
        
        # Regex improved to be more robust for typical config list structures
        match = re.search(r'^\s*(?:[a-zA-Z_][a-zA-Z0-9_]*\s*=\s*)?(\[.*\])\s*$', content, re.DOTALL | re.MULTILINE)
        
        list_candidate_str = None
        if match:
            list_candidate_str = match.group(1) if match.lastindex == 1 else match.group(0)
        
        if not list_candidate_str: # Fallback if specific regex fails, try to find first complete list
            start_idx = content.find('[')
            end_idx = content.rfind(']')
            if not (start_idx != -1 and end_idx != -1 and end_idx > start_idx):
                print(f"Error: No valid list structure found in {config_filepath}. Content: {content[:200]}...")
                return None
            list_candidate_str = content[start_idx : end_idx+1]
        
        lines_in_list_str = []
        for line in list_candidate_str.splitlines():
            cleaned_line = line.split('#', 1)[0].rstrip()
            lines_in_list_str.append(cleaned_line)
        
        final_list_str_for_eval = "".join(lines_in_list_str)
        final_list_str_for_eval = re.sub(r',\s*(?=\])', '', final_list_str_for_eval) 
        final_list_str_for_eval = re.sub(r',\s*(?=\})', '', final_list_str_for_eval)

        species_configs = ast.literal_eval(final_list_str_for_eval)
        
        if not isinstance(species_configs, list):
            print(f"Error: Parsed configuration from {config_filepath} is not a list.")
            return None
        for item in species_configs:
            if not isinstance(item, dict) or not all(k in item for k in ["id_for_db", "out_peptide_fasta_name", "diamond_output_filename"]):
                print(f"Error: Configuration item is not a valid dictionary or missing required keys: {item}")
                return None
        print(f"Successfully loaded {len(species_configs)} species configurations from {config_filepath}")
        return species_configs
    except Exception as e:
        print(f"Error parsing species configuration file {config_filepath}: {e}")
        # print(f"DEBUG: Tried to parse: '''{final_list_str_for_eval if 'final_list_str_for_eval' in locals() else list_candidate_str if 'list_candidate_str' in locals() else 'N/A'}'''")
        return None

def main():
    parser = argparse.ArgumentParser(description="Produce sets of N 'most non-homologous' genes for validation.")
    parser.add_argument("--species_config_file", required=True, help="Path to the text file containing species configurations.")
    parser.add_argument("--proteome_dir", required=True, help="Base directory for peptide FASTA files.")
    parser.add_argument("--blast_dir", required=True, help="Directory containing DIAMOND/BLASTp output files.")
    parser.add_argument("--output_file", required=True, help="Path to save the output pickle file.")
    parser.add_argument("--num_validation_genes", type=int, default=2400, help="Target number of validation genes per species (default: 2400).")
    parser.add_argument("--target_species_id", type=str, default=None, 
                        help="Optional: Process only this specific species ID from the config file (e.g., 'bnap').")
    
    args = parser.parse_args()

    species_configs_all = load_species_configs_from_file(args.species_config_file)
    if species_configs_all is None:
        print("Exiting due to configuration loading error.")
        return

    species_configs_to_process = []
    if args.target_species_id:
        print(f"Targeting species ID: {args.target_species_id}")
        for config_item in species_configs_all:
            if config_item["id_for_db"] == args.target_species_id:
                species_configs_to_process.append(config_item)
                break 
        if not species_configs_to_process:
            print(f"Error: Target species ID '{args.target_species_id}' not found in the species configuration file.")
            print("Available species IDs in config:")
            for item in species_configs_all:
                print(f"- {item['id_for_db']}")
            return
    else:
        print("Processing all species found in the configuration file.")
        species_configs_to_process = species_configs_all

    all_validation_genes = {}

    if not species_configs_to_process:
        print("No species selected for processing.")
        return

    for config_item in species_configs_to_process:
        species_id = config_item["id_for_db"] 
        peptide_fasta_relative_path = config_item["out_peptide_fasta_name"].strip() 
        blast_filename = config_item["diamond_output_filename"].strip()

        peptide_fasta_file = os.path.join(args.proteome_dir, peptide_fasta_relative_path)
        blast_output_file = os.path.join(args.blast_dir, blast_filename)

        val_genes_for_species = process_species_for_validation_set(
            species_id, 
            peptide_fasta_file, 
            blast_output_file,
            args.num_validation_genes
        )
        all_validation_genes[species_id] = val_genes_for_species
    
    data_to_save = all_validation_genes

    if not data_to_save and args.target_species_id: 
        # Check if the specific target species yielded no data
        if args.target_species_id in data_to_save and not data_to_save[args.target_species_id]:
             print(f"No validation genes were generated for the targeted species {args.target_species_id}. Pickle file will not be saved/updated for this species based on current run.")
        # Or if data_to_save is empty overall (implies target species was not found or processed correctly if it was the only one)
        elif not data_to_save:
             print(f"No validation genes were generated for the targeted species {args.target_species_id}. Pickle file not saved as no data was generated.")
        # Fallback for not saving if there's nothing in data_to_save for the target.
        if not any(data_to_save.values()): # if all lists are empty or dict is empty
            print("Pickle file not saved as no data was generated overall.")
            return

    elif not data_to_save and not args.target_species_id:
        print(f"No validation genes were generated for any species. Pickle file will not be saved.")
        return
    elif not any(data_to_save.values()): # Final check if all processed species resulted in empty lists
        print("No validation genes were generated for any processed species. Pickle file will not be saved.")
        return

    output_pickle_file = args.output_file
    if args.target_species_id and os.path.exists(args.output_file) and len(species_configs_to_process) == 1:
        try:
            with open(args.output_file, 'rb') as pf_read:
                existing_data = pickle.load(pf_read)
            print(f"Loaded existing data from {args.output_file}")
            # Only update if the new data for the target species is non-empty
            if data_to_save.get(args.target_species_id):
                 existing_data.update(data_to_save) 
                 data_to_save = existing_data
                 print(f"Data for {args.target_species_id} will be updated in {args.output_file}.")
            else:
                print(f"No new data generated for {args.target_species_id}, existing pickle file unchanged for this key.")
                data_to_save = existing_data # keep existing data if no new data for target
        except Exception as e:
            print(f"Could not load or update existing pickle file {args.output_file}: {e}. A new file (or overwrite) for processed species will occur if there is new data.")
    
    output_dir = os.path.dirname(output_pickle_file)
    if output_dir and not os.path.exists(output_dir): 
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    print(f"\nSaving validation gene sets to: {output_pickle_file}")
    try:
        with open(output_pickle_file, 'wb') as pf:
            pickle.dump(data_to_save, pf, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"Successfully saved to {os.path.basename(output_pickle_file)}")
    except Exception as e:
        print(f"Error saving pickle file: {e}")

    print("\nScript finished.")

if __name__ == "__main__":
    main()
