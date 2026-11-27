"""
Purpose:
    Classify genes into binary high (1), low (0), or unused (2) expression labels from
    normalized log10(mean TPM + 1) tables and emit class-threshold reports and histograms.

Manuscript (Peleke et al., 2026):
    Defines the immutable `classify_genes()` and `generate_report_and_plot()` helpers
    imported via `importlib` by `make_labels_dPHYlogic_meanTPM.0.py` and
    `make_label_reports_dPHYlogic_meanTPM.0.py`. Rule: class 1 > q-quantile (0.75);
    class 0 = matched count of genes (`<= 50%` zero-TPM, remainder lowest non-zero); class 2 = rest.

Usage & I/O:
    python code/01_data_preparation/02_classify_libaries_dPHY.0.3.py \
        --norm_type mean --file_stem leaf --q_threshold 0.75 --random_seed 2026 [--no_plots]
    - Reads `<BASE_RESULTS_DIR>/<species>/normalized/<species>_<file_stem>_log10_<norm_type>_tpm.csv`
    - Writes labeled CSVs, classification reports, and histograms under `LABELED_OUTPUT_SUBDIR`

How to change or extend:
    Pass different `--q_threshold`, `--norm_type`, or `--random_seed` values on the CLI,
    or update `classify_genes()` to implement a new expression binning rule.
"""
import os
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import seaborn as sns
import random
import argparse
import csv # Already imported, good.

# --- Configuration ---
BASE_RESULTS_DIR = "../data/results/"  # Where species folders are
NORMALIZED_INPUT_SUBDIR = "normalized" # Subdir where normalized input files are
LABELED_OUTPUT_SUBDIR = "labeled" # Subdir for outputs of this script (CSV, reports, plots)
REPORT_SUBDIR = "reports" # Subdir for text reports within LABELED_OUTPUT_SUBDIR
PLOTS_SUBDIR = "plots" # Subdir for plots within LABELED_OUTPUT_SUBDIR

GENE_ID_COL = "gene_id" # Standardized gene ID column name in input and output
TARGET_COL = "target" # Name for the classification target column

# Default classification percentile threshold
DEFAULT_Q_THRESHOLD = 0.8
# Default seed for random sampling
DEFAULT_RANDOM_SEED = 123

def load_normalized_file(filepath):
    """
    Loads a polished normalized TPM file.
    The file is expected to have 'gene_id' and a value column (e.g., log10MeanTPM).
    All original columns will be preserved.
    """
    try:
        df = pd.read_csv(filepath) # Let pandas infer quoting, or specify if needed
        if GENE_ID_COL not in df.columns:
            print(f"  Error: Column '{GENE_ID_COL}' not found in {filepath}. Skipping.")
            return None, None
        
        if df.shape[1] < 2: 
            print(f"  Error: File {filepath} has fewer than 2 columns. Expected '{GENE_ID_COL}' and a value column. Skipping.")
            return None, None
            
        if df.columns[0] == GENE_ID_COL:
            value_col_name = df.columns[1]
        else:
            print(f"  Warning: '{GENE_ID_COL}' is not the first column in {filepath}. Using second column '{df.columns[1]}' as metric. Ensure this is correct.")
            value_col_name = df.columns[1]

        if value_col_name == GENE_ID_COL:
             print(f"  Error: Value column name is the same as GENE_ID_COL ('{value_col_name}') in {filepath}. Skipping.")
             return None, None
        
        if not pd.api.types.is_numeric_dtype(df[value_col_name]):
            try:
                df[value_col_name] = pd.to_numeric(df[value_col_name], errors='raise')
            except ValueError:
                print(f"  Error: Value column '{value_col_name}' in {filepath} could not be converted to numeric. Skipping.")
                return None, None
        
        df.set_index(GENE_ID_COL, inplace=True)
        return df, value_col_name # Returns full DataFrame with gene_id as index

    except FileNotFoundError:
        print(f"  Error: Normalized file not found: {filepath}. Skipping.")
        return None, None
    except pd.errors.EmptyDataError:
        print(f"  Warning: Normalized file is empty: {filepath}. Skipping.")
        return None, None
    except Exception as e:
        print(f"  Error processing normalized file {filepath}: {e}. Skipping.")
        return None, None

def classify_genes(df_with_metric, metric_col_name, q_thresh, random_seed_val):
    df = df_with_metric.copy() 

    if df.empty or metric_col_name not in df.columns:
        print("  Error: DataFrame is empty or metric column is missing for classification.")
        return None, None, None

    upper_q_value = df[metric_col_name].quantile(q_thresh)
    
    df[TARGET_COL] = 2 

    high_expr_mask = df[metric_col_name] > upper_q_value
    df.loc[high_expr_mask, TARGET_COL] = 1
    count_class1 = high_expr_mask.sum()
    
    selected_A_indices = []
    actual_A_len = 0
    if count_class1 > 0:
        genes_metric_zero = df[df[metric_col_name] == 0]
        target_A_max_contribution = count_class1 // 2
        
        if not genes_metric_zero.empty and target_A_max_contribution > 0:
            actual_A_len = min(target_A_max_contribution, len(genes_metric_zero))
            if actual_A_len > 0:
                selected_A_indices = genes_metric_zero.sample(n=actual_A_len, random_state=random_seed_val).index.tolist()
    
    needed_from_B = count_class1 - actual_A_len
    selected_B_indices = []
    if needed_from_B > 0 and count_class1 > 0:
        # Ensure potential_B_genes are not already Class 1 and have metric > 0 and <= upper_q_value
        potential_B_genes_mask = (df[metric_col_name] > 0) & \
                                 (df[metric_col_name] <= upper_q_value) & \
                                 (df[TARGET_COL] != 1) # Not already classified as high expression
        potential_B_genes = df[potential_B_genes_mask]
        
        if not potential_B_genes.empty:
            actual_B_len = min(needed_from_B, len(potential_B_genes))
            if actual_B_len > 0:
                selected_B_indices = potential_B_genes.sort_values(by=metric_col_name).head(actual_B_len).index.tolist()
    
    class0_combined_idx = list(set(selected_A_indices + selected_B_indices))
    
    if class0_combined_idx: 
        # Ensure not overwriting Class 1 genes if any overlap (shouldn't happen with correct logic)
        df.loc[df.index.isin(class0_combined_idx) & (df[TARGET_COL] != 1), TARGET_COL] = 0
            
    return df, selected_A_indices, selected_B_indices


def generate_report_and_plot(df_classified, metric_col_name, species_name, input_file_stem_for_report, q_threshold, 
                             selected_A_idx, selected_B_idx,
                             output_base_labeled_path, generate_plots_flag):
    if df_classified is None or df_classified.empty:
        print("  Skipping report and plot generation due to empty/missing classified data.")
        return

    report_dir = os.path.join(output_base_labeled_path, REPORT_SUBDIR)
    os.makedirs(report_dir, exist_ok=True)
    
    upper_quantile_val = df_classified[metric_col_name].quantile(q_threshold)
    lower_quantile_ref_val = df_classified[metric_col_name].quantile(1 - q_threshold) 
    
    count_class1_final = (df_classified[TARGET_COL] == 1).sum()
    count_class0_final = (df_classified[TARGET_COL] == 0).sum()
    count_class2_final = (df_classified[TARGET_COL] == 2).sum()

    genes_actually_zero_metric_in_class_0 = 0
    if selected_A_idx: 
        valid_A_idx = df_classified.index.intersection(selected_A_idx) 
        if not valid_A_idx.empty:
            genes_actually_zero_metric_in_class_0 = (df_classified.loc[valid_A_idx, metric_col_name] == 0).sum()

    max_val_sel_B = np.nan
    if selected_B_idx: 
        valid_B_idx = df_classified.index.intersection(selected_B_idx)
        if not valid_B_idx.empty:
            max_val_sel_B = df_classified.loc[valid_B_idx, metric_col_name].max()

    count_filtered_class2_specific = 0
    # Class 2 genes are those not in Class 0 or Class 1.
    # Filtered Class 2 specific: Metric > max_val_sel_B (if B exists, else > 0) AND Metric < upper_quantile_val
    class2_mask = df_classified[TARGET_COL] == 2
    lower_bound_class2 = 0
    if pd.notna(max_val_sel_B):
        lower_bound_class2 = max_val_sel_B
    
    count_filtered_class2_specific = (class2_mask & 
                                     (df_classified[metric_col_name] > lower_bound_class2) &
                                     (df_classified[metric_col_name] < upper_quantile_val)).sum()


    report_base_name = f"{species_name}_{input_file_stem_for_report}_Q{q_threshold}"
    report_filename = f"{report_base_name}_classification_report.txt"
    report_filepath = os.path.join(report_dir, report_filename)

    # ... (rest of the reporting logic remains the same) ...
    with open(report_filepath, "w") as f:
        f.write(f"EXPRESSION CLASS REPORT for: {species_name} (using input stem: {input_file_stem_for_report})\n")
        f.write(f"Metric Column Used: {metric_col_name}\n")
        f.write(f"Q Threshold for Upper Percentile: {q_threshold*100:.1f}%\n")
        f.write(f"Calculated Upper Cutoff ({metric_col_name}): {upper_quantile_val:.4f}\n")
        lqr_val_formatted = "N/A"
        if pd.notna(lower_quantile_ref_val):
            try: lqr_val_formatted = f"{float(lower_quantile_ref_val):.4f}"
            except (ValueError, TypeError): lqr_val_formatted = str(lower_quantile_ref_val)
        elif pd.isna(lower_quantile_ref_val): lqr_val_formatted = "NaN"
        percentile_description_text = f"{(1 - q_threshold) * 100:.1f}%th percentile"
        f.write(f"Reference Lower Cutoff ({percentile_description_text}): {lqr_val_formatted}\n")
        f.write(f"Total Genes Classified as Class 1 (High Expression): {count_class1_final}\n")
        f.write(f"Total Genes Classified as Class 0 (Low/No Expression): {count_class0_final}\n")
        len_A = len(selected_A_idx) if selected_A_idx is not None else 0
        len_B = len(selected_B_idx) if selected_B_idx is not None else 0
        f.write(f"  - From 'metric == 0' pool (Subgroup A): {len_A} (of which {genes_actually_zero_metric_in_class_0} had metric strictly 0)\n")
        max_val_sel_B_formatted = "N/A"
        if pd.notna(max_val_sel_B):
            try: max_val_sel_B_formatted = f"{float(max_val_sel_B):.4f}"
            except (ValueError, TypeError): max_val_sel_B_formatted = str(max_val_sel_B)
        elif pd.isna(max_val_sel_B): max_val_sel_B_formatted = "NaN"
        if selected_B_idx is not None and len(selected_B_idx) > 0 :
            f.write(f"  - From 'lowest positive metric' pool (Subgroup B): {len_B}\n")
            f.write(f"Lower Percentile Cutoff after balancing (Max value in Subgroup B): {max_val_sel_B_formatted}\n")
        else:
            f.write(f"  - From 'lowest positive metric' pool (Subgroup B): 0\n")
            f.write(f"Lower Percentile Cutoff after balancing (Max value in Subgroup B): N/A (Subgroup B empty or not formed)\n")
        f.write(f"Total Genes Classified as Class 2 (Medium Expression): {count_class2_final}\n")
        f.write(f"Genes in Class 2 between effective lower and upper cutoffs: {count_filtered_class2_specific}\n")

    print(f"    Report saved to: {report_filepath}")

    if generate_plots_flag:
        # ... (plotting logic remains the same) ...
        plot_dir = os.path.join(output_base_labeled_path, PLOTS_SUBDIR)
        os.makedirs(plot_dir, exist_ok=True)
        plot_filename = f"{report_base_name}_classification_histogram.png"
        plot_filepath = os.path.join(plot_dir, plot_filename)
        plt.figure(figsize=(12, 7))
        plot_df = df_classified.copy()
        hist_plot_ax = sns.histplot(data=plot_df, x=metric_col_name, hue=TARGET_COL,
                                    palette={0: "darkcyan", 1: "darkorange", 2: "grey"},
                                    multiple="stack", bins=100, log_scale=(False, True))
        plt.yscale('log')
        current_ylim = plt.gca().get_ylim()
        plt.ylim(max(1, current_ylim[0]), current_ylim[1] * 1.1)
        plt.title(f"Distribution of {metric_col_name} for {species_name}\n(Input: {input_file_stem_for_report}, Q={q_threshold})")
        plt.xlabel(metric_col_name)
        plt.ylabel("Number of Genes (Log Scale)")
        legend = hist_plot_ax.get_legend()
        if legend:
            legend.set_title("Expression Class")
            new_labels = []
            label_map = {'0': '0: Low/No', '1': '1: High', '2': '2: Medium'}
            float_label_map = {'0.0': '0: Low/No', '1.0': '1: High', '2.0': '2: Medium'}
            for t in legend.texts:
                text_val = t.get_text()
                if text_val in label_map: new_labels.append(label_map[text_val])
                elif text_val in float_label_map: new_labels.append(float_label_map[text_val])
                else: new_labels.append(text_val)
            if len(new_labels) == len(legend.texts):
                for i, t in enumerate(legend.texts): t.set_text(new_labels[i])
            else: print("Warning: Legend label mapping was inconsistent. Using default legend labels.")
        else: print("Warning: Could not retrieve legend from Seaborn plot to customize.")
        plt.tight_layout()
        plt.savefig(plot_filepath)
        plt.close()
        print(f"    Plot saved to: {plot_filepath}")

    else:
        print("    Plot generation skipped as per --no_plots flag.")


def main():
    parser = argparse.ArgumentParser(description="Classify genes based on normalized TPM values.")
    # ... (args parsing remains the same) ...
    parser.add_argument("--norm_type", type=str, required=True, choices=["mean", "median", "max"],
                        help="Type of normalization used for the input files (e.g., 'mean' for '_log10_mean_tpm.csv').")
    parser.add_argument("--file_stem", type=str, required=True,
                        help="Stem of the input file (e.g., 'all_tissues', 'leaf'). This is the part between species_name_ and _log10_normtype_tpm.csv.")
    parser.add_argument("--q_threshold", type=float, default=DEFAULT_Q_THRESHOLD,
                        help=f"Quantile threshold for high expression classification (default: {DEFAULT_Q_THRESHOLD}).")
    parser.add_argument("--random_seed", type=int, default=DEFAULT_RANDOM_SEED,
                        help=f"Random seed for reproducibility (default: {DEFAULT_RANDOM_SEED}).")
    parser.add_argument("--no_plots", action="store_true", 
                        help="If set, plot generation will be skipped.")
    
    args = parser.parse_args()

    random.seed(args.random_seed)
    np.random.seed(args.random_seed)
    generate_plots = not args.no_plots 

    print(f"Starting gene classification script with Q_THRESHOLD = {args.q_threshold}, NormType = {args.norm_type}, FileStem = {args.file_stem}, GeneratePlots = {generate_plots}...")

    try:
        species_dirs = [
            d for d in os.listdir(BASE_RESULTS_DIR)
            if os.path.isdir(os.path.join(BASE_RESULTS_DIR, d)) and \
               d not in [LABELED_OUTPUT_SUBDIR, NORMALIZED_INPUT_SUBDIR, "polished_normalized", REPORT_SUBDIR, PLOTS_SUBDIR] # Assuming these are not species names
        ]
    except FileNotFoundError:
        print(f"Error: Base results directory not found: {BASE_RESULTS_DIR}.")
        return
    if not species_dirs:
        print(f"No species directories found in {BASE_RESULTS_DIR}. Exiting.")
        return

    for species_name in species_dirs:
        print(f"\nProcessing species: {species_name}")
        
        input_filename_pattern = f"{species_name}_{args.file_stem}_log10_{args.norm_type}_tpm.csv"
        normalized_filepath = os.path.join(BASE_RESULTS_DIR, species_name, NORMALIZED_INPUT_SUBDIR, input_filename_pattern)
        
        if not os.path.exists(normalized_filepath):
            print(f"  Normalized file not found: {normalized_filepath}. Skipping this species/file combination.")
            continue

        df_with_metric, metric_col_name = load_normalized_file(normalized_filepath)
        if df_with_metric is None:
            continue
            
        df_classified, sel_A_idx, sel_B_idx = classify_genes(df_with_metric, metric_col_name, args.q_threshold, args.random_seed)
        
        if df_classified is None:
            print(f"  Classification failed for {normalized_filepath}. Skipping output generation.")
            continue

        species_labeled_output_path = os.path.join(BASE_RESULTS_DIR, species_name, LABELED_OUTPUT_SUBDIR)
        os.makedirs(species_labeled_output_path, exist_ok=True)
        
        output_csv_basename = f"{species_name}_{args.file_stem}_log10_{args.norm_type}_tpm_Q{args.q_threshold}_labeled.csv"
        output_csv_filepath = os.path.join(species_labeled_output_path, output_csv_basename)
        
        df_to_save = df_classified.reset_index()
        
        # Refined column ordering
        # Ensures gene_id is first, target is last, metric_col_name is second to last
        # All other columns from the original input are preserved in between.
        
        current_cols = list(df_to_save.columns)
        
        # Start with gene_id
        ordered_cols = [GENE_ID_COL]
        
        # Add other columns, excluding those we will place specifically
        other_original_cols = [
            col for col in current_cols 
            if col not in [GENE_ID_COL, metric_col_name, TARGET_COL]
        ]
        ordered_cols.extend(other_original_cols)
        
        # Add metric_col_name if it's not already GENE_ID_COL (it shouldn't be) 
        # and ensure it's present before target.
        # This also handles if metric_col_name was one of the 'other_original_cols'.
        # To ensure it's distinctly placed before TARGET_COL and not duplicated:
        if metric_col_name in ordered_cols: # if it was part of other_original_cols
            ordered_cols.remove(metric_col_name)
        ordered_cols.append(metric_col_name)
        
        ordered_cols.append(TARGET_COL)

        # Ensure all original columns are present in the final list, prevent duplicates
        final_columns = []
        seen = set()
        for col in ordered_cols:
            if col in current_cols and col not in seen: # ensure col actually exists in df_to_save
                final_columns.append(col)
                seen.add(col)
        # Add any remaining columns from df_to_save not captured (should be none if logic is correct)
        for col in current_cols:
            if col not in seen:
                final_columns.append(col)
        
        df_to_save = df_to_save[final_columns]
        
        try:
            # --- MODIFICATION: Change quoting to QUOTE_NONNUMERIC ---
            df_to_save.to_csv(output_csv_filepath, index=False, quoting=csv.QUOTE_NONNUMERIC)
            print(f"  Labeled data saved to: {output_csv_filepath}")
        except Exception as e:
            print(f"  Error saving labeled CSV to {output_csv_filepath}: {e}")

        input_file_stem_for_report = f"{args.file_stem}_log10_{args.norm_type}_tpm"
        generate_report_and_plot(df_classified, metric_col_name, species_name, input_file_stem_for_report, 
                                 args.q_threshold, sel_A_idx, sel_B_idx,
                                 species_labeled_output_path, generate_plots) 

    print("\nGene classification script finished.")

if __name__ == "__main__":
    main()
