#!/usr/bin/env python3

import argparse
from scipy.stats import spearmanr, rankdata, norm
from statsmodels.stats.multitest import multipletests
import pandas as pd
import numpy as np
import os
import logging
from datetime import datetime

# Define argument parser
def parse_args():
    parser = argparse.ArgumentParser(description="TU-TU Correlation Analysis")
    # Add the arguments you want to be available to the script
    parser.add_argument('--asv_file', required=True, help='Path to the ASV table file')
    parser.add_argument('--metadata_file', required=True, help='Path to the metadata file')
    parser.add_argument('--output_dir', required=True, help='Directory to save the output')
    parser.add_argument('--corr_threshold_high', required=True, type=float, help='High correlation threshold')
    parser.add_argument('--corr_threshold_low', required=True, type=float, help='Low correlation threshold')
    parser.add_argument('--pval_threshold', required=True, type=float, help='P-value threshold')
    parser.add_argument('--column', required=True, help='Column name in metadata for grouping')
    parser.add_argument('--treatment_col', required=True, help='Column name in metadata for treatment')
    parser.add_argument('--setA', required=True, help='Comma-separated list of samples for Set A')
    parser.add_argument('--setB_before', required=True, help='Comma-separated list of samples for Set B (Before)')
    parser.add_argument('--setB_after', required=True, help='Comma-separated list of samples for Set B (After)')
    
    return parser.parse_args()

# Setup logger
def setup_logger(output_dir, column):
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(output_dir, f'corr_otu_v_foldincr_{column}_{current_time}.log')
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    console.setFormatter(formatter)
    logger = logging.getLogger()
    if not logger.hasHandlers():
        logger.addHandler(console)

# Load ASV table
def load_asv_table(asv_file):
    with open(asv_file, 'r') as f:
        first_line = f.readline().strip()
        skiprows = [0] if first_line == "# Constructed from biom file" else []
    return pd.read_csv(asv_file, sep='\t', index_col=0, skiprows=skiprows)

# Load metadata
def load_metadata(metadata_file):
    return pd.read_csv(metadata_file, sep='\t', index_col=0)

def average_counts(asv_table, metadata, column, treatment_col):
    metadata = metadata.dropna(subset=[column])
    metadata = metadata.loc[metadata.index.intersection(asv_table.columns)]
    averaged_counts = {}
    grouped = metadata.groupby([column, treatment_col])
    for (group, treatment), samples in grouped:        
        relevant_samples = asv_table[samples.index]
        # Check if relevant_samples is empty BEFORE computing the mean
        if relevant_samples.empty:
            logging.warning(f"No valid samples for group {group}, treatment {treatment}. Skipping.")
            continue
        averaged_counts[(group, treatment)] = relevant_samples.mean(axis=1)
    # If no valid groups were found, return an empty DataFrame
    if not averaged_counts:
        logging.error("No valid samples found after filtering. Please check your inputs.")
        return pd.DataFrame()
    return pd.DataFrame(averaged_counts)

def calculate_setB_ratio(setB_before, setB_after):
    setB_before.columns = setB_before.columns.get_level_values(1)
    setB_after.columns = setB_after.columns.get_level_values(1)
    common_columns = setB_before.columns.intersection(setB_after.columns)
    setB_before_common = setB_before[common_columns]
    setB_after_common = setB_after[common_columns]
    ratio = setB_after_common / setB_before_common.replace(0, np.nan)
    return ratio

def calculate_spearman_matrix(setA_values, setB_values):
    setA_array = setA_values.to_numpy()
    setB_array = setB_values.to_numpy()
    if setA_array.size == 0 or setB_array.size == 0:
        raise ValueError("Empty input array provided for correlation calculation.")
    setA_ranks = np.apply_along_axis(rankdata, 1, setA_array)
    setB_ranks = np.apply_along_axis(rankdata, 1, setB_array)
    setA_ranks -= np.mean(setA_ranks, axis=1, keepdims=True)
    setB_ranks -= np.mean(setB_ranks, axis=1, keepdims=True)
    covariance_matrix = np.dot(setA_ranks, setB_ranks.T)
    std_A = np.sqrt(np.sum(setA_ranks ** 2, axis=1))
    std_B = np.sqrt(np.sum(setB_ranks ** 2, axis=1))
    std_A[std_A == 0] = np.nan
    std_B[std_B == 0] = np.nan
    denominator_matrix = np.outer(std_A, std_B)
    correlation_matrix = covariance_matrix / denominator_matrix
    correlation_matrix = np.nan_to_num(correlation_matrix)
    n = setA_array.shape[1]
    with np.errstate(divide='ignore', invalid='ignore'):
        t_stat = correlation_matrix * np.sqrt((n - 2) / (1 - correlation_matrix ** 2))
    p_values = 2 * (1 - norm.cdf(np.abs(t_stat)))
    return correlation_matrix, p_values

def run_analysis(setA_values, setB_values, otu_ids_A, output_file):
    setA_numeric_columns = setA_values.columns.get_level_values(1)
    matching_columns = set(setA_numeric_columns) & set(setB_values.columns)
    if not matching_columns:
        raise ValueError("No matching columns found between setA_values and setB_values.")
    # Subset both datasets to include only matching columns
    setA_matched = setA_values.loc[:, setA_values.columns.get_level_values(1).isin(matching_columns)]
    setB_matched = setB_values[list(matching_columns)]
    # Extract column identifiers (level 1 values)
    setA_col_ids = setA_matched.columns.get_level_values(1)
    setB_col_ids = setB_matched.columns
    # Count occurrences of each unique value in setA_col_ids
    setA_counts = pd.Series(setA_col_ids).value_counts()
    # Expand setB_matched by repeating rows to match occurrences in setA
    expanded_B_rows = []
    expanded_B_columns = []
    # Rebuild setB_expanded to match setA_matched
    setB_expanded = pd.DataFrame(index=setB_matched.index)
    if isinstance(setA_matched.columns, pd.MultiIndex):
        for col in setA_matched.columns.get_level_values(1):
            # Find all matching columns in setB_matched (can be multiple)
            matching_cols = setB_matched.columns[setB_matched.columns == col]
            if len(matching_cols) == 0:
                raise ValueError(f"Column {col} in setA_matched has no match in setB_matched")
            # Duplicate each matching column to align with setA_matched
            for i in range((setA_matched.columns.get_level_values(1) == col).sum()):
                setB_expanded[col] = setB_matched[matching_cols].iloc[:, 0]  # Take first match
    else:
        # If no MultiIndex, match directly
        common_cols = setA_matched.columns.intersection(setB_matched.columns)
        setB_expanded = setB_matched[common_cols]
    # Ensure column order matches
    setB_expanded = setB_expanded[setA_matched.columns.get_level_values(1)]
    # Check if shapes match
    assert setA_matched.shape == setB_expanded.shape, f"Shape mismatch: {setA_matched.shape} vs {setB_expanded.shape}"
    # Perform Spearman correlation
    correlation_matrix, p_values = calculate_spearman_matrix(setA_matched, setB_expanded)
    # Ensure p_values is 1D and handle them correctly
    p_values = p_values.flatten()  # Flatten to ensure it's 1D
    # Perform FDR correction
    _, corrected_p_values, _, _ = multipletests(p_values, method='fdr_bh')
    # Since p_values is 1D, we don't need to reshape it
    corrected_p_values = corrected_p_values.flatten()
    # Results collection
    results = []
    for i, otu1_id in enumerate(otu_ids_A):
        corr = correlation_matrix[i, 0]  # Assuming correlation_matrix is 2D, but with 1 column
        p_value = p_values[i]  # Use only one index to access p_values
        fdr = corrected_p_values[i]  # Use only one index to access corrected p-values
        if not np.isnan(corr):
            results.append([
                otu1_id,
                corr,
                p_value,
                fdr,
                ','.join(map(str, setA_values.iloc[i])),
                ','.join(map(str, setB_values.iloc[i])),
                ','.join(map(str, setA_values.columns.get_level_values(1)))
            ])
    # Create DataFrame and write to output
    results_df = pd.DataFrame(results, columns=['TU', 'Correlation Coefficient', 'P-Value', 'FDR', 'TU_Counts_SetA', 'TU_Counts_SetB_Ratio', 'Sample_IDs'])
    results_df.to_csv(output_file, sep='\t', index=False)
    return results_df


def main():
    args = parse_args()
    setup_logger(args.output_dir,args.column)
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    logging.info("Starting TU-TU correlation analysis.")
    asv_table = load_asv_table(args.asv_file)
    metadata = load_metadata(args.metadata_file)
    logging.info("Loaded ASV table and metadata.")
    averaged_counts = average_counts(asv_table, metadata, args.column, args.treatment_col)
    if averaged_counts.empty:
        logging.error("No valid samples found after filtering. Please check your inputs.")
        return
    setA_samples = [x.strip() for x in args.setA.split(',')]
    setB_before_samples = [x.strip() for x in args.setB_before.split(',')]
    setB_after_samples = [x.strip() for x in args.setB_after.split(',')]
    setA_values = averaged_counts.loc[:, averaged_counts.columns.get_level_values(0).isin(setA_samples)]
    setB_before = averaged_counts.loc[:, averaged_counts.columns.get_level_values(0).isin(setB_before_samples)]
    setB_after = averaged_counts.loc[:, averaged_counts.columns.get_level_values(0).isin(setB_after_samples)]
    
    # drop any that aren't the same between before and after
    setB_values = calculate_setB_ratio(setB_before, setB_after)
    otu_ids_A = setA_values.index
    output_file = os.path.join(args.output_dir, f"otu_correlation_results_{args.column}.txt")
    run_analysis(setA_values, setB_values, otu_ids_A, output_file)
    logging.info("Correlation analysis completed.")

if __name__ == "__main__":
    main()
