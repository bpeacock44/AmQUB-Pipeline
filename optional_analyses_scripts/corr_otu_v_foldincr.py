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
def setup_logger(output_dir column):
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(output_dir, f'corr_otu_v_foldincr_{column}_{current_time}_log.log')
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

def average_counts(asv_table, metadata, group_col, treatment_col):
    metadata = metadata.dropna(subset=[group_col])
    metadata = metadata.loc[metadata.index.intersection(asv_table.columns)]
    averaged_counts = {}
    grouped = metadata.groupby([group_col, treatment_col])
    for (group, treatment), samples in grouped:
        relevant_samples = asv_table[samples.index]
        if not relevant_samples.empty:
            averaged_counts[(group, treatment)] = relevant_samples.mean(axis=1)
    return pd.DataFrame(averaged_counts)

def calculate_setB_ratio(setB_before, setB_after):
    setB_before.columns = setB_before.columns.get_level_values(1)
    setB_after.columns = setB_after.columns.get_level_values(1)
    common_columns = setB_before.columns.intersection(setB_after.columns)
    setB_before_common = setB_before[common_columns]
    setB_after_common = setB_after[common_columns]
    ratio = setB_after_common / setB_before_common
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
    correlation_matrix, p_value_matrix = calculate_spearman_matrix(setA_values, setB_values)
    p_values = p_value_matrix.flatten()
    _, corrected_p_values, _, _ = multipletests(p_values, method='fdr_bh')
    corrected_p_values = corrected_p_values.reshape(p_value_matrix.shape)
    results = []
    for i, otu1_id in enumerate(otu_ids_A):
        corr = correlation_matrix[i, 0]
        p_value = p_value_matrix[i, 0]
        fdr = corrected_p_values[i, 0]
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
    results_df = pd.DataFrame(results, columns=['TU', 'Correlation Coefficient', 'P-Value', 'FDR', 'TU_Counts_SetA', 'TU_Counts_SetB_Ratio', 'Sample_IDs'])
    results_df.to_csv(output_file, sep='\t', index=False)
    return results_df

def main():
    args = parse_args()
    setup_logger(args.output_dir,args.column)
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
    setB_ratio = calculate_setB_ratio(setB_before, setB_after)

    otu_ids_A = setA_values.index

    output_file = os.path.join(args.output_dir, f"otu_correlation_results_{args.column}.txt")

    run_analysis(setA_values, setB_ratio, otu_ids_A, output_file)
    logging.info("Correlation analysis completed.")


if __name__ == "__main__":
    main()
