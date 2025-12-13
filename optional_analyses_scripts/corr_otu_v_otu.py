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
    parser.add_argument('--asv_file', required=True, help='Path to the ASV table file')
    parser.add_argument('--metadata_file', required=True, help='Path to the metadata file')
    parser.add_argument('--output_dir', required=True, help='Directory to save the output')
    parser.add_argument('--corr_threshold_high', required=True, type=float, help='High correlation threshold')
    parser.add_argument('--corr_threshold_low', required=True, type=float, help='Low correlation threshold')
    parser.add_argument('--pval_threshold', required=True, type=float, help='P-value threshold')
    parser.add_argument('--column', required=True, help='Column name in metadata for grouping')
    parser.add_argument('--treatment_col', help='Column name in metadata for treatment (optional)')
    parser.add_argument('--setA', required=True, help='Comma-separated list of samples for Set A')
    parser.add_argument('--setB', required=True, help='Comma-separated list of samples for Set B')
    parser.add_argument('--min_avg_abundance', required=False, type=float, default=0.0,
                        help='Minimum average abundance across the columns used for correlation. OTUs with mean < this will be removed.')
    return parser.parse_args()

def setup_logger(output_dir, column):
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(output_dir, f'corr_otu_v_otu_{column}_{current_time}.log')
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

def load_asv_table(asv_file):
    with open(asv_file, 'r') as f:
        first_line = f.readline().strip()
        skiprows = [0] if first_line == "# Constructed from biom file" else []
    return pd.read_csv(asv_file, sep='\t', index_col=0, skiprows=skiprows)

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

def calculate_spearman_matrix(setA_values, setB_values):
    setA_array = setA_values.to_numpy()
    setB_array = setB_values.to_numpy()
    n = setA_array.shape[1]
    if n < 3:
        raise ValueError(f"Spearman correlation requires at least 3 samples (n={n}).")
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
    with np.errstate(divide='ignore', invalid='ignore'):
        t_stat = correlation_matrix * np.sqrt((n - 2) / (1 - correlation_matrix ** 2))
    p_values = 2 * (1 - norm.cdf(np.abs(t_stat)))
    return correlation_matrix, p_values

def run_analysis(setA_values, setB_values, otu_ids_A, otu_ids_B, output_file,
                 corr_threshold_high=None, corr_threshold_low=None, pval_threshold=None):
    """
    Compute spearman correlations between setA_values and setB_values and write filtered results to output_file.

    Behavior:
    - If otu_ids_A and otu_ids_B are identical (same length and same order), only the upper triangle (i < j)
      is kept to avoid duplicate pairs and self-correlations.
    - If the OTU sets differ, all pairwise comparisons are kept, but exact self-pairs (same OTU ID in both positions)
      are removed.
    """

    # Compute Spearman correlation matrix and p-values
    correlation_matrix, p_value_matrix = calculate_spearman_matrix(setA_values, setB_values)

    nA, nB = correlation_matrix.shape

    # Create index grids (i for rows in A, j for cols in B)
    i_idx, j_idx = np.meshgrid(np.arange(nA), np.arange(nB), indexing='ij')

    # Flatten matrices
    corr_flat = correlation_matrix.ravel()
    pval_flat = p_value_matrix.ravel()
    otu1_flat = np.array(otu_ids_A)[i_idx.ravel()]
    otu2_flat = np.array(otu_ids_B)[j_idx.ravel()]

    # Determine if the two OTU lists are identical (same length AND same order)
    identical_sets = (nA == nB) and np.array_equal(np.array(otu_ids_A), np.array(otu_ids_B))

    # Upper-triangle mask: if comparing a set to itself (identical_sets), keep only i < j
    # Otherwise keep all pairs (we will still remove exact self-pairs below)
    if identical_sets:
        mask_upper = i_idx.ravel() < j_idx.ravel()
    else:
        mask_upper = np.ones_like(corr_flat, dtype=bool)

    # Always remove exact self-pairs where the OTU ID strings are identical
    mask_no_self = otu1_flat != otu2_flat
    
    # Threshold masks
    
    # Start with all False, because we'll OR positive/negative filters into it
    mask_corr = np.zeros_like(corr_flat, dtype=bool)
    
    # Strong positive correlations
    if corr_threshold_high is not None:
        mask_corr |= (corr_flat >= corr_threshold_high)
    # Strong negative correlations
    if corr_threshold_low is not None:
        mask_corr |= (corr_flat <= corr_threshold_low)
    # If no corr thresholds were provided, allow all correlations
    if (corr_threshold_high is None) and (corr_threshold_low is None):
        mask_corr = np.ones_like(corr_flat, dtype=bool)
    # Apply p-value threshold (this must be AND)
    if pval_threshold is not None:
        mask_corr &= (pval_flat <= pval_threshold)
    
    mask_thresh = mask_corr

    # Combine masks and remove NaNs
    final_mask = mask_upper & mask_no_self & mask_thresh & (~np.isnan(corr_flat))

    # indices of kept flattened entries (for constructing count strings)
    kept_i = i_idx.ravel()[final_mask]
    kept_j = j_idx.ravel()[final_mask]

    # Build final results DataFrame
    results_df = pd.DataFrame({
        'TU1': otu1_flat[final_mask],
        'TU2': otu2_flat[final_mask],
        'Correlation Coefficient': corr_flat[final_mask],
        'P-Value': pval_flat[final_mask],
        # TU1_Counts and TU2_Counts are stored as comma-separated strings of the sample values
        'TU1_Counts': [','.join(map(str, setA_values.iloc[i].tolist())) for i in kept_i],
        'TU2_Counts': [','.join(map(str, setB_values.iloc[j].tolist())) for j in kept_j],
        'Sample_IDs': ','.join(map(str, setA_values.columns))
    })

    # Write results
    results_df.to_csv(output_file, sep='\t', index=False)

    return results_df



def main():
    args = parse_args()

    setA_samples = [x.strip() for x in args.setA.split(',')]
    setB_samples = [x.strip() for x in args.setB.split(',')]

    if args.treatment_col is None and setA_samples != setB_samples:
        logging.error("No treatment_col specified and Set A and Set B are different. Exiting analysis.")
        return

    setup_logger(args.output_dir,args.column)
    logging.info("Starting TU-TU correlation analysis.")
    asv_table = load_asv_table(args.asv_file)
    metadata = load_metadata(args.metadata_file)
    logging.info("Loaded ASV table and metadata.")

    if args.treatment_col:
        averaged_counts = average_counts(asv_table, metadata, args.column, args.treatment_col)
        if averaged_counts.empty:
            logging.error("No valid samples found after filtering. Please check your inputs.")
            #return
        setA_values = averaged_counts.loc[:, averaged_counts.columns.get_level_values(0).isin(setA_samples)]
        setB_values = averaged_counts.loc[:, averaged_counts.columns.get_level_values(0).isin(setB_samples)]
        setA_values = setA_values.T.groupby(level=1).mean().T
        setB_values = setB_values.T.groupby(level=1).mean().T
    else:
        # Directly correlate without averaging if treatment_col is not used
        drop_cols = [col for col in asv_table.columns if col.lower() in ["sequence", "taxonomy"]]
        asv_table = asv_table.drop(columns=drop_cols)
        if asv_table.empty:
            logging.error("No valid samples found after filtering. Please check your inputs.")
            #return
        setA_values = asv_table
        setB_values = asv_table
    common_treatments = set(setA_values.columns).intersection(setB_values.columns)

    if len(common_treatments) < 3:
        logging.error(f"Not enough samples for correlation (n={len(common_treatments)}). Need at least 3. Exiting.")
        return

    setA_values = setA_values[list(common_treatments)]
    setB_values = setB_values[list(common_treatments)]

    if args.min_avg_abundance and args.min_avg_abundance > 0.0:
        # compute mean across columns (these are the columns used in correlation)
        meanA = setA_values.mean(axis=1)
        meanB = setB_values.mean(axis=1)

        # log counts before filtering
        logging.info(f"OTUs before abundance filter: SetA={setA_values.shape[0]}, SetB={setB_values.shape[0]}")

        # Keep OTUs with mean >= threshold (applied separately to each set)
        setA_values = setA_values.loc[meanA >= args.min_avg_abundance]
        setB_values = setB_values.loc[meanB >= args.min_avg_abundance]

        logging.info(f"Applied min_avg_abundance = {args.min_avg_abundance}. OTUs after filter: SetA={setA_values.shape[0]}, SetB={setB_values.shape[0]}")

        if setA_values.empty or setB_values.empty:
            logging.error("After applying min_avg_abundance filter, Set A or Set B is empty. Exiting.")
            return

    # remove constant rows
    setA_values = setA_values.loc[~(setA_values.nunique(axis=1) == 1)]
    setB_values = setB_values.loc[~(setB_values.nunique(axis=1) == 1)]

    otu_ids_A = setA_values.index
    otu_ids_B = setB_values.index

    output_file = os.path.join(args.output_dir, f"otu_correlation_results_{args.column}.txt")

    results_df = run_analysis(
        setA_values, setB_values, otu_ids_A, otu_ids_B, output_file,
        corr_threshold_high=args.corr_threshold_high,
        corr_threshold_low=args.corr_threshold_low,
        pval_threshold=args.pval_threshold
    )
    
    logging.info("Correlation analysis completed.")

if __name__ == "__main__":
    main()
