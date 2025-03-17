import argparse
from scipy.stats import spearmanr, rankdata, norm
from statsmodels.stats.multitest import multipletests
import pandas as pd
import numpy as np
import os
import logging

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
    parser.add_argument('--treatment_col', help='Column name in metadata for treatment (optional)')
    parser.add_argument('--setA', required=True, help='Comma-separated list of samples for Set A')
    parser.add_argument('--setB', required=True, help='Comma-separated list of samples for Set B')
    return parser.parse_args()

def setup_logger(output_dir):
    log_file = os.path.join(output_dir, 'correlation_analysis.log')
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

def run_analysis(setA_values, setB_values, otu_ids_A, otu_ids_B, output_file):
    correlation_matrix, p_value_matrix = calculate_spearman_matrix(setA_values, setB_values)
    p_values = p_value_matrix.flatten()
    p_values_no_nans = p_values[~np.isnan(p_values)]
    _, corrected_p_values, _, _ = multipletests(p_values_no_nans, method='fdr_bh')
    corrected_p_values_full = np.full_like(p_values, np.nan)
    corrected_p_values_full[~np.isnan(p_values)] = corrected_p_values
    corrected_p_values_full = corrected_p_values_full.reshape(p_value_matrix.shape)
    results = []
    for i, otu1_id in enumerate(otu_ids_A):
        for j, otu2_id in enumerate(otu_ids_B):
            corr = correlation_matrix[i, j]
            p_value = p_value_matrix[i, j]
            fdr = corrected_p_values_full[i, j]
            if not np.isnan(corr):
                results.append([
                    otu1_id,
                    otu2_id,
                    corr,
                    p_value,
                    fdr,
                    ','.join(map(str, setA_values.iloc[i])),
                    ','.join(map(str, setB_values.iloc[j])),
                    ','.join(map(str, setA_values.columns))
                ])
    results_df = pd.DataFrame(results, columns=['TU1', 'TU2', 'Correlation Coefficient', 'P-Value', 'FDR', 'TU1_Counts', 'TU2_Counts', 'Sample_IDs'])
    results_df.to_csv(output_file, sep='\t', index=False)
    return results_df

def main():
    args = parse_args()

    # Check if treatment_col is missing and setA/setB are different
    setA_samples = [x.strip() for x in args.setA.split(',')]
    setB_samples = [x.strip() for x in args.setB.split(',')]

    if args.treatment_col is None and setA_samples != setB_samples:
        logging.error("No treatment_col specified and Set A and Set B are different. Exiting analysis.")
        return

    setup_logger(args.output_dir)
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
    if not common_treatments:
        logging.error("No common treatments found between Set A and Set B. Exiting.")
        return

    setA_values = setA_values[list(common_treatments)]
    setB_values = setB_values[list(common_treatments)]

    setA_values = setA_values.loc[~(setA_values.nunique(axis=1) == 1)]
    setB_values = setB_values.loc[~(setB_values.nunique(axis=1) == 1)]

    otu_ids_A = setA_values.index
    otu_ids_B = setB_values.index

    output_file = os.path.join(args.output_dir, f"otu_correlation_results_{args.column}.txt")

    results_df = run_analysis(setA_values, setB_values, otu_ids_A, otu_ids_B, output_file)

    logging.info("Correlation analysis completed.")

if __name__ == "__main__":
    main()
