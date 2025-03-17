import argparse
from scipy.stats import spearmanr
import pandas as pd
import numpy as np
import os
import logging

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
    logging.getLogger().addHandler(console)

def load_asv_table(asv_file):
    with open(asv_file, 'r') as f:
        first_line = f.readline().strip()
        skiprows = [0] if first_line == "# Constructed from biom file" else []
    return pd.read_csv(asv_file, sep='\t', index_col=0, skiprows=skiprows)

def load_metadata(metadata_file):
    return pd.read_csv(metadata_file, sep='\t', index_col=0)

def main(asv_file, metadata_file, output_dir, column):
    setup_logger(output_dir)
    logging.info("Starting TU-Correlation analysis.")

    asv_table = load_asv_table(asv_file)
    metadata = load_metadata(metadata_file)

    common_samples = asv_table.columns.intersection(metadata.index)
    valid_samples = metadata.index[metadata[column].notna()]  # Only samples with values in the specified column
    filtered_samples = common_samples.intersection(valid_samples)

    filtered_asv_table = asv_table[filtered_samples]

    # Filter out TUs with constant counts across all valid samples
    constant_tus = []
    filtered_asv_table_no_constant = []

    for otu_id, otu_counts in filtered_asv_table.iterrows():
        if otu_counts.nunique() == 1:  # If the counts are constant
            constant_tus.append(otu_id)
        else:
            filtered_asv_table_no_constant.append(otu_counts)

    if constant_tus:
        logging.info(f"Removed TUs with constant counts: {', '.join(constant_tus)}")

    # Rebuild filtered ASV table without constant TUs
    filtered_asv_table = pd.DataFrame(filtered_asv_table_no_constant, columns=filtered_samples)

    results = []

    otu_ids = filtered_asv_table.index
    otu_count = len(otu_ids)

    # Now correlate each TU with the specified column (instead of pairwise TU correlation)
    for otu_id in otu_ids:
        otu_counts = filtered_asv_table.loc[otu_id].astype(float)
        corr, p_value = spearmanr(otu_counts, metadata.loc[filtered_samples, column])

        if np.isnan(corr) or np.isnan(p_value):
            continue

        results.append({
            'TU': otu_id,
            'Correlation Coefficient': corr,
            'P-Value': p_value,
            'TU_Counts': ','.join(map(str, otu_counts)),
            'Sample_IDs': ','.join(filtered_samples)
        })

    # Create the results DataFrame
    results_df = pd.DataFrame(results)
    results_df = results_df[['TU', 'Correlation Coefficient', 'P-Value', 'TU_Counts', 'Sample_IDs']]

    output_file = os.path.join(output_dir, f"correlation_results_{column}.txt")
    results_df.to_csv(output_file, sep='\t', index=False)

    logging.info(f"Results saved to {output_file}.")
    logging.info("Correlation analysis completed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TU correlation with a metadata column.")
    parser.add_argument('--asv_file', required=True, help='Path to ASV table file')
    parser.add_argument('--metadata_file', required=True, help='Path to metadata file')
    parser.add_argument('--output_dir', required=True, help='Directory to save output files')
    parser.add_argument('--column', required=True, help='Column showing which samples to correlate.')

    args = parser.parse_args()
    main(args.asv_file, args.metadata_file, args.output_dir, args.column)
