#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import argparse
import logging
import datetime

def setup_logger(output_dir):
    os.makedirs(output_dir, exist_ok=True)
    current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_file = os.path.join(output_dir, f'prism_maker_{current_time}.log')
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

parser = argparse.ArgumentParser(description="Generate PRISM-compatible OTU abundance summaries.")
parser.add_argument('--file_path', required=True, help="Path to the ASV table file (normalized).")
parser.add_argument('--map_file', required=True, help="Path to the mapping file.")
parser.add_argument('--column', required=True, help="Column in the mapping file to group by (e.g., 'SoilNumber').")
parser.add_argument('--output_dir', required=True, help="Directory to save output files.")
parser.add_argument('--num_taxa', type=int, default=3, help="Number of top taxa to include (default: 3).")
parser.add_argument('--output_mode', choices=['avg_all', 'split'], required=True,
                    help="Specify output mode: 'avg_all' for one summary row or 'split' for separate files.")
args = parser.parse_args()

setup_logger(args.output_dir)

# Extract an ID from file_path for naming
asv_table_id = os.path.splitext(os.path.basename(args.file_path))[0]

# Read ASV table
logging.info(f"Reading ASV table from: {args.file_path}")
with open(args.file_path) as f:
    first_line = f.readline().strip()
skip_rows = 1 if first_line == "# Constructed from biom file" else 0
df = pd.read_csv(args.file_path, sep='\t', skiprows=skip_rows, comment='~')
df.rename(columns={df.columns[0]: "ASV_ID"}, inplace=True)

if 'taxonomy' in df.columns:
    df['taxonomy'] = df['taxonomy'] + '_' + df['ASV_ID']
    df.drop(columns=['ASV_ID'], inplace=True)
else:
    df.rename(columns={'ASV_ID': 'taxonomy'}, inplace=True)

# Read map file
logging.info(f"Reading mapping file from: {args.map_file}")
map_df = pd.read_csv(args.map_file, sep='\t', comment='~')
map_df.rename(columns={'#SampleID': 'sampleID'}, inplace=True)
map_df = map_df[[args.column, 'sampleID']]
treatments = map_df[args.column].dropna().unique()

# Transpose ASV table
df_transposed = df.set_index('taxonomy').T
df_transposed.reset_index(inplace=True)
df_transposed.rename(columns={'index': 'sampleID'}, inplace=True)

# Merge ASV and metadata
merged_df = pd.merge(map_df, df_transposed, on='sampleID')
merged_df = merged_df[~pd.isna(merged_df[args.column])]

# Determine top N taxa per treatment
treatment_top_taxa_dict = {}
unique_top_taxa = set()
treatment_avg_abundance_dict = {}

for treatment in treatments:
    treatment_df = merged_df[merged_df[args.column] == treatment]
    treatment_abundance = treatment_df.drop(columns=['sampleID', args.column])
    avg_abundance = treatment_abundance.mean(axis=0).reset_index()
    avg_abundance.columns = ['taxonomy', 'average_abundance']
    avg_abundance = avg_abundance.sort_values(by='average_abundance', ascending=False)
    top_taxa = avg_abundance.head(args.num_taxa)['taxonomy'].tolist()
    treatment_top_taxa_dict[treatment] = top_taxa
    unique_top_taxa.update(top_taxa)
    treatment_avg_abundance_dict[treatment] = avg_abundance.set_index('taxonomy')

final_top_taxa = sorted(unique_top_taxa)

# Generate output
if args.output_mode == 'avg_all':
    treatment_means = pd.DataFrame([
        treatment_avg_abundance_dict[t].loc[:, 'average_abundance']
        for t in treatments
    ]).fillna(0)

    avg_row = {'treatment': 'avg_abun'}
    avg_row.update({otu: treatment_means[otu].mean() for otu in final_top_taxa})
    avg_row['other'] = treatment_means.drop(columns=final_top_taxa, errors='ignore').sum(axis=1).mean()

    avg_df = pd.DataFrame([avg_row])
    row_sum = avg_df.drop(columns='treatment').sum(axis=1).values[0]
    if row_sum > 0:
        avg_df.iloc[:, 1:] = (avg_df.drop(columns='treatment') / row_sum) * 100

    filename = f"prism.column_{args.column}.numtaxa_{args.num_taxa}.outputmode_avg_all.tsv"
    output_path = os.path.join(args.output_dir, filename)
    logging.info(f"Saving average-only output to: {output_path}")
    avg_df.to_csv(output_path, sep='\t', index=False)

elif args.output_mode == 'split':
    for treatment in treatments:
        avg_abundance = treatment_avg_abundance_dict[treatment]['average_abundance']
        row = {'treatment': treatment}
        row.update({otu: avg_abundance.get(otu, 0) for otu in final_top_taxa})
        row['other'] = avg_abundance.drop(index=final_top_taxa, errors='ignore').sum()
        treatment_df = pd.DataFrame([row])
        row_sum = treatment_df.drop(columns='treatment').sum(axis=1).values[0]
        if row_sum > 0:
            treatment_df.iloc[:, 1:] = (treatment_df.drop(columns='treatment') / row_sum) * 100

        filename = f"prism.column_{args.column}.numtaxa_{args.num_taxa}.outputmode_split.{treatment}.tsv"
        output_path = os.path.join(args.output_dir, filename)
        logging.info(f"Saving per-treatment file for {treatment} to: {output_path}")
        treatment_df.to_csv(output_path, sep='\t', index=False)

logging.info("Script complete.")
