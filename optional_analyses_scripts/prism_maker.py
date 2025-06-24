#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import argparse
import logging
from collections import Counter
import datetime

"""
Description:
This script processes an ASV/OTU table and a mapping file to generate average abundance 
data per treatment group and outputs a consolidated PRISM-compatible file. 
...
"""

# Set up the logger
def setup_logger(output_dir):
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    # Get the current date and time in the format "YYYY-MM-DD_HH-MM-SS"
    current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    # Create a log file name with the timestamp
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

# Set up argument parsing
parser = argparse.ArgumentParser(description="Process ASV table and mapping file to generate average abundance data.")
parser.add_argument('--file_path', required=True, help="Path to the ASV table file (normalized).")
parser.add_argument('--map_file', required=True, help="Path to the mapping file.")
parser.add_argument('--column', required=True, help="Column in the mapping file to group by (e.g., 'SoilNumber').")
parser.add_argument('--output_avg_abundance_directory', required=True, help="Directory to save average abundance files.")
parser.add_argument('--output_prism_file', required=True, help="Path for the split output (PRISM TSV file).")
parser.add_argument('--num_taxa', type=int, default=3, help="Number of top taxa to include (default: 3).")

# Parse the arguments
args = parser.parse_args()

# Setup logger
setup_logger(args.output_avg_abundance_directory)

# Read in data file (ASV table)
logging.info(f"Reading ASV table from: {args.file_path}")
with open(args.file_path) as f:
    first_line = f.readline().strip()

# Decide whether to skip the first row
skip_rows = 1 if first_line == "# Constructed from biom file" else 0

# Load the DataFrame, skipping the first line if necessary
df = pd.read_csv(args.file_path, sep='\t', skiprows=skip_rows, comment='~')
df.rename(columns={df.columns[0]: "ASV_ID"}, inplace=True)

# Check if 'taxonomy' exists and modify accordingly
if 'taxonomy' in df.columns:
    df['taxonomy'] = df['taxonomy'] + '_' + df['ASV_ID']
    df.drop(columns=['ASV_ID'], inplace=True)
else:
    df.rename(columns={'ASV_ID': 'taxonomy'}, inplace=True)

# Read map file (treatment information)
logging.info(f"Reading mapping file from: {args.map_file}")
map_df = pd.read_csv(args.map_file, sep='\t', comment='~')
map_df.rename(columns={'#SampleID': 'sampleID'}, inplace=True)
map_df = map_df[[args.column, 'sampleID']]  # Adjust based on the column name you want to use

# Get unique treatment values from the specified column
treatments = map_df[args.column].unique()

# Filter out any 'nan' values from treatments (Remove NaN entries)
treatments = treatments[~pd.isna(treatments)]  # Remove NaN entries

# Ensure the output directory exists
os.makedirs(args.output_avg_abundance_directory, exist_ok=True)
logging.info(f"Output directory created: {args.output_avg_abundance_directory}")

# Transpose the ASV table to have sample IDs as rows (aligned with mapping file)
df_transposed = df.set_index('taxonomy').T  # Taxonomy becomes column names, sampleIDs are rows
df_transposed.reset_index(inplace=True)  # Sample IDs move to a column named "index"
df_transposed.rename(columns={'index': 'sampleID'}, inplace=True)  # Rename the index column to 'sampleID'

# Merge with the mapping file and filter out rows where treatment is NaN
logging.info("Merging ASV data with treatment information from the map file")
merged_df = pd.merge(map_df, df_transposed, on='sampleID')  # Align treatments with ASV data
merged_df = merged_df[~pd.isna(merged_df[args.column])]  # Remove rows where treatment is NaN

# Identify top OTUs per treatment and collect unique taxa
treatment_top_taxa_dict = {}  # Dictionary to store top taxa per treatment
unique_top_taxa = set()  # Set to store all unique top taxa across treatments
treatment_avg_abundance_dict = {}  # Initialize dictionary to store average abundances per treatment

for treatment in treatments:
    logging.info(f"Processing treatment: {treatment}")
    treatment_df = merged_df[merged_df[args.column] == treatment]
    treatment_abundance = treatment_df.drop(columns=['sampleID', args.column])
    avg_abundance = treatment_abundance.mean(axis=0).reset_index()
    avg_abundance.columns = ['taxonomy', 'average_abundance']
    avg_abundance = avg_abundance.sort_values(by='average_abundance', ascending=False)
    # Get top N taxa for this treatment
    top_taxa = avg_abundance.head(args.num_taxa)['taxonomy'].tolist()
    treatment_top_taxa_dict[treatment] = top_taxa
    # Add these top taxa to the global unique set
    unique_top_taxa.update(top_taxa)

# Convert the set of unique top taxa into a sorted list
final_top_taxa = sorted(unique_top_taxa)  # Sorting for consistency

# Step 1.5: Compute and store the average abundance data for each treatment
for treatment in treatments:
    treatment_df = merged_df[merged_df[args.column] == treatment]
    treatment_abundance = treatment_df.drop(columns=['sampleID', args.column])
    avg_abundance = treatment_abundance.mean(axis=0).reset_index()
    avg_abundance.columns = ['taxonomy', 'average_abundance']
    avg_abundance = avg_abundance.set_index('taxonomy')
    # Store this in the dictionary
    treatment_avg_abundance_dict[treatment] = avg_abundance
    
# Step 2: Prepare the output DataFrame with dynamic columns
columns = ['treatment'] + final_top_taxa + ['other']
output_df = pd.DataFrame(columns=columns)

# Ensure correct column types
output_df = output_df.astype({col: float for col in columns if col != 'treatment'})
output_df['treatment'] = output_df['treatment'].astype(str)

# Step 3: Calculate average abundances for each treatment
for treatment in treatments:
    avg_abundance = treatment_avg_abundance_dict[treatment].set_index('taxonomy')['average_abundance']
    
    row = {'treatment': treatment}
    row.update({otu: avg_abundance.get(otu, 0) for otu in final_top_taxa})  # Fill missing taxa with 0
    row['other'] = avg_abundance.drop(index=final_top_taxa, errors='ignore').sum()  # Sum other taxa
    
    output_df = pd.concat([output_df, pd.DataFrame([row])], ignore_index=True)

# Step 4: Calculate the "all" row (average across treatments)
treatment_means = pd.DataFrame([
    treatment_avg_abundance_dict[treatment].set_index('taxonomy')['average_abundance']
    for treatment in treatments
]).fillna(0)

overall_row = {'treatment': 'avg_abun'}
overall_row.update({otu: treatment_means[otu].mean() for otu in final_top_taxa})
overall_row['other'] = treatment_means.drop(columns=final_top_taxa, errors='ignore').sum(axis=1).mean()

output_df = pd.concat([pd.DataFrame([overall_row]), output_df], ignore_index=True)

# Step 5: Convert values to proportions
for idx, row in output_df.iterrows():
    row_sum = row.drop(['treatment']).sum()
    if row_sum > 0:
        output_df.iloc[idx, 1:] = (row.drop(['treatment']) / row_sum) * 100

# Step 6: Save the updated DataFrame
logging.info(f"Saving the output PRISM file to: {args.output_prism_file}")
output_df.to_csv(args.output_prism_file, sep='\t', index=False)

logging.info("Processing complete. Results saved.")

