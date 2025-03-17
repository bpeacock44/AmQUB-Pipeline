#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
from collections import Counter
import argparse

"""
Description:
This script processes an ASV/OTU table and a mapping file to generate average abundance 
data per treatment group and outputs a consolidated PRISM-compatible file. 

The main steps include:
1. Reading the ASV table and mapping file.
2. Calculating average abundances for each treatment group.
3. Identifying the top taxa (ASVs/OTUs) based on abundance.
4. Generating an "all" row for overall average abundance across all treatments.
5. Saving average abundances for each treatment and an overall output file in 
   PRISM-compatible format.

Usage:
This script accepts command-line arguments for paths to the input files, output 
directory, and other parameters. Run the script as follows:

    ./process_asv.py --file_path <path_to_asv_table> \
        --map_file <path_to_mapping_file> \
        --column <column_name_for_grouping> \
        --output_avg_abundance_directory <output_directory_for_abundances> \
        --output_prism_file <path_to_prism_output> \
        --num_taxa <number_of_top_taxa>

Example:
    ./process_asv.py --file_path asv_table.txt \
        --map_file map_file.txt \
        --column SoilNumber \
        --output_avg_abundance_directory ./output/avg_abundances \
        --output_prism_file ./output/prism_output.tsv \
        --num_taxa 5

Input Notes:
1. The ASV/OTU table (file_path) should be normalized before beginning. 
2. Currently the input mapping file (map_file) is expected to have "#SampleID" as the sample ID column header. 
3. The column indicated by "column" should contain values only in rows of data that is to be included. 
    If all values in column are the same, then all samples indicated will be considered to be part of 
    the same treatment. If there are different values, samples with the same value will be aggregated 
    (averaged). The final output prism file will have a row for each treatment as well as an "all" row 
    for average abundance across all taxa.
4. Num taxa will determine how many top taxa per treatment to include in the plot. 
    (e.g. if there are 2 treatments and num_taxa is set to 5, then the top 5 taxa from each treatment will be 
    excluded from the "other" category and given their own column. If there are duplicates between the two treatments
    then they will only be included once. 

"""

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

# Read in data file (ASV table)
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
map_df = pd.read_csv(args.map_file, sep='\t', comment='~')
map_df.rename(columns={'#SampleID': 'sampleID'}, inplace=True)
map_df = map_df[[args.column, 'sampleID']]  # Adjust based on the column name you want to use

# Get unique treatment values from the specified column
treatments = map_df[args.column].unique()

# Ensure the output directory exists
os.makedirs(args.output_avg_abundance_directory, exist_ok=True)

# Transpose the ASV table to have sample IDs as rows (aligned with mapping file)
df_transposed = df.set_index('taxonomy').T  # Taxonomy becomes column names, sampleIDs are rows
df_transposed.reset_index(inplace=True)  # Sample IDs move to a column named "index"
df_transposed.rename(columns={'index': 'sampleID'}, inplace=True)  # Rename the index column to 'sampleID'

# Merge with the mapping file
merged_df = pd.merge(map_df, df_transposed, on='sampleID')  # Align treatments with ASV data

# Step 1: Identify top OTUs per treatment using Counter
top_otus_counter = Counter()

# Prepare a dictionary to store treatment average abundance DataFrames
treatment_avg_abundance_dict = {}

for treatment in treatments:
    # Filter for samples belonging to this treatment
    treatment_df = merged_df[merged_df[args.column] == treatment]
    # Drop non-numeric columns (like 'sampleID' and treatment column)
    treatment_abundance = treatment_df.drop(columns=['sampleID', args.column])
    # Calculate mean abundance for each taxonomy across samples in the treatment
    avg_abundance = treatment_abundance.mean(axis=0).reset_index()
    # Rename columns for clarity
    avg_abundance.columns = ['taxonomy', 'average_abundance']
    # Sort by abundance in descending order
    avg_abundance = avg_abundance.sort_values(by='average_abundance', ascending=False)
    # Save to dictionary
    treatment_avg_abundance_dict[treatment] = avg_abundance
    # Save to file
    output_file = os.path.join(args.output_avg_abundance_directory, f"{treatment}_avg_abundance.tsv")
    avg_abundance.to_csv(output_file, sep='\t', index=False)
    # Update top OTUs counter
    top_otus_counter.update(avg_abundance.head(args.num_taxa)['taxonomy'])

# Extract the top OTUs (those that appear most across treatments)
top_otus = [otu for otu, _ in top_otus_counter.most_common()]

# Step 2: Prepare the output DataFrame with predefined column names and types
columns = ['treatment'] + top_otus + ['other']
output_df = pd.DataFrame(columns=columns)

# Ensure correct column types
output_df = output_df.astype({col: float for col in columns if col != 'treatment'})
output_df['treatment'] = output_df['treatment'].astype(str)

# Step 3: Calculate average abundances for each treatment and "all"
for treatment in treatments:
    # Get average abundance for the treatment from the dictionary
    avg_abundance = treatment_avg_abundance_dict[treatment].set_index('taxonomy')['average_abundance']
    # Separate abundances for top OTUs and others
    row = {'treatment': treatment}
    row.update({otu: avg_abundance.get(otu, 0) for otu in top_otus})  # Top OTUs
    row['other'] = avg_abundance.drop(index=top_otus).sum()  # "Other" category
    # Append the row to the output DataFrame
    output_df = pd.concat([output_df, pd.DataFrame([row])], ignore_index=True)

# Step 4: Add the "all" row (average of treatment averages)
# Create a DataFrame to store treatment averages for each OTU
treatment_means = pd.DataFrame([
    treatment_avg_abundance_dict[treatment].set_index('taxonomy')['average_abundance']
    for treatment in treatments
]).fillna(0)  # Fill NaN for OTUs missing in some treatments

# Calculate the average of averages for the "all" row
overall_row = {'treatment': 'all'}
overall_row.update({otu: treatment_means[otu].mean() for otu in top_otus})
overall_row['other'] = treatment_means.drop(columns=top_otus).sum(axis=1).mean()  # Average of "other"

# Add the "all" row to the output DataFrame
output_df = pd.concat([pd.DataFrame([overall_row]), output_df], ignore_index=True)

# Step 5: Sort OTUs based on overall abundance (row "all")
all_row = output_df.loc[output_df['treatment'] == 'all']
sorted_otus = [otu for otu in all_row.iloc[0, 1:-1].sort_values(ascending=False).index]
sorted_columns = ['treatment'] + sorted_otus + ['other']

# Reorder columns in the output DataFrame
output_df = output_df[sorted_columns]

# Step 6: Convert values to proportions (percentage of the row sum)
for idx, row in output_df.iterrows():
    row_sum = row.drop(['treatment']).sum()  # Include "other" in the sum
    if row_sum > 0:
        output_df.iloc[idx, 1:] = (row.drop(['treatment']) / row_sum) * 100

# Step 7: Save the updated DataFrame to a file
output_df.to_csv(args.output_prism_file, sep='\t', index=False)
