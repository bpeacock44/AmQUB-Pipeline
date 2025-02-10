#!/usr/bin/env python

import pandas as pd
import sys

"""
This script merges two OTU tables that contain some of the same samples. It ensures that all columns 
from both tables are retained in the merged output, filling missing values with 0 for columns 
that don't exist in one of the input tables. The output is saved as a new tab-delimited file.

Usage:
    python align_otu.py <input_file1> <input_file2> <output_file>

Arguments:
    input_file1  - Path to the first OTU table
    input_file2  - Path to the second OTU table
    output_file  - Path where the merged table will be saved

Example:
    python align_otu.py pre-existing_otu_table_00.txt unbinned_otu_table_00.txt merged_otu_table.txt
"""

def main(file1, file2, output_file):
    # Load both OTU tables as dataframes
    otu1 = pd.read_csv(file1, sep="\t")
    otu2 = pd.read_csv(file2, sep="\t")

    # Get the combined set of columns from both dataframes (including "#OTU ID")
    all_columns = sorted(set(otu1.columns) | set(otu2.columns))

    # Reindex both dataframes to align with the combined columns, filling missing values with 0
    otu1 = otu1.reindex(columns=all_columns, fill_value=0)
    otu2 = otu2.reindex(columns=all_columns, fill_value=0)

    # Concatenate the two dataframes
    merged_otu = pd.concat([otu1, otu2], axis=0)

    # Save the merged table to the output file
    merged_otu.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python align_otu.py <input_file1> <input_file2> <output_file>")
        sys.exit(1)

    input_file1 = sys.argv[1]
    input_file2 = sys.argv[2]
    output_file = sys.argv[3]

    main(input_file1, input_file2, output_file)