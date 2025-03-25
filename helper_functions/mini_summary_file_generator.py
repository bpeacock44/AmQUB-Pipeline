#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse

def read_fasta(fa):
    # Read the fasta file and return a dataframe with 'ID' and 'sequence'
    ids = []
    sequences = []
    with open(fa, 'r') as f:
        header = None
        seq = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    ids.append(header)
                    sequences.append(seq)
                header = line[1:]  # Remove the '>' character
                seq = ''
            else:
                seq += line
        # Add the last sequence
        if header is not None:
            ids.append(header)
            sequences.append(seq)
    
    # Return a DataFrame with IDs and formatted sequences
    return pd.DataFrame({
        'ID': ids,
        'sequence': [f'>{header}#{seq}' for header, seq in zip(ids, sequences)]
    })

def process_data_and_write_summary(norm, fa, tax, tax_c, output_file):
    # Read the main data file without recognizing any comment characters
    df = pd.read_csv(norm, sep='\t', comment=None, header=1)
    df.rename(columns={df.columns[0]: 'ID'}, inplace=True)
        
    # Read primary taxonomy file
    tax_df = pd.read_csv(tax, sep='\t', comment=None, header=0)
    tax_df.columns = ['ID', 'taxonomy', 'blast_bitscore', 'blast_per_ID', 'blast_per_qcov']
    df = pd.merge(df, tax_df[['ID', 'taxonomy', 'blast_per_ID', 'blast_per_qcov']], on='ID', how='left')
    
    # Check if tax_c is provided and process it if available
    if tax_c and tax_c.lower() != 'null':
        # Read classifier taxonomy file if given
        tax_df_C = pd.read_csv(tax_c, sep='\t', comment=None, header=0)
        tax_df_C.columns = ['ID', 'classifier_taxonomy', 'classifier_confidence']
        df = pd.merge(df, tax_df_C, on='ID', how='left')
        # Include c_taxonomy and confidence in excluded columns
        excluded_columns = ['ID', 'taxonomy', 'blast_per_ID', 'blast_per_qcov', 'classifier_taxonomy', 'classifier_confidence', 'sequence']
    else:
        # Exclude c_taxonomy and confidence if tax_c is not provided
        excluded_columns = ['ID', 'taxonomy', 'blast_per_ID', 'blast_per_qcov', 'sequence']
    
    # Calculate mean abundance using only numeric columns
    numeric_cols = df.select_dtypes(include=[np.number]).columns.difference(excluded_columns)
    df['avg_abun'] = df[numeric_cols].mean(axis=1, skipna=True)

    # Remove numeric columns after 'avg_abun' has been created
    df.drop(columns=numeric_cols, inplace=True)

    # Read and process the fasta file
    fasta_df = read_fasta(fa)

    # Merge the fasta data with the main dataframe by 'ID'
    df = pd.merge(df, fasta_df[['ID', 'sequence']], on='ID', how='left')

    # Write dataframe to the specified output file
    df.to_csv(output_file, sep='\t', index=False, quoting=False)

def main():
    # Create an argument parser
    parser = argparse.ArgumentParser(description='Process data and write summary to file.')

    # Define command-line arguments
    parser.add_argument('norm', type=str, help='Path to the normalized data file')
    parser.add_argument('fa', type=str, help='Path to the fasta file')
    parser.add_argument('tax', type=str, help='Path to the blast taxonomy file')
    parser.add_argument('tax_c', type=str, nargs='?', default=None, help='Path to the classifier taxonomy file (optional, use "NULL" to skip)')
    parser.add_argument('output_file', type=str, help='Path to the output file')

    # Parse command-line arguments
    args = parser.parse_args()

    # Call the function with the arguments
    process_data_and_write_summary(args.norm, args.fa, args.tax, args.tax_c, args.output_file)

if __name__ == '__main__':
    main()
