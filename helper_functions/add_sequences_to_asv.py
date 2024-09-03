#!/usr/bin/env python

import pandas as pd
import argparse
from Bio import SeqIO

def load_asv_table(tbl_fp):
    """
    Loads the ASV table from a file, skipping the first metadata line and setting the correct headers.
    """
    # Read the first few lines to determine the correct header
    with open(tbl_fp, 'r') as file:
        lines = file.readlines()

    # Skip the metadata line and use the second line as the header
    header_line = lines[1].strip()
    headers = header_line.split('\t')

    # Use the remaining lines to create the DataFrame
    tbl = pd.read_csv(tbl_fp, sep='\t', skiprows=2, header=None, names=headers, index_col=0)

    return tbl

def load_fasta_file(fasta_fp):
    """
    Loads a FASTA file into a dictionary with sequence IDs as keys.
    """
    fasta_seqs = {}
    for record in SeqIO.parse(fasta_fp, "fasta"):
        seq_id = record.id
        sequence = str(record.seq)
        fasta_seqs[seq_id] = sequence
    # Convert dictionary to DataFrame for easier handling
    fasta_df = pd.DataFrame(list(fasta_seqs.items()), columns=['ID', 'Sequence'])
    fasta_df.set_index('ID', inplace=True)
    return fasta_df

def add_sequences_to_asv_table(tbl_fp, fasta_fp, out_fp):
    """
    Adds sequences from a FASTA file to an ASV table and saves the result.
    """
    # Load the ASV table
    tbl = load_asv_table(tbl_fp)

    # Load the FASTA sequences
    fasta_seqs = load_fasta_file(fasta_fp)

    # Reorder sequences to match the ASV table row order
    fS = fasta_seqs.loc[tbl.index]

    # Check if row names agree
    if not all(tbl.index == fS.index):
        raise ValueError("Row names of ASV table and FASTA sequences do not match.")

    # Update sequence identifiers to remove any abundances, if necessary
    fS['ID'] = fS.index

    # Create combined sequence entries with the format: ">ID#Sequence"
    combined_seqs = fS.apply(lambda row: f">{row['ID']}#{row['Sequence']}", axis=1)

    # Add the combined sequences to the ASV table
    tbl['Sequence'] = combined_seqs.values

    # Save the modified ASV table
    tbl.to_csv(out_fp, sep='\t', index=True, header=True, quoting=0)

def main():
    parser = argparse.ArgumentParser(description="Add sequences from a FASTA file to an ASV table.")
    parser.add_argument('tbl_fp', type=str, help='Path to the ASV table file')
    parser.add_argument('fasta_fp', type=str, help='Path to the FASTA file')
    parser.add_argument('out_fp', type=str, help='Output file path for the resulting ASV table with sequences')
    
    args = parser.parse_args()
    
    add_sequences_to_asv_table(args.tbl_fp, args.fasta_fp, args.out_fp)

if __name__ == '__main__':
    main()
