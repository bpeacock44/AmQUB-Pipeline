#!/usr/bin/env python

import argparse
from Bio import SeqIO

def read_asv_counts(asv_table_fp):
    asv_counts = {}
    with open(asv_table_fp, 'r') as file:
        lines = file.readlines()
        header = lines[0].strip().split()
        for line in lines[1:]:
            parts = line.strip().split()
            asv_id = parts[0]
            counts = map(int, parts[1:])
            total_count = sum(counts)
            asv_counts[asv_id] = total_count
    return asv_counts

def update_fasta_headers(fasta_fp, asv_counts_fp, updated_fasta_fp):
    asv_counts = read_asv_counts(asv_counts_fp)
    
    with open(fasta_fp, 'r') as fasta_file, open(updated_fasta_fp, 'w') as updated_fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            asv_id = record.id
            if asv_id in asv_counts:
                total_count = asv_counts[asv_id]
                new_id = f"{asv_id}_{total_count}"
                record.id = new_id
                record.description = ""  # Remove the description part
            SeqIO.write(record, updated_fasta_file, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Add total counts to FASTA headers.")
    parser.add_argument('asv_table', type=str, help='Path to the ASV table file')
    parser.add_argument('fasta_file', type=str, help='Path to the FASTA file')
    parser.add_argument('output_file', type=str, help='Path to the output FASTA file')
    
    args = parser.parse_args()
    
    update_fasta_headers(args.fasta_file, args.asv_table, args.output_file)

if __name__ == '__main__':
    main()
