#!/usr/bin/env python

import argparse
from Bio import SeqIO

"""
This script updates the headers of sequences in a FASTA file using counts from an ASV table. 
It adds the total count for each ASV to the FASTA header, sorts the sequences by count in descending order, 
and optionally modifies the headers with a custom prefix and counter.

Usage:
    python update_fasta_headers.py <table> <fasta_file> <output_file> [--modify-headers] [--typ <prefix>]

Arguments:
    table      - Path to the table file containing counts for each taxonomic unit.
    fasta_file     - Path to the input FASTA file with sequence records.
    output_file    - Path where the updated FASTA file will be saved.
    --modify-headers  (optional) - If set, modifies the headers with a custom prefix and counter.
    --typ <prefix>    (optional) - Type to append to the header prefix (required with --modify-headers).

Example:
    python add_counts_to_fasta.py asv_table_00.txt otus.fasta updated_otus.fasta --modify-headers --typ Otu
"""

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

def update_fasta_headers(fasta_fp, asv_counts_fp, updated_fasta_fp, modify_headers=False, typ=None):
    asv_counts = read_asv_counts(asv_counts_fp)
    
    # Read all records first
    records = list(SeqIO.parse(fasta_fp, "fasta"))

    # Add counts to headers
    for record in records:
        asv_id = record.id
        if asv_id in asv_counts:
            total_count = asv_counts[asv_id]
            record.id = f"{asv_id}_{total_count}"  # Add total count to the header

    # Sort the records based on the total count (descending order)
    records.sort(key=lambda record: asv_counts.get(record.id.rsplit('_', 1)[0], 0), reverse=True)

    with open(updated_fasta_fp, 'w') as updated_fasta_file:
        counter = 1
        for record in records:
            new_id = record.id

            # Modify the header if the flag is set
            if modify_headers:
                if typ:
                    # Add m{typ}{counter} prefix and remove the last _number part
                    new_id = f"m{typ}{counter}_{new_id.rsplit('_', 1)[0]}"
                    counter += 1  # Increment the counter

            # Strip the description part and update the record
            record.id = new_id
            record.description = ""  # Remove the description part

            SeqIO.write(record, updated_fasta_file, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Add total counts to FASTA headers and optionally modify headers.")
    parser.add_argument('table', type=str, help='Path to the table file')
    parser.add_argument('fasta_file', type=str, help='Path to the FASTA file')
    parser.add_argument('output_file', type=str, help='Path to the output FASTA file')
    parser.add_argument('--modify-headers', action='store_true', help="Modify headers with a prefix and number")
    parser.add_argument('--typ', type=str, help="Type to append to the header prefix (required with --modify-headers)", default=None)

    args = parser.parse_args()

    if args.modify_headers and not args.typ:
        print("Error: --typ is required when using --modify-headers")
        exit(1)

    update_fasta_headers(args.fasta_file, args.table, args.output_file, modify_headers=args.modify_headers, typ=args.typ)

if __name__ == '__main__':
    main()
