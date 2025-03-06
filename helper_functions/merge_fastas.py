#!/usr/bin/env python3
from Bio import SeqIO
import sys

"""
This script removes exact duplicate sequences from a FASTA file while preserving reverse complements as separate entries.
It merges all sequence IDs for identical sequences into a single header.

Usage:
    python deduplicate_fasta.py <input_file>

Arguments:
    input_file - Path to the input FASTA file.

Output:
    A new FASTA file named <input_file>_unique.fa will be created, containing only unique sequences.

How it works:
    - Exact duplicates are merged (case-insensitive), but reverse complements are retained as distinct sequences.
    - All IDs for identical sequences are combined into a single FASTA header.
"""

# Input and output file paths
input_file = sys.argv[1]
output_file = f"{input_file.rsplit('.', 1)[0]}_unique.fa"

def normalize_sequence(sequence):
    return "".join(str(sequence).upper().split())

def deduplicate_sequences(input_file, output_file):
    sequence_dict = {}  # Dictionary to store unique sequences and their IDs

    # Parse the input FASTA file
    for record in SeqIO.parse(input_file, "fasta"):
        normalized_sequence = normalize_sequence(record.seq)

        if normalized_sequence in sequence_dict:
            # Merge IDs if the sequence already exists
            sequence_dict[normalized_sequence].append(record.id)
        else:
            # Add new unique sequence and its ID
            sequence_dict[normalized_sequence] = [record.id]

    # Write unique sequences with merged IDs to the output file
    with open(output_file, "w") as out_fasta:
        for sequence, ids in sequence_dict.items():
            merged_ids = ",".join(ids)
            out_fasta.write(f">{merged_ids}\n{sequence}\n")

    print(f"Processed file saved as {output_file}")

# Run the deduplication process
deduplicate_sequences(input_file, output_file)
