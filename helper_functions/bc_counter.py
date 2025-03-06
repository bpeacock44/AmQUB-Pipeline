#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict

def count_barcodes(mapping_file, fq_file, output_file):
    """
    Count the occurrences of barcode sequences from a mapping file in a FASTQ file.

    Args:
        mapping_file (str): Path to the mapping file containing barcode sequences and sample IDs.
        fq_file (str): Path to the FASTQ file where barcode sequences will be searched.
        output_file (str): Path to the output file where results will be saved in tab-delimited format.
    
    The output file will contain:
        - Columns for sample_ID, barcode, and read_count.
    """
    # Step 1: Extract barcode sequences and sample IDs
    barcodes = {}
    with open(mapping_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            barcodes[row['BarcodeSequence']] = row['#SampleID']

    # Step 2: Count barcode occurrences
    barcode_counts = defaultdict(int)
    with open(fq_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line in barcodes:
                barcode_counts[line] += 1

    # Step 3: Write to output file
    with open(output_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['sample_ID', 'barcode', 'read_count'])
        for barcode, count in barcode_counts.items():
            writer.writerow([barcodes[barcode], barcode, count])

    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Count barcode occurrences in a FASTQ file.')
    parser.add_argument('mapping_file', help='Path to the mapping file (e.g., JB117s_map.txt)')
    parser.add_argument('fq_file', help='Path to the FASTQ file (e.g., JB117s_BC.M0.fq)')
    parser.add_argument('output_file', help='Path to the output file (e.g., barcode_counts.tsv)')
    
    args = parser.parse_args()
    
    count_barcodes(args.mapping_file, args.fq_file, args.output_file)
