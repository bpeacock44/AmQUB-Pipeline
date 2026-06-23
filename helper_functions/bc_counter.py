#!/usr/bin/env python3
import argparse
import csv
import sys
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
        - One row for every sample in the mapping file, including samples with a count of 0,
          in mapping-file order.
    """
    # Step 1: Extract barcode sequences and sample IDs (preserving mapping-file order)
    barcodes = {}        # barcode sequence -> sample ID
    sample_order = []    # list of (sample_ID, barcode) in mapping-file order
    with open(mapping_file, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')

        #make sure the mapping file actually has the columns we need
        required = ('#SampleID', 'BarcodeSequence')
        fields = reader.fieldnames or []
        missing = [col for col in required if col not in fields]
        if missing:
            sys.exit("** Error ** Mapping file ["+mapping_file+"] is missing required "
                     +"column(s): "+", ".join(missing))

        for row in reader:
            barcodes[row['BarcodeSequence']] = row['#SampleID']
            sample_order.append((row['#SampleID'], row['BarcodeSequence']))

    # Step 2: Count barcode occurrences (only the sequence line of each FASTQ record)
    barcode_counts = defaultdict(int)
    with open(fq_file, 'r') as f:
        for line_number, line in enumerate(f):
            #the barcode is on the sequence line (line 2 of each 4-line record)
            if line_number % 4 == 1:
                line = line.strip()
                if line in barcodes:
                    barcode_counts[line] += 1

    # Step 3: Write to output file (every sample, including zero counts, in mapping-file order)
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['sample_ID', 'barcode', 'read_count'])
        for sample_id, barcode in sample_order:
            writer.writerow([sample_id, barcode, barcode_counts[barcode]])

    print(f"Results saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Count barcode occurrences in a FASTQ file.')
    parser.add_argument('mapping_file', help='Path to the mapping file (e.g., JB117s_map.txt)')
    parser.add_argument('fq_file', help='Path to the FASTQ file (e.g., JB117s_BC.M0.fq)')
    parser.add_argument('output_file', help='Path to the output file (e.g., barcode_counts.tsv)')

    args = parser.parse_args()

    count_barcodes(args.mapping_file, args.fq_file, args.output_file)