#!/usr/bin/env python3
import sys

def filter_blast_output(env_accessions_file, blast_input_file, blast_output_file):
    # Load the list of taxids to filter out
    with open(env_accessions_file, "r") as f:
        filter_set = set(line.strip() for line in f if line.strip())

    # Process the blast output file
    with open(blast_input_file, "r") as infile, open(blast_output_file, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                # Preserve comment lines
                outfile.write(line)
                continue
            
            columns = line.strip().split("\t")
            if len(columns) >= 8:
                # Split the 8th column by ';' and check if any match the filter set
                accessions = columns[7].split(";")
                if not any(accession in filter_set for accession in accessions):
                    outfile.write(line)
            else:
                # In case the line has fewer than 8 columns, write it as is
                outfile.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python filter_blast_output.py <env_accessions_file> <blast_input_file> <blast_output_file>")
        sys.exit(1)

    env_accessions_file = sys.argv[1]
    blast_input_file = sys.argv[2]
    blast_output_file = sys.argv[3]

    filter_blast_output(env_accessions_file, blast_input_file, blast_output_file)

