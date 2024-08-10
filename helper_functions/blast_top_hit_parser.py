#!/usr/bin/env python3

"""
blast_top_hit_parser.py
Loads a BLAST output file and assigns taxonomy to the ASVs/query sequences.
BLAST output must have the following outfmt:
OPTS="qseqid sseqid pident length mismatch evalue bitscore staxids stitle qcovs"
-outfmt "7 $OPTS"

We recommend setting max sequences to 5000. 

Usage:
    blast_top_hit_parser.py -i <file> -o <file>
    blast_top_hit_parser.py -h
Arguments:
    -i <file>     filepath to BLAST output
    -o <output>   parsed BLAST output file
Related scripts:
1) assign_LCA_via_blast.py (uses the output from blast_top_hit_parser.py to generate taxonomic assignments)
"""

import re
import argparse
from collections import defaultdict
import os

def parse_blast_output(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    query_hits = defaultdict(lambda: defaultdict(lambda: {'count': 0, 'parts': None}))
    current_query = None
    max_bitscore = 0

    for line in lines:
        line = line.strip()

        if line.startswith('# Query:'):
            match = re.match(r'# Query: (\S+)(?: (\d+))?', line)
            if match:
                query_id = match.group(1)
                query_number = match.group(2) if match.group(2) else ''
                if current_query:
                    for tax_id in query_hits[current_query]:
                        if query_hits[current_query][tax_id]['parts'] is None:
                            continue
                        if float(query_hits[current_query][tax_id]['parts'][6]) == max_bitscore:
                            query_hits[current_query][tax_id]['count'] = 1

                current_query = f"{query_id}_{query_number}"
                max_bitscore = 0
        elif not line.startswith('#'):
            parts = line.split('\t')
            if len(parts) < 10:
                continue

            query_id = parts[0]
            bitscore = float(parts[6])
            tax_ids = parts[7]
            query_coverage = parts[9]

            if ";" in tax_ids:
                continue  # Skip hits with multiple tax IDs

            tax_id = tax_ids.split(';')[-1]

            if query_id == current_query.split('_')[0]:
                if bitscore > max_bitscore:
                    max_bitscore = bitscore
                    query_hits[current_query] = defaultdict(lambda: {'count': 0, 'parts': None})
                if bitscore == max_bitscore:
                    if tax_id not in query_hits[current_query]:
                        query_hits[current_query][tax_id] = {'count': 1, 'parts': parts}
                    else:
                        query_hits[current_query][tax_id]['count'] += 1

    if current_query:
        for tax_id in query_hits[current_query]:
            if query_hits[current_query][tax_id]['parts'] is None:
                continue
            if float(query_hits[current_query][tax_id]['parts'][6]) == max_bitscore:
                query_hits[current_query][tax_id]['count'] = 1

    return query_hits

def format_output(query_hits):
    output_lines = []

    for query_id, hits in query_hits.items():
        hit_groups = []

        for tax_id, hit_data in hits.items():
            bitscore = hit_data['parts'][6]
            percent_id = hit_data['parts'][2]
            query_coverage = hit_data['parts'][9]
            tax_id_count = hit_data['count']
            hit_groups.append(f"({bitscore}|{percent_id}|{tax_id}|{tax_id_count}|{query_coverage})")

        formatted_hits = ' '.join(hit_groups)
        output_lines.append(f"{query_id}\t[{formatted_hits}]")

    return '\n'.join(output_lines)

def main():
    parser = argparse.ArgumentParser(description='Process BLAST output and format it.')
    parser.add_argument('-i', '--input', required=True, type=str, help='Path to the input BLAST output file')
    parser.add_argument('-o', '--output', required=True, type=str, help='Path to the output file')
    
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"Error: Input file '{args.input}' does not exist.")
        return

    query_hits = parse_blast_output(args.input)
    formatted_output = format_output(query_hits)

    with open(args.output, 'w') as file:
        file.write(formatted_output)

    print(f"Formatted output saved to '{args.output}'")

if __name__ == "__main__":
    main()
