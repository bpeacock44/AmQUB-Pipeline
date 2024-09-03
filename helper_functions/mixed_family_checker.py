#!/usr/bin/env python3

from collections import defaultdict
from Bio import Entrez
from urllib.error import HTTPError
import time
import argparse

def parse_blast_output(file_path):
    queries = defaultdict(list)
    current_query = None
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("# Query:"):
                if current_query:
                    yield current_query, queries[current_query]
                current_query = line.strip().split(":")[-1].strip().split()[0]
                queries[current_query] = []
            elif not line.startswith("#") and current_query:
                fields = line.strip().split("\t")
                tax_id = fields[7].split(";")[0]
                if len(queries[current_query]) < 10:
                    queries[current_query].append(tax_id)
    if current_query:
        yield current_query, queries[current_query]

def fetch_taxonomy(tax_ids, email, batch_size=200):
    tax_ids = [tid for tid in tax_ids if tid]
    Entrez.email = email
    tax_cache = {}
    attempts = 3

    for i in range(0, len(tax_ids), batch_size):
        batch = tax_ids[i:i + batch_size]
        while attempts > 0:
            try:
                handle = Entrez.efetch(db="taxonomy", id=",".join(batch), retmode="xml")
                records = Entrez.read(handle)
                for record in records:
                    tax_id = record['TaxId']
                    lineage = record.get('LineageEx', [])
                    family = next((entry['ScientificName'] for entry in lineage if entry['Rank'] == 'family'), None)
                    tax_cache[tax_id] = family
                break
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    attempts -= 1
                    time.sleep(2 ** (3 - attempts))
                else:
                    break
            except Exception as e:
                break
    return tax_cache

def process_identifiers(identifiers, tax_cache, output_file):
    with open(output_file, "w") as file:
        for identifier, tax_ids in identifiers.items():
            families = set()
            for tax_id in tax_ids:
                family = tax_cache.get(tax_id)
                if family:
                    families.add(family)
                if len(families) > 1:
                    file.write(identifier + "\n")
                    break

def main():
    parser = argparse.ArgumentParser(description='Process BLAST output file.')
    parser.add_argument('blastout_file', type=str, help='Path to BLAST output file')
    parser.add_argument('--email', type=str, help='Your email address for Entrez')
    parser.add_argument('--output', type=str, default="mixed_family_checker_out.txt", help='Path to output file (default: mixed_family_checker_out.txt)')
    args = parser.parse_args()

    blast_output_file = args.blastout_file
    email = args.email
    output_file = args.output

    # Step 1: Parse BLAST output and collect Tax IDs
    identifiers = {}
    for query, tax_ids in parse_blast_output(blast_output_file):
        identifiers[query] = tax_ids

    # Step 2: Fetch taxonomy data
    all_tax_ids = {tax_id for tax_ids in identifiers.values() for tax_id in tax_ids}
    tax_cache = fetch_taxonomy(all_tax_ids, email)

    # Step 3: Process identifiers to identify queries with multiple families
    process_identifiers(identifiers, tax_cache, output_file)

if __name__ == "__main__":
    main()
