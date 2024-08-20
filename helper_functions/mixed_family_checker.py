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
                    queries[current_query] = []
                current_query = line.strip().split(":")[-1].strip().split()[0]
            elif not line.startswith("#") and current_query:
                fields = line.strip().split("\t")
                tax_id = fields[7].split(";")[0]
                queries[current_query].append(tax_id)
    if current_query:
        yield current_query, queries[current_query]

# Initialize tax_cache and failed_requests
tax_cache = {}
failed_requests = set()

def fetch_taxonomy(tax_ids, email):
    # Fetch multiple taxonomies at once
    tax_ids = [tid for tid in tax_ids if tid not in tax_cache and tid not in failed_requests]
    if not tax_ids:
        return

    Entrez.email = email
    try:
        handle = Entrez.efetch(db="taxonomy", id=",".join(tax_ids), retmode="xml")
        records = Entrez.read(handle)
        for record in records:
            tax_id = record['TaxId']
            lineage = record.get('LineageEx', [])
            family = next((entry['ScientificName'] for entry in lineage if entry['Rank'] == 'family'), None)
            tax_cache[tax_id] = family
    except HTTPError as err:
        if 500 <= err.code <= 599 or err.code == 400:
            print("Received error from server for tax_ids {}: {}".format(tax_ids, err))
            failed_requests.update(tax_ids)
        else:
            print("Error from server for tax_ids {}: {}".format(tax_ids, err))
    except Exception as e:
        print("Error fetching taxonomy for tax_ids {}: {}".format(tax_ids, e))

def process_identifiers(identifiers, output_file, email):
    to_fetch = set()
    for tax_ids in identifiers.values():
        to_fetch.update(tax_ids)
    
    fetch_taxonomy(to_fetch, email)  # Fetch taxonomy info in batch
    
    first_family = None
    with open(output_file, "a") as file:
        for identifier, tax_ids in identifiers.items():
            family_count = 0
            for tax_id in tax_ids:
                family = tax_cache.get(tax_id)
                if family:
                    if first_family is None:
                        first_family = family
                        family_count = 1
                    elif family != first_family:
                        file.write(identifier + "\n")
                        break
                    else:
                        family_count += 1
                    if family_count >= 10:
                        break

def main():
    parser = argparse.ArgumentParser(description='Process BLAST output file.')
    parser.add_argument('blastout_file', type=str, help='Path to BLAST output file')
    parser.add_argument('--email', type=str, help='Your email address for Entrez')
    args = parser.parse_args()
    output_file = "mixed_family_checker_out.txt"

    blast_output_file = args.blastout_file
    email = args.email
    for query, tax_ids in parse_blast_output(blast_output_file):
        process_identifiers({query: tax_ids}, output_file, email)

if __name__ == "__main__":
    main()
