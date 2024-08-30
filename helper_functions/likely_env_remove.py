#!/usr/bin/env python3

import sys
import time
from Bio import Entrez

def retrieve_taxonomy(query, output_file, email, max_attempts=3, retmax=100000):
    """
    Retrieves likely environmental sample taxonomic IDs.

    :param query: The search query for the taxonomy database.
    :param output_file: The path to the file where the results will be saved.
    :param email: The email address for NCBI Entrez requests.
    :param max_attempts: The maximum number of attempts to retrieve the data.
    :param retmax: Maximum number of IDs to retrieve in one request (default is 100,000).
    """
    Entrez.email = email
    attempt = 1
    success = False

    while attempt <= max_attempts:
        print(f"Attempt {attempt}: Retrieving taxonomy for {output_file}")
        try:
            # Use esearch to find taxonomy IDs with retmax to match the large output of bash script
            search_handle = Entrez.esearch(db="taxonomy", term=query, retmax=retmax)
            search_results = Entrez.read(search_handle)
            search_handle.close()

            # Get the list of taxonomy IDs
            tax_ids = search_results.get("IdList", [])

            # Write the taxonomy IDs to the output file
            with open(output_file, "w") as f:
                for tax_id in tax_ids:
                    f.write(f"{tax_id}\n")

            # Check if sufficient results were retrieved
            if tax_ids:
                success = True
                print(f"Successfully retrieved {len(tax_ids)} taxonomy IDs for {output_file}")
                break

        except Exception as e:
            print(f"Attempt {attempt} failed with error: {e}. Retrying in 5 seconds...")
            time.sleep(5)
            attempt += 1

    if not success:
        print(f"Failed to retrieve taxonomy for {output_file} after {max_attempts} attempts. "
              f"This may be due to requesting too many in quick succession. "
              f"Run part 4 again with the s flag to skip the BLAST.")
        sys.exit(1)

def filter_blast_output(env_accessions_file, blast_input_file, blast_output_file):
    """
    Filter BLAST output based on a list of taxonomy IDs from the env_accessions_file.
    
    :param env_accessions_file: The file containing taxids to filter out.
    :param blast_input_file: The BLAST input file.
    :param blast_output_file: The filtered BLAST output file.
    """
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
        print("Usage: python combined_script.py <email> <blast_input_file> <blast_output_file>")
        sys.exit(1)

    email = sys.argv[1]
    blast_input_file = sys.argv[2]
    blast_output_file = sys.argv[3]
    env_accessions_file = "likely_env_taxids_removed.txt"  # Output file for taxonomy retrieval

    # Define the query
    query = "\"environmental samples\"[subtree] OR \"Environmental Samples\"[subtree] OR " \
            "\"unclassified\"[subtree] OR \"Unclassified\"[subtree] OR " \
            "\"uncultured\"[subtree] OR \"Uncultured\"[subtree]"

    # Retrieve taxonomy and save it to file
    retrieve_taxonomy(query, env_accessions_file, email)

    # Step 2: Filter BLAST output
    filter_blast_output(env_accessions_file, blast_input_file, blast_output_file)
