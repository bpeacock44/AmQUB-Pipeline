#!/usr/bin/env python3
import os
import sys

def load_taxonomy_file(file_path, min_fields):
    """Load a taxonomy file into a dictionary."""
    taxonomy_dict = {}
    with open(file_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= min_fields:
                asv_id = fields[0]
                taxonomy = fields[1]
                taxonomy_dict[asv_id] = taxonomy
    return taxonomy_dict

def main(blastout_file, taxonomy_file, blast_assignments_file, solid_hovi_file, putative_hovi_file):
    # Step 1: Extract ASV IDs from BLAST output (first column only)
    with open(blastout_file, 'r') as f:
        asv_ids = set(line.strip().split('\t')[0] for line in f if line.strip())
    # Step 2: Load taxonomy data into dictionaries
    taxonomy_dict = load_taxonomy_file(taxonomy_file, min_fields=3)
    blast_dict = load_taxonomy_file(blast_assignments_file, min_fields=5)
    # Step 3: Initialize sets for output
    solid_hovi_asvs = set()
    putative_hovi_asvs = set()
    # Define partial and exact match lists
    partial_keywords = [
        "g__unclassified_Eukaryota", "g__unclassified_Ascomycota", 
        "g__Ascomycota_gen_Incertae_sedis", "g__Fungi_phy_Incertae_sedis"
    ]
    exact_matches = [
        "k__Eukaryota;p__Ascomycota", "k__Eukaryota", 
        "k__Fungi", "k__Fungi;p__Ascomycota"
    ]
    # Orbiliales-related keywords
    orbiliales_keywords = ["Hyalorbilia", "Orbiliales", "Orbiliaceae"]

    # Step 4: Apply filtering criteria
    def matches_exact_or_partial(taxonomy):
        """Check if a taxonomy matches either exact or partial criteria."""
        return (
            taxonomy in exact_matches or
            any(keyword in taxonomy for keyword in partial_keywords)
        )

    for asv_id in asv_ids:
        taxonomy1 = taxonomy_dict.get(asv_id, "")
        taxonomy2 = blast_dict.get(asv_id, "")
        # Check for "oviparasitica"
        if "oviparasitica" in taxonomy1 or "oviparasitica" in taxonomy2:
            solid_hovi_asvs.add(asv_id)
        # Check for Orbiliales-related terms
        if any(keyword in taxonomy1 for keyword in orbiliales_keywords) or \
           any(keyword in taxonomy2 for keyword in orbiliales_keywords):
            putative_hovi_asvs.add(asv_id)
            continue  # If Orbiliales criteria met, skip further checks for this ASV
        # Check if both taxonomies match either exact or partial criteria
        if matches_exact_or_partial(taxonomy1) and matches_exact_or_partial(taxonomy2):
            putative_hovi_asvs.add(asv_id)

    # Step 5: Write results to output files
    with open(solid_hovi_file, 'w') as f:
        for asv in sorted(solid_hovi_asvs):
            f.write(f"{asv}\n")
    with open(putative_hovi_file, 'w') as f:
        for asv in sorted(putative_hovi_asvs):
            f.write(f"{asv}\n")
    print("Filtering complete. Results saved to output files.")

# Ensure the script is run with the correct number of arguments
if __name__ == "__main__":
    if len(sys.argv) != 6:
        print(f"Usage: {sys.argv[0]} <blastout_file> <taxonomy_file> <blast_assignments_file> <solid_hovi_file> <putative_hovi_file>")
        sys.exit(1)

    # Assign command-line arguments to variables
    blastout_file = sys.argv[1]
    taxonomy_file = sys.argv[2]
    blast_assignments_file = sys.argv[3]
    solid_hovi_file = sys.argv[4]
    putative_hovi_file = sys.argv[5]

    # Call the main function
    main(blastout_file, taxonomy_file, blast_assignments_file, solid_hovi_file, putative_hovi_file)
