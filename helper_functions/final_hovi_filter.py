#!/usr/bin/env python3 
import re
import argparse
from collections import defaultdict

# Define search/exclude terms
TERMS = {"Hyalorbilia", "Brachyphoris", "Dactylella"}

def extract_genus(stitle):
    """Extract the genus (first word) from the stitle and strip special characters."""
    genus = stitle.split()[0] if stitle else ""
    # Remove any non-alphabetic characters from the genus
    genus = re.sub(r'[^a-zA-Z]', '', genus)
    return genus

def process_blast_output(blast_file):
    """Process the BLAST output to find matching qseqids."""
    qseqid_hits = defaultdict(list)
    filtered_qseqids = set()

    # Read the BLAST output
    with open(blast_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue

            qseqid, stitle = fields[0], fields[-2]
            qseqid_hits[qseqid].append(stitle)
    
    print("Parsed BLAST output. Processing qseqids...")

    for qseqid, hits in qseqid_hits.items():
        print(f"\nProcessing qseqid: {qseqid}")
        print(f"Total hits: {len(hits)}")

        # Get the first 50 hits
        top_hits = hits[:50]

        # Condition 1: Search terms in stitle
        if any(any(term in stitle for term in TERMS) for stitle in top_hits):
            print(f"qseqid {qseqid} included due to search terms in titles.")
            filtered_qseqids.add(qseqid)
            continue

        # Condition 2: Fewer than 10 hits
        if len(top_hits) < 10:
            print(f"qseqid {qseqid} included due to fewer than 10 hits.")
            filtered_qseqids.add(qseqid)
            continue

        # Condition 3: Fewer than 6 hits in top 10 share the same genus
        genus_counts = defaultdict(int)
        for stitle in top_hits[:10]:
            genus = extract_genus(stitle)
            genus_counts[genus] += 1

        print(f"Genus counts in top 10 hits: {dict(genus_counts)}")

        predominant_genus, predominant_count = max(genus_counts.items(), key=lambda x: x[1], default=(None, 0))
        print(f"Predominant genus: {predominant_genus}, Count: {predominant_count}")

        if predominant_count < 6 and predominant_genus not in TERMS:
            print(f"qseqid {qseqid} included due to genus count criteria.")
            filtered_qseqids.add(qseqid)

    return filtered_qseqids

def save_qseqids(qseqids, output_file):
    """Save the filtered qseqids to a file."""
    print(f"Saving {len(qseqids)} filtered qseqids to {output_file}.")
    with open(output_file, "w") as f:
        for qseqid in sorted(qseqids):
            f.write(f"{qseqid}\n")

def main():
    """Main function to process the BLAST output and filter qseqids."""
    parser = argparse.ArgumentParser(description="Filter qseqids from BLAST output based on criteria.")
    parser.add_argument("blast_file", help="Path to the BLAST output file.")
    parser.add_argument("output_file", help="Path to save the filtered qseqids.")
    args = parser.parse_args()

    filtered_qseqids = process_blast_output(args.blast_file)
    save_qseqids(filtered_qseqids, args.output_file)
    print(f"Filtered qseqids saved to {args.output_file}")

if __name__ == "__main__":
    main()
