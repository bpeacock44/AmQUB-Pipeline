#!/bin/bash

set -euo pipefail

# Default values
QUERY=""
SUBJECT=""

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --query)
            QUERY="$2"
            shift 2
            ;;
        --subject)
            SUBJECT="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: $0 --query <query_file> --subject <subject_file>"
            exit 1
            ;;
    esac
done

# Check required inputs
if [[ -z "$QUERY" || -z "$SUBJECT" ]]; then
    echo "Error: --query and --subject are required."
    echo "Usage: $0 --query <query_file> --subject <subject_file>"
    exit 1
fi

OPTS="qseqid sseqid pident length mismatch evalue bitscore staxids stitle qcovs"
NUM_THREADS=$(nproc)

# Convert FASTQ to FASTA
convert_fastq_to_fasta() {
    local file="$1"
    local out="${file%.*}.fasta"
    seqkit fq2fa "$file" -o "$out"
    echo "$out"
}

# Convert if necessary
[[ "$QUERY" == *.fastq || "$QUERY" == *.fq ]] && echo "Converting query to fasta." && QUERY=$(convert_fastq_to_fasta "$QUERY")
[[ "$SUBJECT" == *.fastq || "$SUBJECT" == *.fq ]] && echo "Converting subject to fasta." && SUBJECT=$(convert_fastq_to_fasta "$SUBJECT")

# Base names for output
BASE1=$(basename "$QUERY" | sed 's/\.[^.]*$//')
BASE2=$(basename "$SUBJECT" | sed 's/\.[^.]*$//')

# Create BLAST DB directory
DB_DIR="blastdb_dir"
mkdir -p "$DB_DIR"
DB_PATH="$DB_DIR/blastdb_subject"
echo "Creating BLAST database for subject. Storing in ${DB_PATH}."

# Make BLAST database
makeblastdb -in "$SUBJECT" -dbtype nucl -out "$DB_PATH" -parse_seqids

# Output file
OUTFILE="blast_${BASE1}_against_${BASE2}.blastout"

echo "Running BLAST with ${NUM_THREADS} threads."

# Run BLAST
blastn -query "$QUERY" -db "$DB_PATH" -outfmt "6 $OPTS" -num_threads "$NUM_THREADS" > "$OUTFILE"

echo "BLAST complete. Output saved to: $OUTFILE"
