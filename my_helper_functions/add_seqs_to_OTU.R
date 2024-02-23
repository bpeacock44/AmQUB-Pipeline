#!/usr/bin/env Rscript

# Check if command line arguments are provided
if (length(commandArgs(trailingOnly = TRUE)) != 2) {
  stop("Usage: add_seqs_to_OTU_v2.R <input_text_otu_table> <output_text_otu_table>")
}

# Get command line arguments
otbl_fp <- commandArgs(trailingOnly = TRUE)[1]
outfa_fp <- commandArgs(trailingOnly = TRUE)[2]

# Load necessary functions
source('/sw/paul_helper_scripts/pipeline_helper_functions.R')

# Add sequences to otu table
add_sequences_to_otu_table(otbl_fp, "rep_set/seqs_chimera_filtered_otus.fasta", outfa_fp)

