#!/usr/bin/env Rscript

# Check if command line arguments are provided
if (length(commandArgs(trailingOnly = TRUE)) != 2) {
  stop("Usage: add_seqs_to_ASV_v2.R <input_text_ASV_table> <output_text_ASV_table>")
}

# Get command line arguments
otbl_fp <- commandArgs(trailingOnly = TRUE)[1]
outfa_fp <- commandArgs(trailingOnly = TRUE)[2]

# Load necessary functions
source('/home/bpeacock_ucr_edu/real_projects/PN94_singularity_of_microbiome_pipeline/targeted_microbiome_via_blast/helper_functions/pipeline_helper_functions.R')

# Add sequences to ASV table
add_sequences_to_ASV_table(otbl_fp, "rep_set/seqs_chimera_filtered_ASVs.fasta", outfa_fp)

