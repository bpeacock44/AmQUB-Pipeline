
#This script adds abundance counts to a fasta file
#Usage:
# Rscript add_counts_to_fasta_seqs.R

source('/sw/paul_helper_scripts/pipeline_helper_functions.R')

otblfp="ASV_table_01.txt"
fastafp="ASVs.fa"
outfp="seqs_chimera_filtered_ASVs.fasta"

cat("Loading [ASV_table_00.txt]\n")
add_counts_to_fasta_sequences(otblfp, fastafp, outfp)

#load ASV table
cat("Saved [seqs_chimera_filtered_ASVs.fasta]\n")

