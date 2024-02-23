
#This script adds abundance counts to a fasta file
#Usage:
# Rscript add_counts_to_fasta_seqs.R

source('/sw/paul_helper_scripts/pipeline_helper_functions.R')

otblfp="otu_table_01.txt"
fastafp="zotus.fa"
outfp="seqs_chimera_filtered_otus.fasta"

cat("Loading [otu_table_00.txt]\n")
add_counts_to_fasta_sequences(otblfp, fastafp, outfp)

#load otu table
cat("Saved [seqs_chimera_filtered_otus.fasta]\n")

