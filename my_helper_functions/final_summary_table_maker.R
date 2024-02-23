# Usage: Rscript final_summary_table_maker.R
# This script creates an excel file of the OTU table with detailed information about 
# taxonomic assignments and other relevant information for OTUs. 
# Files required to run this script are defined below. 

# Load required libraries
library(dplyr)
library(writexl)

# Define file paths
file_path <- "otu_table_03_add_seqs_norm.txt"
file_path3 <- "rep_set/assgntax/seqs_chimera_filtered_tax_assignments.txt"
file_path2 <- "rep_set/assgntax/nf_seqs_chimera_filtered_tax_assignments.txt"
file_path4 <- "otu_table_03_add_seqs.txt"
file_path_top10 <- "rep_set/top_10_family_checker_out.txt"

# Read the main data file
df <- read.table(file_path, sep = "\t", comment.char = "", header = TRUE)
colnames(df)[1] <- "OTU_ID"
df$taxonomy <- NULL

# Read primary taxonomy file
tax_df <- read.table(file_path3, sep = "\t", comment.char = "", header = TRUE)
colnames(tax_df) <- c("OTU_ID", "taxonomy", "bitscore", "per_ID", "per_qcov", "size")
df <- merge(df, tax_df[c(1:2,4:5)], by = "OTU_ID", all.x = TRUE)

# Read alternate taxonomy files
tax_df_alt <- read.table(file_path2, sep = "\t", comment.char = "", header = TRUE)
colnames(tax_df_alt) <- c("OTU_ID", "nf_taxonomy", "nf_bitscore", "nf_per_ID", "nf_per_qcov", "size")
df <- merge(df, tax_df_alt[c(1:2,4:5)], by = "OTU_ID", all.x = TRUE)

# Calculate mean abundance
excluded_columns <- c("OTU_ID", "nf_taxonomy", "nf_per_ID", "nf_per_qcov", "sequence", "taxonomy", "per_ID", "per_qcov")
df$avg_abun <- rowMeans(df[, !names(df) %in% excluded_columns], na.rm = TRUE)

# Flag OTUs with mixed families in top 10
top10_df <- read.table(file_path_top10)
df <- mutate(df, mixed_fam_top_10 = "no")
top10_identifiers <- top10_df[[1]]
matches <- df$OTU_ID %in% top10_identifiers
df$mixed_fam_top_10[matches] <- "yes"

# Reorder columns
df <- df[, c(1, (length(df) - 8):length(df), 2:(length(df) - 9))]

# Add row with total raw counts
raw <- read.table(file_path4, sep = "\t", comment.char = "", header = TRUE)
raw[,1] <- NULL
raw$taxonomy <- NULL
raw$sequence <- NULL
sums <- colSums(raw)
nd_columns <- rep("nd", 10)
new_vector <- c(nd_columns, sums)
names(new_vector)[1:10] <- colnames(df)[1:10]
df <- rbind(new_vector, df)

# Write dataframe to an Excel file
write_xlsx(df, "otu_summary_table.xlsx")
