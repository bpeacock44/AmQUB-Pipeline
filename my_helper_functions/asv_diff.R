library(readr)
library(dplyr)

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 6) {
  stop("Usage: Rscript ASV_diff.R diff_path map_path ASV_path column_of_interest, var_a, var_b")
}

# Assign arguments to variables
diff_path <- args[1]
map_path <- args[2]
ASV_path <- args[3]
column_of_interest <- args[4]
var_a <- args[5]
var_b <- args[6]

# Read in files
diff <- read.table(diff_path, header = TRUE, sep = "\t", comment.char = "~")
map <- read_tsv(map_path, comment = "~") %>%
  rename(sampleID = `#SampleID`)
ASV <- read.table(ASV_path, header = TRUE, sep = "\t", comment.char = "~", skip = 1)
ASV_ids <- ASV$"X.ASV.ID"

# Extract sample IDs corresponding to variable of choice
a_samples <- map$sampleID[map[[column_of_interest]] == var_a]
b_samples <- map$sampleID[map[[column_of_interest]] == var_b]

# Subset ASV table based on samples and calculate row means
rm_a <- rowMeans(ASV[, a_samples, drop = FALSE])
rm_b <- rowMeans(ASV[, b_samples, drop = FALSE])

# Merge into final dataframe
col_a_name <- paste0(var_a,"_avg")
col_b_name <- paste0(var_b,"_avg")
avgs <- cbind(ASV_ids,rm_a,rm_b)
colnames(avgs) <- c("ASV",col_a_name,col_b_name)
final <- merge(diff,avgs,by="ASV",all.x=TRUE)

# Sort and save
final <- final %>% arrange(ASV)
write.table(final, file = "differential_abundance_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)