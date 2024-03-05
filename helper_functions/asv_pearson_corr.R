library(readr)
library(dplyr)

# Check if there are enough command-line arguments
if (length(commandArgs(trailingOnly = TRUE)) < 3) {
  stop("Usage: Rscript asv_pearson_corr.R <file_path> <map_file> <rating_column>")
}

# Get command-line arguments
file_path <- commandArgs(trailingOnly = TRUE)[1]
map_file <- commandArgs(trailingOnly = TRUE)[2]
rating_column <- commandArgs(trailingOnly = TRUE)[3]

# read in data file and store ASV names
df <- read.table(file_path, header = TRUE, sep = "\t", comment.char = "~", skip = 1)
ASV_ids <- df$"X.ASV.ID"
df <- within(df, { "X.ASV.ID" <- NULL; "taxonomy" <- NULL })

# Read map file using readr functions
map <- read_tsv(map_file, comment = "~") %>%
  rename(sampleID = `#SampleID`) %>%
  select({{ rating_column }}, sampleID)

df2 <- df %>% filter_all(any_vars(!is.na(.)))
df3 <- data.frame(t(df2), stringsAsFactors = FALSE )
df3$sampleID <- rownames(df3)
df4 <- merge(map,df3,by="sampleID",all.y=TRUE)
df4 <- na.omit(df4)
avg_values <- colMeans(df4[,3:ncol(df4)])

# Apply correlation and p-value functions
results <- lapply(df4[, 3:ncol(df4)], function(col) {
  cor_val <- cor(col, df4$Rating)
  p_val <- cor.test(col, df4$Rating)$p.value
  fdr_val <- p.adjust(p_val, method = "fdr", n = ncol(df4) - 2)
  data.frame(cor_coef = cor_val, p_value = p_val, FDR = fdr_val)
})

# Combine results into a data frame
corres_data <- do.call(rbind, results)
final <- cbind(ASV_ids,corres_data,avg_values)
colnames(final) <- c("ASV","cor_coef","p_value","FDR","avg_count")

final <- final %>% arrange(ASV)
write.table(final, file = "correlation_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)