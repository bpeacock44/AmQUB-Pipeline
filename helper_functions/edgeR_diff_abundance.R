#!/usr/bin/env Rscript

library(optparse)
library(edgeR)

# Define command line options
option_list <- list(
  make_option(c("-t", "--asv_table"), type="character", default=NULL,
              help="Path to the ASV table", metavar="character"),
  make_option(c("-m", "--mapping_file"), type="character", default=NULL,
              help="Path to the mapping file", metavar="character"),
  make_option(c("-d", "--diff_col"), type="character", default=NULL,
              help="Column name in mapping file to use for differential expression", metavar="character"),
  make_option(c("-v1", "--var_1"), type="character", default=NULL,
              help="First variable level", metavar="character"),
  make_option(c("-v2", "--var_2"), type="character", default=NULL,
              help="Second variable level", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default="output_results.txt",
              help="Output file name [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Function to run edgeR analysis
run_edgeR <- function(asv_table, mapping_file, diff_col, var_1, var_2, output_file) {
  
  # Load ASV table
  asv_data <- read.table(asv_table, header=TRUE, comment.char="", skip=1, sep="\t", row.names=1)
  
  # Load mapping file
  map_data <- read.table(mapping_file, header=TRUE, comment.char="", row.names=1, sep="\t")
  
  # Subset mapping file based on factors var_1 and var_2
  subset_map <- map_data[map_data[[diff_col]] %in% c(var_1, var_2), ]
  
  # Subset ASV table based on samples in the subsetted mapping file
  subset_asv <- asv_data[, colnames(asv_data) %in% rownames(subset_map)]
  
  # Ensure that the diff_col is a factor and has levels var_1 and var_2
  subset_map[[diff_col]] <- factor(subset_map[[diff_col]], levels = c(var_1, var_2))
  
  # Create DGEList object
  dge <- DGEList(counts = subset_asv, group = subset_map[[diff_col]])
  
  # Perform edgeR analysis
  design <- model.matrix(~0 + subset_map[[diff_col]])
  colnames(design) <- levels(subset_map[[diff_col]])
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  contrast <- makeContrasts(contrasts = paste(var_1, var_2, sep = "-"), levels=design)
  res <- glmLRT(fit, contrast=contrast)
  
  # Write all results to output file
  write.table(topTags(res, n = nrow(res$table))$table, file=output_file, sep="\t", quote=FALSE, row.names=TRUE)
  
  # Print summary
  summary(res)
}

# Check if mandatory arguments are provided
if (is.null(opt$asv_table) || is.null(opt$mapping_file) || is.null(opt$diff_col) || is.null(opt$var_1) || is.null(opt$var_2)) {
  print_help(opt_parser)
  stop("Missing argument(s)", call.=FALSE)
}

# Run edgeR analysis with provided options
run_edgeR(opt$asv_table, opt$mapping_file, opt$diff_col, opt$var_1, opt$var_2, opt$output_file)
