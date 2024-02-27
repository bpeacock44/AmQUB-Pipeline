#!/usr/bin/env Rscript
setwd("asvs")

# Source the helper script
source('/sw/paul_helper_scripts/pipeline_helper_functions.R')

# Load the ASV table
tbl <- loadQIIMEotutable("asv_table_00.txt")

# Sort columns alphabetically
tbl <- sortQIIMEotutable(tbl, sortby="col", normalize_sort=FALSE)

# Sort rows descending by rowSums
tbl <- sortQIIMEotutable(tbl, sortby="row", normalize_sort=FALSE)

# Save the table to a new file
write.table(tbl, file="asv_table_01.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

