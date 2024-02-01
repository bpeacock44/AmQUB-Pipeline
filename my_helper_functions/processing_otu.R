#!/usr/bin/env Rscript
setwd("zotus")

# Source the helper script
source('/sw/paul_helper_scripts/pipeline_helper_functions.R')

# Load the OTU table
tbl <- loadQIIMEotutable("otu_table_00.txt")

# Sort columns alphabetically
tbl <- sortQIIMEotutable(tbl, sortby="col", normalize_sort=FALSE)

# Sort rows descending by rowSums
tbl <- sortQIIMEotutable(tbl, sortby="row", normalize_sort=FALSE)

# Save the table to a new file
write.table(tbl, file="otu_table_01.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

