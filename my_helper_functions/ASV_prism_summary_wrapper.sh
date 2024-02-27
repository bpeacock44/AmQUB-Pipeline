#!/bin/bash

# Choose a column in your mapping file (the merged version from all your data) that contains the treatments you want to group samples by.
# Input will be the column followed by a comma-delimited list of treatments you want to include.
# Also put the number of taxa you want in the final output. (i.e. The number that would be included in a bar plot.)
# e.g. if I have a column "tissue" with various types of tissues in it, I might run this script like this:
# ./ASV_prism_summary_wrapper.sh tissue stem_whole,leaf_scrapings,leaf_whole 20
# Note that spelling and capitalization must be the same!

module load r

# Check if all required arguments are provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <ASV table path> <mapping file> <column name> <comma-delimited treatment list> <number of taxa desired in final output>"
    exit 1
fi

# Assign command line arguments to variables
file="$1"
map="$2"
col="$3"
treatments="$4"
num="$5"

Rscript ASV_prism_summary.R $file $map $col $treatments $num

