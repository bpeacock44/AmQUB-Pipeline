#!/bin/bash

# Usage: ./otu_prism_summary_wrapper.sh 
# Choose a column in your mapping file (the merged version from all your data) that contains numbers you want to correlate your OTUs by.
# e.g. if I have a column "Rating" with disease ratings in it, I will run this script like this:
# ./diff_test.sh Rating
# Note that spelling and capitalization must be the same!

echo " - -- --- ---- ---- --- -- -"
echo "Loading Qiime1 and helper functions"
echo " - -- --- ---- ---- --- -- -"
source /sw/anaconda3/etc/profile.d/conda.sh
conda activate qiime1
source ~/paul_helper_scripts/qiime_shell_helper_functions.sh
module load r

echo
echo " - -- --- ---- ---- --- -- -"
echo "Checking arguments and samples"
echo " - -- --- ---- ---- --- -- -"
# Check if all required arguments are provided
if [ "$#" -ne 1 ]; then
            echo "Usage: $0 <corr column>"
                exit 1
fi

# Assign command line arguments to variables
col="$1"

# format James wants:
#       Microbe 1   Microbe 2   Microbe 3
#Treatment Group 1  42.03   17.61   1.52
#Treatment Group 2  14.66   17.90   4.09
#Treatment Group 3  17.28   16.85   2.

#Lump treatment groups together (based on same column as edgeR) - maybe user lists treatments they want included.
#Transpose
#Determine top 20 taxa across all (number may change)
#Combine other taxa to "other" category
#Change taxa names so they just contain the rank