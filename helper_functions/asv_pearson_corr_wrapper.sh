#!/bin/bash

# Usage: ./asv_pearson_corr_wrapper.sh <column>
# Choose a column in your mapping file (the merged version from all your data) that contains numbers you want to correlate your ASVs by.
# e.g. if I have a column "Rating" with disease ratings in it, I will run this script like this:
# ./diff_test.sh Rating
# Note that spelling and capitalization must be the same!

set -e

echo " - -- --- ---- ---- --- -- -"
echo "Loading Qiime1 and helper functions"
echo " - -- --- ---- ---- --- -- -"
source /sw/anaconda3/etc/profile.d/conda.sh
conda activate qiime1
source ~/home/bpeacock_ucr_edu/real_projects/PN94_singularity_of_microbiome_pipeline/targeted_microbiome_via_blast/helper_functions/qiime_shell_helper_functions.sh
module load r

echo
echo " - -- --- ---- ---- --- -- -"
echo "Checking arguments and samples"
echo " - -- --- ---- ---- --- -- -"
# Check if all required arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <asv_table> <mapping_file> <corr column>"
    exit 1
fi

# Assign command line arguments to variables
asv="$1"
map="$2"
col="$3"
txtv="${asv%.biom}.txt"

# Warn about missing samples
if [ ! -e "$txtv" ]; then
    # Run biom2txt only if $txtv doesn't already exist
    biom2txt "$asv" "$txtv"
fi

# Extract the second row of txtv
second_row=$(awk 'NR==2' "$txtv")

# Compare each element of the second row with the first column of map
IFS=$'\t' read -r -a elements <<< "$second_row"
for element in "${elements[@]}"; do
    if [[ -n $element && $element != "#ASV ID" && $element != "taxonomy" ]]; then
        grep -q "^$element$" <(awk 'NR > 1 {print $1}' "$map")
        if [[ $? -eq 1 ]]; then
            echo "${element} is present in ${txtv} but not in ${map}. This means it will not be included in the analyis." 
        fi
    fi
done
# Compare each element of the first column of map with the elements of the second row of txtv
while read -r element; do
    if [[ -n $element ]]; then
        grep -q "$element" <(echo "$second_row" | tr '\t' '\n')
        if [[ $? -eq 1 ]]; then
            echo "${element} is present in ${map} but not in ${txtv}. This means it will not be included in the analyis."
        fi
    fi
done < <(tail -n +2 "$map" | awk '{print $1}')

echo
echo " - -- --- ---- ---- --- -- -"
echo "Running correlation analysis"
echo " - -- --- ---- ---- --- -- -"
Rscript asv_pearson_corr.R $txtv $map $col

echo
echo " - -- --- ---- ---- --- -- -"
echo "correlation_results.txt has been successfully generated."
echo " - -- --- ---- ---- --- -- -"