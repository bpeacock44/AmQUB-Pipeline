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
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <corr column>"
    exit 1
fi

# Assign command line arguments to variables
col="$1"

# Warn about missing samples
biom2txt asv_table_02_add_taxa_norm.nc.biom asv_table_02_add_taxa_norm.nc.txt
file1="asv_table_02_add_taxa_norm.nc.txt"
file2="merged_map.txt"

# Extract the second row of file1
second_row=$(awk 'NR==2' "$file1")

# Compare each element of the second row with the first column of file2
IFS=$'\t' read -r -a elements <<< "$second_row"
for element in "${elements[@]}"; do
    if [[ -n $element && $element != "#ASV ID" && $element != "taxonomy" ]]; then
        grep -q "^$element$" <(awk 'NR > 1 {print $1}' "$file2")
        if [[ $? -eq 1 ]]; then
            echo "${element} is present in ${file1} but not in ${file2}. This means it will not be included in the analyis." 
        fi
    fi
done
# Compare each element of the first column of file2 with the elements of the second row of file1
while read -r element; do
    if [[ -n $element ]]; then
        grep -q "$element" <(echo "$second_row" | tr '\t' '\n')
        if [[ $? -eq 1 ]]; then
            echo "${element} is present in ${file2} but not in ${file1}. This means it will not be included in the analyis."
        fi
    fi
done < <(tail -n +2 "$file2" | awk '{print $1}')

echo
echo " - -- --- ---- ---- --- -- -"
echo "Running Correlation Analysis"
echo " - -- --- ---- ---- --- -- -"
Rscript asv_pearson_corr.R asv_table_02_add_taxa_norm.nc.txt merged_map.txt $col

echo
echo " - -- --- ---- ---- --- -- -"
echo "correlation_results.txt has been successfully generated."
echo " - -- --- ---- ---- --- -- -"