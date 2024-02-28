#!/bin/bash

# Usage: ./asv_diff_abun_wrapper.sh <column> <factor 1> <factor 2>
# Choose a column in your mapping file (the merged version from all your data) that contains factors you want to compare.
# e.g. if I have a column "Tissue" and I want to compare "stems" to "leaves" in that column, I will run this script like this:
# ./diff_test.sh Tissue stems leaves
# Note that spelling and capitalization must be the same!

echo " - -- --- ---- ---- --- -- -"
echo "Loading Qiime1 and helper functions"
echo " - -- --- ---- ---- --- -- -"
source /sw/anaconda3/etc/profile.d/conda.sh
conda activate qiime1
source ~/paul_helper_scripts/qiime_shell_helper_functions.sh
module load r

# LD_LIBRARY_PATH modification
LD_LIBRARY_PATH=/home/bpeacock_ucr_edu/.conda/envs/qiime1/lib:$LD_LIBRARY_PATH

echo
echo " - -- --- ---- ---- --- -- -"
echo "Checking arguments and samples"
echo " - -- --- ---- ---- --- -- -"
# Check if all required arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <column> <factor 1> <factor 2>"
    exit 1
fi

# Assign command line arguments to variables
diff_col="$1"
var_a="$2"
var_b="$3"

# Warn about missing samples
biom2txt asv_table_02_add_taxa.nc.biom asv_table_02_add_taxa.nc.txt
file1="asv_table_02_add_taxa.nc.txt"
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
echo "Running differential abundance analysis"
echo " - -- --- ---- ---- --- -- -"
# Check if diff_col exists in the header of merged_map.txt
if awk -F'\t' -v dc="$diff_col" 'NR==1{for(i=1;i<=NF;i++) if($i==dc) {dc_flag=1}; if(dc_flag) exit 0; exit 1}' merged_map.txt; then
    # Check if both var_a and var_b are in diff_col
    tissue_col=$(awk -F'\t' 'NR==1 {for(i=1; i<=NF; i++) {if($i=="Tissue") {print i; exit}}}' merged_map.txt)
    if awk -F'\t' -v dc="$tissue_col" -v va="$var_a" -v vb="$var_b" 'BEGIN{found_a=0; found_b=0} {if($dc==va) found_a=1; if($dc==vb) found_b=1} END{if(found_a && found_b) exit 0; exit 1}' merged_map.txt; then
        # Run differential_abundance.py
        echo "Running differential_abundance.py"
        differential_abundance.py -i asv_table_02_add_taxa.nc.biom -o diff_ASVs.txt -m merged_map.txt -a DESeq2_nbinom -c "$diff_col" -x "$var_a" -y "$var_b" -d

        # Run Rscript asv_diff.R
        echo "Adding average normalized counts to the output."
        Rscript asv_diff.R diff_ASVs.txt merged_map.txt asv_table_02_add_taxa_norm.nc.txt "$diff_col" "$var_a" "$var_b"
    else
        echo "Error: One or both of the specified factors do not exist in $diff_col column of merged_map.txt."
        exit 1
    fi
else
    echo "Error: The specified column $diff_col does not exist in merged_map.txt."
    exit 1
fi

echo
echo " - -- --- ---- ---- --- -- -"
echo "differential_abundance_results.txt and diff_ASVs_diagnostic_plots.pdf have been successfully generated using the differential_abundance.py program from Qiime1."
echo " - -- --- ---- ---- --- -- -"