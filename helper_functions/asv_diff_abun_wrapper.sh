#!/bin/bash

# Usage: ./asv_diff_abun_wrapper.sh <column> <factor 1> <factor 2>
# Choose a column in your mapping file (the merged version from all your data) that contains factors you want to compare.
# e.g. if I have a column "Tissue" and I want to compare "stems" to "leaves" in that column, I will run this script like this:
# ./diff_test.sh Tissue stems leaves
# Note that spelling and capitalization must be the same!

set -e

echo " - -- --- ---- ---- --- -- -"
echo "Loading Qiime1 and helper functions"
echo " - -- --- ---- ---- --- -- -"
source /sw/anaconda3/etc/profile.d/conda.sh
conda activate qiime1
source ~/home/bpeacock_ucr_edu/real_projects/PN94_singularity_of_microbiome_pipeline/targeted_microbiome_via_blast/helper_functions/qiime_shell_helper_functions.sh
module load r

# LD_LIBRARY_PATH modification
LD_LIBRARY_PATH=/home/bpeacock_ucr_edu/.conda/envs/qiime1/lib:$LD_LIBRARY_PATH

echo
echo " - -- --- ---- ---- --- -- -"
echo "Checking arguments and samples"
echo " - -- --- ---- ---- --- -- -"
# Check if all required arguments are provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <asv_table> <mapping_file> <column> <factor 1> <factor 2>"
    exit 1
fi

# Assign command line arguments to variables
asv="$1"
map="$2"
diff_col="$3"
var_a="$4"
var_b="$5"
txtv="${asv%.biom}.txt"

if [ ! -e "$txtv" ]; then
    # Run biom2txt only if $txtv doesn't already exist
    biom2txt "$asv" "$txtv"
fi

# Warn about missing samples

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
echo "Running differential abundance analysis"
echo " - -- --- ---- ---- --- -- -"
# Check if diff_col exists in the header of ${map}
if awk -F'\t' -v dc="$diff_col" 'NR==1{for(i=1;i<=NF;i++) if($i==dc) {dc_flag=1}; if(dc_flag) exit 0; exit 1}' ${map}; then
    # Check if both var_a and var_b are in diff_col
    tissue_col=$(awk -F'\t' 'NR==1 {for(i=1; i<=NF; i++) {if($i=="Tissue") {print i; exit}}}' ${map})
    if awk -F'\t' -v dc="$tissue_col" -v va="$var_a" -v vb="$var_b" 'BEGIN{found_a=0; found_b=0} {if($dc==va) found_a=1; if($dc==vb) found_b=1} END{if(found_a && found_b) exit 0; exit 1}' ${map}; then
        # Run differential_abundance.py
        echo "Running differential_abundance.py"
        differential_abundance.py -i ${asv} -o diff_ASVs.txt -m ${map} -a DESeq2_nbinom -c "$diff_col" -x "$var_a" -y "$var_b" -d

        # Run Rscript asv_diff.R
        echo "Adding average normalized counts to the output."
        Rscript asv_diff.R diff_ASVs.txt ${map} asv_table_02_add_taxa_norm.nc.txt "$diff_col" "$var_a" "$var_b"
    else
        echo "Error: One or both of the specified factors do not exist in $diff_col column of ${map}."
        exit 1
    fi
else
    echo "Error: The specified column $diff_col does not exist in ${map}."
    exit 1
fi

echo
echo " - -- --- ---- ---- --- -- -"
echo "differential_abundance_results.txt and diff_ASVs_diagnostic_plots.pdf have been successfully generated using the differential_abundance.py program from Qiime1."
echo " - -- --- ---- ---- --- -- -"