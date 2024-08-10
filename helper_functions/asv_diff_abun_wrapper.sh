#!/bin/bash

set -e

echo
echo " - -- --- ---- ---- --- -- -"
echo "Checking arguments and samples"
echo " - -- --- ---- ---- --- -- -"

# Assign command line arguments to variables
edgeR_option=
DESeq2_option=

while getopts ":t:m:c:v1:v2:RD" opt; do
    case $opt in
        t) asv="$OPTARG";;
        m) map="$OPTARG";;
        c) diff_col="$OPTARG";;
        1) var_1="$OPTARG";;
        2) var_2="$OPTARG";;
        R) edgeR_option=1;;
        D) DESeq2_option=1;;
        \?) echo "Invalid option: -$OPTARG" >&2
            exit 1;;
        :) echo "Option -$OPTARG requires an argument." >&2
            exit 1;;
    esac
done

# Check if exactly one of -R or -D is provided
#if [[ -z "$edgeR_option" && -z "$DESeq2_option" ]]; then
    #echo "Error: Either -R (edgeR) or -D (DESeq2) must be specified." >&2
    #exit 1
#elif [[ -n "$edgeR_option" && -n "$DESeq2_option" ]]; then
    #echo "Error: Only one of -R (edgeR) or -D (DESeq2) can be specified, not both." >&2
    #exit 1
#fi

edgeR_option=1

# Check if all required arguments are provided
if [[ -z "$asv" || -z "$map" || -z "$diff_col" || -z "$var_1" || -z "$var_2" ]]; then
    echo "Error: Missing required arguments. Please provide -t, -m, -c, -1, -2." >&2
    exit 1
fi

txtv="${asv%.biom}.txt"

# Check if diff_col exists in the header of ${map}
if awk -F'\t' -v dc="$diff_col" 'NR==1{for(i=1;i<=NF;i++) if($i==dc) {dc_flag=1}; if(dc_flag) exit 0; exit 1}' "${map}"; then
    # Check if both var_1 and var_2 are in diff_col
    tissue_col=$(awk -F'\t' 'NR==1 {for(i=1; i<=NF; i++) {if($i=="'"${diff_col}"'") {print i; exit}}}' "${map}")
    if awk -F'\t' -v dc="$tissue_col" -v va="$var_1" -v vb="$var_2" 'BEGIN{found_a=0; found_b=0} {if($dc==va) found_a=1; if($dc==vb) found_b=1} END{if(found_a && found_b) exit 0; exit 1}' "${map}"; then
        echo "Column and variables check passed."
    else
        echo "Error: One or both of the specified factors do not exist in $diff_col column of ${map}."
        exit 1
    fi
else
    echo "Error: The specified column $diff_col does not exist in ${map}."
    exit 1
fi

if [[ -n "$edgeR_option" ]]; then
    echo
    echo " - -- --- ---- ---- --- -- -"
    echo "Running edgeR differential abundance analysis"
    echo " - -- --- ---- ---- --- -- -"

    # Assuming edgeR_diff_abundance.py is your Python script for edgeR analysis
    edgeR_diff_abundance.R -t "${txtv}" -o edgeR_${txtv}_${var_1}.v.${var_2}.txt -m "${map}" -d "${diff_col}" -v1 "${var_1}" -v2 "${var_2}"
fi

if [[ -n "$DESeq2_option" ]]; then
    echo "DESeq2 not currently available for this container."
    #echo
    #echo " - -- --- ---- ---- --- -- -"
    #echo "Running DESeq2 differential abundance analysis"
    #echo " - -- --- ---- ---- --- -- -"
#
    ## Assuming differential_abundance.py is your script for DESeq2 analysis
    #differential_abundance.py -i "${asv}" -o deseq2_${asv}_${var_1}.v.${var_2}.txt -m "${map}" -a DESeq2_nbinom -c "${diff_col}" -x "${var_1}" -y "${var_2}" -d
fi
