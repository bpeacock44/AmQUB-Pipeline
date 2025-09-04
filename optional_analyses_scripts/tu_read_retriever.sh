#!/bin/bash
# See mbio_tutorial.md for further guidance!

### USAGE ###
# This script accepts two types of input:
# 1. Direct command-line arguments:
#    tu_read_retriever.sh -i output_dir -p 0.97 -t otu_list.txt --strategy2

## Required Flags
# -i: The output directory generated in part 3 that you also ran part 4 on.
# -p: Percent identity threshold to use for clustering (this HAS to match your previous run!) Default is 97%.
# -t: A text file that contains a list of the TUs you want to retrieve the reads from.

## Optional Flags:
# --strategy2: Process the output from strategy 2 instead of default strategy.
# --strategy3: Process the output from strategy 3 instead of default strategy. (Requires -o flag.)
# -o: A fasta file of the pre-existing OTUs that were used originally to run Strategy 3. (Requires --strategy3 flag.)

# 2. A parameter template file:
#    tu_read_retriever.sh params.csv
#    Where params.csv contains the following rows, comma delimited, no white space between.
#    The labels at the beginning of each row should be the same as below.
#    The following shows what all possible rows:

#        Part 3 Output Folder To Process,output_dir
#        TU List,otu_list.txt
#        Percent Identity, 0.97
#        Process Strategy 2,true
#        Process Strategy 3,true
#        Strategy 3 Pre-Existing OTUs,pre_existing.fa
#
#Any optional line can be left out of the file if you want to use default settings.

# Set error handling
set -e  # Exit on any error

# Function to parse parameter file
parse_parameter_file() {
    local param_file="$1"
    while IFS= read -r line; do
        case "$line" in
            "Part 3 Output Folder To Process,"*) output_dir="${line#*,}" ;;
            "TU List,"*) tupath="${line#*,}" ;;
            "Percent Identity,"*) TBLID="${line#*,}" ;;
            "Process Strategy 2,"*) STR2_simp="${line#*,}" ;;
            "Process Strategy 3,"*) STR3_simp="${line#*,}" ;;
            "Strategy 3 Pre-Existing OTUs,"*) pre="${line#*,}" ;;
        esac
    done < "$param_file"
}

# Initialize default values
pid=0.97
STR2_simp=false
STR3_simp=false

echo "
        ┌── ===
 ┌──────┤
 │      └── ooo
─┤
 │ ┌── [A m Q U B]  
 └─┤
   └──── ~<>~   

TU Read Retriever

 - -- --- ---- ---- --- -- -"   

# Check if the first argument is a parameter file
if [[ -f "$1" ]]; then
    echo "Reading parameters from file: $1
 - -- --- ---- ---- --- -- -"
    parse_parameter_file "$1"
    shift  # Remove the parameter file argument
else
    # Parse command-line arguments
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -i|--output) output_dir="$2"; shift 2 ;;
            -t|--tu) tupath="$2"; shift 2 ;;
            -p|--per_id) TBLID="$2"; shift 2 ;;
            -o|--str3_pre) pre="$2"; shift 2 ;;
            --nostrategy1) skSTR1=true; shift ;;
            --strategy2) STR2_simp=true; shift ;;
            --strategy3) STR3_simp=true; shift ;;
            *) echo "Unknown option: $1" >&2; exit 1 ;;
        esac
    done
fi

# Check for mandatory arguments
if [ -z "$output_dir" ] || [ -z "$tupath" ] ; then
    echo "Usage: $0 -i <output directory from part 3> -t <path to TU list> [--strategy2 --strategy3 -p <path to pre-existing OTUs>]"
    exit 1
fi

# Make sure output directory and blast database exists
if [[ ! -d "$output_dir" ]]; then
    echo "Error: Output directory '$output_dir' does not exist."
    exit 1
fi

if [[ ! -e "$tupath" ]]; then
    echo "Error: File '$tupath' not found."
    exit 1
fi

if [[ -n "$pre" && ! -e "$pre" ]]; then
    echo "Error: File '$pre' not found."
    exit 1
fi

# Check if the output directory contains the directories "otus" and/or "asvs"
if [[ -d "${output_dir}/otus" && -d "${output_dir}/asvs" ]]; then
    echo "Error: Both 'otus' and 'asvs' directories exist in ${output_dir}. Exiting."
    exit 1
elif [[ -d "${output_dir}/otus" ]]; then
    typ="otu"
elif [[ -d "${output_dir}/asvs" ]]; then
    typ="asv"
else
    echo "Error: Neither 'otus' nor 'asvs' directories exist in ${output_dir}. Exiting."
    exit 1
fi

# check if strategy 2 and 3 have been specified.
if [ "$STR2_simp" = true ] && [ "$STR3_simp" = true ]; then
    echo "Please only process one strategy at a time. To process default strategy, do not use strategy 2 or 3 flags. Exiting."
    exit 1
else
    if [ "$STR2_simp" = true ] || [ "$STR3_simp" = true ]; then
        skSTR1=true
    else
        skSTR1=false
    fi
fi

# detect if STRATEGY 3 and STRATEGY 2 were made and process if they are present.
if [ "$skSTR1" == false ]; then
    DIRS=("${output_dir}/${typ}s")
else 
    DIRS=()
fi

# Check if Strategy 2 should be processed
if [[ "$STR2_simp" == true ]]; then
    if [[ ! -d "${output_dir}/${typ}s/STRATEGY2" ]]; then
        echo "Error: Strategy 2 flag was set, but no STRATEGY2 directory found."
        exit 1
    fi
    DIRS+=("${output_dir}/${typ}s/STRATEGY2/${typ}s")
fi

# Check if Strategy 3 should be processed
if [[ "$STR3_simp" == true ]]; then
    if [[ ! -d "${output_dir}/${typ}s/STRATEGY3" ]]; then
        echo "Error: Strategy 3 flag was set, but no STRATEGY3 directory found."
        exit 1
    fi
    if [[ -z "$pre" ]]; then
        echo "Error: Strategy 3 flag was set, but you didn't specify a pre-existing OTU file (-o flag)."
        exit 1
    fi
    DIRS+=("${output_dir}/${typ}s/STRATEGY3/${typ}s")
fi

# initiate log
timestamp="$(date +"%y%m%d_%H:%M")"
output_file="${output_dir}/tu_read_retriever_${timestamp}.log"
exec > >(tee -ia "$output_file") 2>&1

# log header
echo "Log file for TU Read Retriever. Processed the following arguments:
Output directory to be processed: ${output_dir}
List of TUs to Retrieve: ${tupath}"| tee /dev/tty
if [ "$skSTR1" != false ]; then 
    echo "Strategy 1 (default) results are not being processed."
fi

if [ "$STR2_simp" != false ]; then 
    echo "Strategy 2 results will be processed."
fi

if [ "$STR3_simp" != false ]; then 
    echo "Strategy 3 results will be processed using pre-existing OTUs located at ${pre}."
fi

echo " - -- --- ---- ---- --- -- -"

### Strategy 1 (default)
if [ "$skSTR1" == false ]; then
    # check for and if necessary generate binning reports
    BINNING_REPORT="${output_dir}/${typ}s/read_binning_report.txt"
    if [[ ! -f "$BINNING_REPORT" ]]; then
        echo "Generating read binning report for Strategy 1..."
        usearch --otutab "${output_dir}/combined.fq" -quiet \
            -zotus "${output_dir}/${typ}s/${typ}s.fa" \
            -otutabout "${output_dir}/${typ}s/${typ}_table_00.txt" \
            -id "$TBLID" \
            -mapout "$BINNING_REPORT"
    else
        echo "Strategy 1 binning report already exists: $BINNING_REPORT"
    fi

    # extract reads
    READFA_DIR="${output_dir}/${typ}s/${typ}_read_fastas"
    mkdir -p "$READFA_DIR"  

    echo "Extracting reads for TUs in Strategy 1 using seqkit..."
    while IFS= read -r tu; do
        echo "Extracting reads for ${tu}."
        grep -w "$tu" "$BINNING_REPORT" | cut -f1 > tmp_ids.txt
        num=$(wc -l < tmp_ids.txt)
        echo "There were ${num} reads matching ${tu}."
        if [[ -s tmp_ids.txt ]]; then
            seqkit grep -f tmp_ids.txt "${output_dir}/combined.fq" | seqkit fq2fa > "${READFA_DIR}/${tu}.fasta"
        else
            echo "No reads found for $tu in Strategy 1."
        fi
    done < "$tupath"
    rm -f tmp_ids.txt
fi

### Strategy 2
if [ "$STR2_simp" == true ]; then
    # check for and if necessary generate binning reports
    sub_outdir="${output_dir}/${typ}s/STRATEGY2"
    BINNING_REPORT="${sub_outdir}/${typ}s/read_binning_report.txt"
    if [[ ! -f "$BINNING_REPORT" ]]; then
        echo "Generating read binning report for Strategy 2..."
        usearch --otutab "${output_dir}/combined.fq" -quiet \
            -zotus "${sub_outdir}/${typ}s/${typ}s.fa" \
            -otutabout "${sub_outdir}/${typ}s/${typ}_table_00.txt" \
            -id "$TBLID" \
            -mapout "$BINNING_REPORT"
    else
        echo "Strategy 2 binning report already exists: $BINNING_REPORT"
    fi
    READFA_DIR="${sub_outdir}/${typ}s/${typ}_read_fastas"
    mkdir -p "$READFA_DIR"
    
    echo "Extracting reads for TUs in Strategy 2 using seqkit..."
    while IFS= read -r tu; do
        echo "Extracting reads for ${tu}."
        grep -w "$tu" "$BINNING_REPORT" | cut -f1 > tmp_ids.txt
        num=$(wc -l < tmp_ids.txt)
        echo "There were ${num} reads matching ${tu}."
        if [[ -s tmp_ids.txt ]]; then
            seqkit grep -f tmp_ids.txt "${output_dir}/combined.fq" | seqkit fq2fa > "${READFA_DIR}/${tu}.fasta"
        else
            echo "No reads found for $tu in Strategy 1."
        fi
    done < "$tupath"
    rm -f tmp_ids.txt
fi

### Strategy 3
if [ "$STR3_simp" == true ]; then
    # check for and if necessary generate binning reports
    sub_outdir="${output_dir}/${typ}s/STRATEGY3"
    # Pre-existing binning
    PRE_BINNING="${sub_outdir}/${typ}s/pre-existing_read_binning_report.txt"
    if [[ ! -f "$PRE_BINNING" ]]; then
        echo "Generating pre-existing read binning report for Strategy 3..."
        usearch --otutab "${output_dir}/combined.fq" -quiet \
            -zotus "${pre}" \
            -otutabout "${sub_outdir}/${typ}s/pre-existing_${typ}_table_00.txt" \
            -id "$TBLID" \
            -mapout "$PRE_BINNING"
    else
        echo "Pre-existing binning report already exists: $PRE_BINNING"
    fi

    # Unbinned reads
    UNBINNING_REPORT="${sub_outdir}/${typ}s/unbinned_read_binning_report.txt"
    if [[ ! -f "$UNBINNING_REPORT" ]]; then
        echo "Generating unbinned read binning report for Strategy 3..."
        usearch --otutab "${sub_outdir}/${typ}s/reads_unbinned_in_pre-exisiting.fq" -quiet \
            -zotus "${sub_outdir}/${typ}s/unbinned_${typ}s.fa" \
            -otutabout "${sub_outdir}/${typ}s/unbinned_${typ}_table_00.txt" \
            -id "$TBLID" \
            -mapout "$UNBINNING_REPORT"
    else
        echo "Unbinned read binning report already exists: $UNBINNING_REPORT"
    fi
    READFA_DIR="${sub_outdir}/${typ}s/${typ}_read_fastas"
    mkdir -p "$READFA_DIR"
    
    echo "Extracting reads for TUs in Strategy 3 using seqkit..."
    while IFS= read -r tu; do
        echo "Extracting reads for ${tu}."
        grep -w "$tu" "$PRE_BINNING" | cut -f1 > tmp_ids.txt
        grep -w "$tu" "$UNBINNING_REPORT" | cut -f1 >> tmp_ids.txt
        sort -u tmp_ids.txt > tmp_ids_unique.txt
        num=$(wc -l < tmp_ids.txt)
        echo "There were ${num} reads matching ${tu}."
        if [[ -s tmp_ids.txt ]]; then
            seqkit grep -f tmp_ids.txt "${output_dir}/combined.fq" | seqkit fq2fa > "${READFA_DIR}/${tu}.fasta"
        else
            echo "No reads found for $tu in Strategy 1."
        fi
    done < "$tupath"

    rm -f tmp_ids.txt tmp_ids_unique.txt  
fi
