#!/bin/bash

### USAGE ###
# This script expects to be given at least two aguments:
# -d: a working directory, which contains a one folder for each of your fastq files named by ID
# -j: all the IDs you intend to process in a comma-delimited list (ID1,ID2,ID3,etc.)

# Optional argument:
# -m: number of mismatched bases (OPTIONAL - if you want to convert barcode with given number of mismatches into perfect match barcodes)

# Examples:
# mbio_part1.sh -d /path/to/dir -j JB141,JB143 
# mbio_part1.sh -d /path/to/dir -j JB141,JB143 -m 1

### INPUT ###
# Each folder needs to contain a fastq file named by ID followed by "_L1P1.fq" and an appropriately 
#       formatted mapping file followed by "_map.txt"
# For example, the directory indicated contains a folder called JB141 and it contains the files 
#       "JB141_L1P1.fq" and "JB141_map.txt."



# ########## MAPPING FILE FORMAT ########## 
# MAP FILES MUST CONTAIN AT LEAST TWO COLUMNS:  
# 1) SampleID (with a # before as in #SampleID) - these are the IDs you associate with each sample
# 2) BarcodeSequence - the barcodes for each sample
# 3) SampleType - you can use this to filter later on - 
# e.g. removing controls before data analysis, removing samples that aren't relevant, etc.
# 4) Any characteristics you want to us later on to run differential analysis or correlation, etc.

# Other metadata can also be included in additional columns as desired. Here is a short example:

#SampleID    BarcodeSequence    SampleType    PlatePosition   Library    TubeLabel    Contents    DateTaken
#B001.110    CTCGACTACTGA    SAMPLE    A1    JB110    1    Psyllid 1-6    2/28/19
#B002.110    TGACCAGTAGTC    SAMPLE    A2    JB110    2    Psyllid 7-12    2/28/19
#B003.110    GCGATTAGGTCG    IGNORE    A3    JB110    3    Psyllid 13-18    2/28/19
#PCR_CONTROL   ACATGGCCTAAT    CONTROL    A4    JB110    NA    NA    NA
# ########## MAPPING FILE FORMAT ########## 



# <> # TO DO:
# <> # conda qiime1 activation??

# CODE FOLLOWS HERE #

set -e

# Custom error handler function
error_handler() {
    local error_message=$1
    echo "Error on line $error_message" | tee /dev/tty
}

# Trap errors and call the error handler
trap 'error_handler "$BASH_COMMAND"' ERR

# ARGUMENTS
while getopts ":d:j:m:" opt; do
  case $opt in
    d) DIR="$OPTARG"
    ;;
    j) IFS=',' read -ra JBS <<< "$OPTARG"
    ;;
    m) mmatchnum="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

# Check for mandatory arguments
if [ -z "$DIR" ] || [ -z "$JBS" ]; then
    echo "Usage: $0 -d <directory_path> -j <data_ID> [-m <optional mismatch number>]"
    exit 1
fi

# Check if mmatchnum is not defined, set it to "0"
if [ -z "$mmatchnum" ]; then
    mmatchnum="0"
fi

echo " - -- --- ---- ---- --- -- -
Checking for input files
 - -- --- ---- ---- --- -- -"

# show your fastq files and map
for JB in "${JBS[@]}"; do
    if [ ! -e "${DIR}/${JB}/${JB}_L1P1.fq" ]; then
        echo "File ${DIR}/${JB}/${JB}_L1P1.fq not found!"
        exit 1
    fi
    echo "File ${DIR}/${JB}/${JB}_L1P1.fq will be processed."

    if [ ! -e "${DIR}/${JB}/${JB}_map.txt" ]; then
        echo "File ${DIR}/${JB}/${JB}_map.txt not found!"
        exit 1
    fi
    echo "File ${DIR}/${JB}/${JB}_map.txt is the corresponding mapping file."
done

# initialize log for errors
timestamp="$(date +"%Y%m%d_%H:%M:%S")"
output_file="${DIR}/part1_${timestamp}.log"
exec > "$output_file" 2>&1

# make log file header
echo " - -- --- ---- ---- --- -- -
Log file for Part 1 of the Microbiome Pipeline. Processing the following arguments:
Working Directory: ${DIR}
Data IDs: ${JBS[@]}
Mismatches if specified: ${mmatchnum}
 - -- --- ---- ---- --- -- -"
echo 
echo " - -- --- ---- ---- --- -- -
Checking barcode collisions...
 - -- --- ---- ---- --- -- -"

#1) check_barcode_collisions.pl
for JB in "${JBS[@]}"; do
    _BC_=$(grep -cP "^[A-Z]" "${DIR}/${JB}/${JB}_map.txt")
    echo "$JB [${_BC_}]"

    # Check for barcode collisions
    check_barcode_collisions.pl -i "${DIR}/${JB}/${JB}_L1P1.fq" -m "${DIR}/${JB}/${JB}_map.txt" -M${mmatchnum} -C -o "${DIR}/${JB}/uFQBC_${JB}_L1P1.fq_BC${_BC_}_M${mmatchnum}.txt" 
done

# must be in qiime1 env for the next one.
source /sw/miniconda3/bin/activate qiime1

for JB in ${JBS[@]}; do
    # Define _BC_ again.
    _BC_=$(grep -cP "^[A-Z]" "${DIR}/${JB}/${JB}_map.txt")

echo 
echo " - -- --- ---- ---- --- -- -
Filtering barcode noncollisions
 - -- --- ---- ---- --- -- -"

# set option for mismatches based on previously selected mismath value
case $mmatchnum in
    0)
        VAR="-m0"
        ;;
    1)
        VAR="-m1"
        ;;
    2)
        VAR="-m12"
        ;;
    3)
        VAR="-m123"
        ;;
    4)
        VAR="-m1234"
        ;;
    5)
        VAR="-m12345"
        ;;
esac

    # Filter barcode non-collisions
    filter_barcode_noncollisions.py -k -i "${DIR}/${JB}/uFQBC_${JB}_L1P1.fq_BC${_BC_}_M${mmatchnum}.txt" $VAR --output_for_fastq_convert > "${DIR}/${JB}/${JB}_M${mmatchnum}.fbncs" 

echo 
echo " - -- --- ---- ---- --- -- -
Converting mismatches to perfect matches 
 - -- --- ---- ---- --- -- -"

    # Convert mismatches to perfect matches and extract barcodes
    fastq_convert_mm2pm_barcodes.py -t read -i "${DIR}/${JB}/${JB}_L1P1.fq" -m "${DIR}/${JB}/${JB}_M${mmatchnum}.fbncs" -o "${DIR}/${JB}/${JB}_A1P1.M${mmatchnum}.fq" 
    extract_barcodes.go -f "${DIR}/${JB}/${JB}_A1P1.M${mmatchnum}.fq" && mv -v "${DIR}/${JB}/barcodes.fastq" "${DIR}/${JB}/${JB}_A1P2.M${mmatchnum}.fq" 

    # Check if files are present
    [[ -e "${DIR}/${JB}/${JB}_A1P1.M${mmatchnum}.fq" ]] && echo "File ${DIR}/${JB}/${JB}_A1P1.M${mmatchnum}.fq was generated." || (echo "Error: File ${DIR}/${JB}/${JB}_A1P1.M${mmatchnum}.fq was not generated!" && exit 1)
    [[ -e "${DIR}/${JB}/${JB}_A1P2.M${mmatchnum}.fq" ]] && echo "File ${DIR}/${JB}/${JB}_A1P2.M${mmatchnum}.fq was generated." || (echo "Error: File ${DIR}/${JB}/${JB}_A1P2.M${mmatchnum}.fq was not generated!" && exit 1)
done

# deactivate qiime1
conda deactivate

# print number of reads per barcode
echo | tee /dev/tty

for JB in ${JBS[@]}; do
echo " - -- --- ---- ---- --- -- -" | tee /dev/tty
    _BC_=$(grep -cP "^[A-Z]" "${DIR}/${JB}/${JB}_map.txt")
    grep tot "${DIR}/${JB}/uFQBC_${JB}_L1P1.fq_BC${_BC_}_M${mmatchnum}.txt" | awk -v OFS='\t' '{print $6, $2, $3}' | sed '1i sample_ID\tbarcode\tread_count' > "${DIR}/${JB}/${JB}_read_counts_M${mmatchnum}.txt" 
echo "The number of reads per sample that resulted from this script for ${JB} can be found in this file: 
"${DIR}/${JB}/${JB}_read_counts_M${mmatchnum}.txt" 

You should check this file, as it may indicate that you should remove or ignore certain samples downstream.
Here are the first few lines:"  | tee /dev/tty
head "${DIR}/${JB}/${JB}_read_counts_M${mmatchnum}.txt" | tee /dev/tty
echo " - -- --- ---- ---- --- -- -" | tee /dev/tty
echo | tee /dev/tty
done

# final message - what is next
echo | tee /dev/tty
echo " - -- --- ---- ---- --- -- -
Part 1 completed. Part 2 will run on each of your data files individually. 
You may also specify a new mapping file for -o if you want to run further analyses 
on a subset of your demultiplexed data.

For example, next you might run: 
part2.sh -d ${DIR} -j "${JBS[0]}" (-o if you want a subset)
 - -- --- ---- ---- --- -- -" | tee /dev/tty