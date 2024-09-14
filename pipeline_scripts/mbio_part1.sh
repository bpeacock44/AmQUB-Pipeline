#!/bin/bash
# See mbio_tutorial.md for further guidance.

### USAGE ###
# This script expects to be given at least two arguments:
# -d: a working directory, which contains a folder for each of your fastq files named by ID
# -j: all the IDs you intend to process in a comma-delimited list (ID1,ID2,ID3,etc.)

# Optional argument:
# -m: number of mismatched bases (OPTIONAL - if you want to convert barcodes with given number of 
#     mismatches into perfect match barcodes)

set -e

# error handling function
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

# Show your fastq files and map
for JB in "${JBS[@]}"; do
    if [ ! -e "${DIR}/${JB}/${JB}_raw.fq" ]; then
        echo "File ${DIR}/${JB}/${JB}_raw.fq not found!"
        exit 1
    fi
    echo "File ${DIR}/${JB}/${JB}_raw.fq will be processed."

    if [ ! -e "${DIR}/${JB}/${JB}_map.txt" ]; then
        echo "File ${DIR}/${JB}/${JB}_map.txt not found!"
        exit 1
    fi
    echo "File ${DIR}/${JB}/${JB}_map.txt is the corresponding mapping file."
done

# Initialize log for errors
timestamp="$(date +"%Y%m%d_%H:%M:%S")"
output_file="${DIR}/part1_${timestamp}.log"
exec > "$output_file" 2>&1

# Make log file header
echo " - -- --- ---- ---- --- -- -
Log file for Part 1 of the Microbiome Pipeline. Processed the following arguments:
Working Directory: ${DIR}
Data IDs: ${JBS[@]}"
if [ "$mmatchnum" -ne 0 ]; then
    echo "Mismatches specified: ${mmatchnum}"
fi
echo " - -- --- ---- ---- --- -- -"
echo 

if [ "$mmatchnum" -ne 0 ]; then

echo " - -- --- ---- ---- --- -- -
Checking barcode collisions...
 - -- --- ---- ---- --- -- -"

    # Check barcode collisions
    for JB in "${JBS[@]}"; do
        _BC_=$(grep -cP "^[A-Z]" "${DIR}/${JB}/${JB}_map.txt")
        echo "$JB [${_BC_}]"
    
        check_barcode_collisions.pl -i "${DIR}/${JB}/${JB}_raw.fq" -m "${DIR}/${JB}/${JB}_map.txt" -M${mmatchnum} -C -o "${DIR}/${JB}/${JB}.BC${_BC_}_M${mmatchnum}.collisions.txt"
    done
    
    for JB in "${JBS[@]}"; do
        _BC_=$(grep -cP "^[A-Z]" "${DIR}/${JB}/${JB}_map.txt")
    
        echo 
echo " - -- --- ---- ---- --- -- -
Filtering barcode noncollisions
 - -- --- ---- ---- --- -- -"
        
        # Set option for mismatches based on previously selected mismatch value
        case $mmatchnum in
            0) VAR="-m0" ;;
            1) VAR="-m1" ;;
            2) VAR="-m12" ;;
            3) VAR="-m123" ;;
            4) VAR="-m1234" ;;
            5) VAR="-m12345" ;;
        esac
        
        # Filter barcode non-collisions
        filter_barcode_noncollisions.py -k -i "${DIR}/${JB}/${JB}.BC${_BC_}_M${mmatchnum}.collisions.txt" $VAR --output_for_fastq_convert > "${DIR}/${JB}/${JB}_M${mmatchnum}.fbncs" 
        
        echo 
echo " - -- --- ---- ---- --- -- -
Converting mismatches to perfect matches 
 - -- --- ---- ---- --- -- -"
    
        fastq_convert_mm2pm_barcodes.py -t read -i "${DIR}/${JB}/${JB}_raw.fq" -m "${DIR}/${JB}/${JB}_M${mmatchnum}.fbncs" -o "${DIR}/${JB}/${JB}.M${mmatchnum}.fq" 
        extract_barcodes.go -f "${DIR}/${JB}/${JB}.M${mmatchnum}.fq" && mv -v "${DIR}/${JB}/barcodes.fastq" "${DIR}/${JB}/${JB}_BC.M${mmatchnum}.fq" 
    
        # Check if files are present
        [[ -e "${DIR}/${JB}/${JB}.M${mmatchnum}.fq" ]] && echo "File ${DIR}/${JB}/${JB}.M${mmatchnum}.fq was generated." || (echo "Error: File ${DIR}/${JB}/${JB}.M${mmatchnum}.fq was not generated!" && exit 1)
        [[ -e "${DIR}/${JB}/${JB}_BC.M${mmatchnum}.fq" ]] && echo "File ${DIR}/${JB}/${JB}_BC.M${mmatchnum}.fq was generated." || (echo "Error: File ${DIR}/${JB}/${JB}_BC.M${mmatchnum}.fq was not generated!" && exit 1)
    done

else
echo " - -- --- ---- ---- --- -- -
Skipping barcode collision check. 
 - -- --- ---- ---- --- -- -"
    for JB in "${JBS[@]}"; do
        extract_barcodes.go -f "${DIR}/${JB}/${JB}_raw.fq" && mv -v "${DIR}/${JB}/barcodes.fastq" "${DIR}/${JB}/${JB}_BC.M${mmatchnum}.fq" 
        ln -sf "${DIR}/${JB}/${JB}_raw.fq" "${DIR}/${JB}/${JB}.M${mmatchnum}.fq"       
    done
fi

# Print number of reads per barcode
echo | tee /dev/tty

for JB in "${JBS[@]}"; do
    echo " - -- --- ---- ---- --- -- -" | tee /dev/tty
    bc_counter.py "${DIR}/${JB}/${JB}_map.txt" "${DIR}/${JB}/${JB}_BC.M${mmatchnum}.fq" "${DIR}/${JB}/${JB}_read_counts_M${mmatchnum}.txt"
    
echo "The number of reads per sample that resulted from this script for ${JB} can be found in this file: 
${DIR}/${JB}/${JB}_read_counts_M${mmatchnum}.txt
    
You should check this file, as it may indicate that you should remove or ignore certain samples downstream.
Here are the first few lines:" | tee /dev/tty
    
    head "${DIR}/${JB}/${JB}_read_counts_M${mmatchnum}.txt" | tee /dev/tty
    echo " - -- --- ---- ---- --- -- -" | tee /dev/tty
    echo | tee /dev/tty
done

# Final message - what is next
echo | tee /dev/tty
echo " - -- --- ---- ---- --- -- -
Part 1 completed. Part 2 will run on each of your data files individually. 
You may also specify a new mapping file for -o if you want to run further analyses 
on a subset of your demultiplexed data.

For example, next you might run: 
part2.sh -d ${DIR} -j \"${JBS[0]}\" (-o if you want a subset)
 - -- --- ---- ---- --- -- -" | tee /dev/tty
