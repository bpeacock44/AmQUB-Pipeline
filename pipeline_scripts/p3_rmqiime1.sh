#!/bin/bash
# See mbio_tutorial.md for further guidance!

### USAGE ###
#This script expects to be given at least 4 arguments:
#-d: a working directory, which contains one folder for each of your fastq files named by ID
#-j: the folders created in the last part that you intend to process in a comma-delimited list 
    #(ID1_subset1_output, ID2_output, ID3_subset2_output, etc.)
#-l: the length you want to trim your reads to. Note ALL files will be trimmed to this length.
#-o: the name of your output directory

#Optional arguments:
#-m: number of mismatches, if using (again, this should have been specified from part1)
#-n: change the minsize of the unoise3 algorithm (default is 8)

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
while getopts ":d:j:l:o:m:n" opt; do
  case $opt in
    d) DIR="$OPTARG"
    ;;
    j) IFS=',' read -ra JBS <<< "$OPTARG"
    ;;
    l) LEN="$OPTARG"
    ;;
    o) OUTDIR="$OPTARG"
    ;;    
    m) mmatchnum="$OPTARG"
    ;;
    n) MIN="$OPTARG"
    ;;     
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

# Check for mandatory arguments
if [ -z "$DIR" ] || [ -z "${JBS[*]}" ] || [ -z "$LEN" ] || [ -z "$OUTDIR" ]; then
    echo "Usage: $0 -d <directory_path> -j <comma-separated output directories> -l <trim length> -o <desired name of output dir> [-m <mismatch number> -n <min size for asv calling>]"
    exit 1
fi

# Check if mmatchnum is not defined or is set to an empty string, set it to "0"
if [ -z "${mmatchnum+x}" ] || [ -z "$mmatchnum" ]; then
    mmatchnum="0"
fi

# Remove "_output" from each element in the array
for ((i=0; i<${#JBS[@]}; i++)); do
  JBS[$i]=${JBS[$i]%%_output}
done

# Check appropriateness fo the trim length.
for JB in "${JBS[@]}"; do
    file="${DIR}/${JB}_output/${JB}_A1P1.M${mmatchnum}.fq"

    # Check if the file exists and is not empty
    if [ -s "$file" ]; then
        # Get the length of the first sequence in the FASTA file
        first_sequence_length=$(head -n 2 "$file" | grep -v '^>' | tr -d '\n' | wc -c)

        # Compare the length with $LEN
        if [ "$first_sequence_length" -lt "$(($LEN - 5))" ]; then
            echo "Your chosen trim length (-l) is inappropriate for the length of your reads."
            exit 1
        fi
    else
        echo "Error while attempting to check read length: File $file doesn't exist or is empty. This file needs to be present!"
        exit 1
    fi

    break
done

echo " - -- --- ---- ---- --- -- -"
echo "Checking for input files"
echo " - -- --- ---- ---- --- -- -"

# show your fastq files and map
for JB in "${JBS[@]}"; do
    if [ ! -e "${DIR}/${JB}_output/${JB}_A1P1.M${mmatchnum}.fq" ]; then
        echo "File ${DIR}/${JB}_output/${JB}_A1P1.M${mmatchnum}.fq not found!"
        exit 1
    fi

    if [ ! -e "${DIR}/${JB}_output/${JB}_A1P2.M${mmatchnum}.fq" ]; then
        echo "File ${DIR}/${JB}_output/${JB}_A1P2.M${mmatchnum}.fq not found!"
        exit 1
    fi
done

# initiate log
timestamp="$(date +"%Y%m%d_%H:%M:%S")"
output_file="${DIR}/part3_${timestamp}.log"
exec > "$output_file" 2>&1

# log header
echo " - -- --- ---- ---- --- -- -"
echo "Log file for Part 3 of the Microbiome Pipeline. Processed the following arguments:
Working Directory: ${DIR}
Data IDs: ${JBS[@]}
Trim Length: ${LEN}
Output Directory: ${OUTDIR}"

if [ "$mmatchnum" -ne 0 ]; then
    echo "Mismatches specified: ${mmatchnum}"
fi

if [ -n "$MIN" ]; then
    echo "UNOISE3 was run with minsize ${mmatchnum}"
fi

echo " - -- --- ---- ---- --- -- -"

echo
echo " - -- --- ---- ---- --- -- -"
echo "Trimming and Filtering Reads"
echo " - -- --- ---- ---- --- -- -"

# set up and go to output directory
output_dir="${DIR}/${OUTDIR}"
mkdir -vp "${DIR}/${OUTDIR}"
echo "This folder contains the results of ${JBS[@]} trimmed at ${LEN}." > "${output_dir}/summary.txt"

# make tax directory
mkdir -vp "${DIR}/${OUTDIR}/tax_dir"
TAXDIR="${DIR}/${OUTDIR}/tax_dir"

#truncate reads at LEN
for JB in ${JBS[@]}; do
  usearch -fastx_truncate "${DIR}/${JB}_output/${JB}_A1P1.M${mmatchnum}.fq" -quiet -trunclen ${LEN} -fastqout "${DIR}/${JB}_output/${JB}_A1P1_${LEN}bp.fq"; # &
done

rm -f "${output_dir}/combined.fq"
if [ ${#JBS[@]} -gt 1 ]; then
  #EITHER combine if you have multiple files
  for JB in ${JBS[@]}; do
  cat "${DIR}/${JB}_output/${JB}_A1P1_${LEN}bp.fq" >> "${output_dir}/combined.fq"
  done
else
  #OR create a symlink if you have only one file
  ln -s "${DIR}/${JB}_output/${JB}_A1P1_${LEN}bp.fq" "${output_dir}/combined.fq"
fi

#maxee quality filtering of demultiplexed/truncated fq files (*** keep THREADS=1 for repeatability ***)
for JB in ${JBS[@]}; do
  usearch -threads 1 -fastq_filter "${DIR}/${JB}_output/${JB}_A1P1_${LEN}bp.fq" -quiet -fastq_maxee 1.0 -fastaout "${DIR}/${JB}_output/${JB}.filtered.fa"
done
echo
echo " - -- --- ---- ---- --- -- -"
echo "Pooling Samples and Creating ASVs"
echo " - -- --- ---- ---- --- -- -"

#sample pooling (https://www.drive5.com/usearch/manual/pool_samples.html)
rm -f "${output_dir}/filtered.fa"
if [ ${#JBS[@]} -gt 1 ]; then
  #EITHER combine if you have multiple files
  for JB in ${JBS[@]}; do
  cat "${DIR}/${JB}_output/${JB}.filtered.fa" >> "${output_dir}/filtered.fa"
  done
else
  #OR create a symlink if you have only one file
  ln -s "${DIR}/${JB}_output/${JB}.filtered.fa" "${output_dir}/filtered.fa"
fi

#find unique sequences
usearch -fastx_uniques "${output_dir}/filtered.fa" -quiet -fastaout "${output_dir}/uniques.fa" -sizeout -relabel Uniq

#make a subdirectory for the asvs
mkdir -vp "${output_dir}/asvs"

if [ -n "$MIN" ]; then
    # Cluster unique sequences into ASVs using the UNOISE3 algorithm with custom minsize
    usearch -unoise3 "${output_dir}/uniques.fa" -quiet -minsize "${MIN}" -zotus "${output_dir}/asvs/asvs.fa"
else
    # Cluster unique sequences into ASVs using the UNOISE3 algorithm with default minsize (8)
    usearch -unoise3 "${output_dir}/uniques.fa" -quiet -zotus "${output_dir}/asvs/asvs.fa"
fi

# Convert '>Zotu' to '>Asv' in the file
sed 's/>Zotu/>Asv/g' "${output_dir}/asvs/asvs.fa" > "${output_dir}/asvs/z.fa"

# Check if the replacement was successful before overwriting
if grep -q '>Asv' "${output_dir}/asvs/z.fa"; then
    # Overwrite the original file if '>Asv' is found
    mv -v "${output_dir}/asvs/z.fa" "${output_dir}/asvs/asvs.fa"
else
    echo "Header replacement failed. Original file not overwritten."
fi
echo
echo " - -- --- ---- ---- --- -- -"
echo "Creating Initial ASV Table"
echo " - -- --- ---- ---- --- -- -"

#create an ASV table ("Input should be reads before quality filtering and before discarding low-abundance unique sequences, e.g. singletons")
usearch --otutab "${output_dir}/combined.fq" -quiet -zotus "${output_dir}/asvs/asvs.fa" -otutabout "${output_dir}/asvs/asv_table_00.txt"
sed -i 's/#OTU/#ASV/g' "${output_dir}/asvs/asv_table_00.txt"

# Run python script to sort qiime table
qiime_table_sorter.py "${output_dir}/asvs/asv_table_00.txt" "${output_dir}/asvs/asv_table_01.txt"

#add counts to ASV file
add_counts_to_fasta.py "${output_dir}/asvs/asv_table_01.txt" "${output_dir}/asvs/asvs.fa" "${output_dir}/asvs/asvs_counts.fa"

echo
echo " - -- --- ---- ---- --- -- -"
echo "It's time to assign taxonomy! You will either do this by BLAST-ing your sequences against a 
local BLAST database or by using Qiime2 to classify them." | tee /dev/tty








