#!/bin/bash

### USAGE ###
# This script expects to be given at least two aguments:
# -d: a working directory, which contains the folder containing the fastq file you want to process.
# -j: a single ID. This script must be run individually on your IDs 
# (This is in contrast to part 1, which was run just once for all.)

# Optional arguments:
# -m: the number of mismatches you want to use. This needs to match the files you generated in part 1.
# -o: a subset ID - if you want to run further analyses on a subset of the samples in your data, 
# you can create a mapping file in the same format as the original with the lines of unwanted 
# samples removed. This file will be named ID_map.subsetID.txt (e.g. JB141_map.Nickels01.txt) and be 
# placed in the same ID folder as the other files are in.

# Examples:
# ./mbio_part2.sh -d /path/to/dir -j JB141 
# ./mbio_part2.sh -d /path/to/dir -j JB141 -o Nickels01 
# ./mbio_part2.sh -d /path/to/dir -j JB141 -m 1
# ./mbio_part2.sh -d /path/to/dir -j JB141 -o Nickels01 -m 1

### INPUT ###
# This script can only be run once the original fastq file (e.g. JB141_L1P1.fq) has been run through part 1, 
#       which is used to find and replace mismatched barcodes with perfect match barcodes. 

# Each folder needs to contain the fastq files resulting from 1a, which are named by ID followed by 
#       _A1P1.M#.fq and _A1P2.M#.fq, as well as a mapping file (either the original or a subset.)

# So, as an example, your working directory might now include:
#       Folder JB141 (containing JB141_A1P1.M0.fq, JB141_A1P2.M0.fq, and JB141_map.Nickels01.txt)
#       (JB141_map.txt should also be present in folder if subset map isn't used.)

# When this code is run, a new directory will be created for your output named either with the unique identifier
#       for your subset, if given (e.g. JB141_Nickels01_output), or it will be named after your regular ID if no unique map 
#       was provided (e.g. JB141_output) 

# <> # TO DO:
# <> # HDIR FILES SHOULD BE IN PATH

# CODE FOLLOWS HERE #

set -e

# ARGUMENTS
while getopts ":d:j:o:m:" opt; do
  case $opt in
    d) DIR="$OPTARG"
    ;;
    j) JB="$OPTARG"
    ;;
    o) OPTION_O="$OPTARG"
    ;;
    m) mmatchnum="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

# Check for mandatory arguments
if [ -z "$DIR" ] || [ -z "$JB" ]; then
    echo "Usage: $0 -d <directory_path> -j <data_ID> [-o <optional alternative_map> -m <optional mismatch number>]"
    exit 1
fi

# Check if mmatchnum is not defined or is set to an empty string, set it to "0"
if [ -z "${mmatchnum+x}" ] || [ -z "$mmatchnum" ]; then
    mmatchnum="0"
fi

echo " - -- --- ---- ---- --- -- -"
echo "Checking for input files"
echo " - -- --- ---- ---- --- -- -"

HDIR=/home/bpeacock_ucr_edu/real_projects/PN94_singularity_of_microbiome_pipeline/targeted_microbiome_via_blast/paul_helper_functions

# show your fastq files 
if [ ! -e "${DIR}/${JB}/${JB}_A1P1.M${mmatchnum}.fq" ]; then
    echo "File ${DIR}/${JB}/${JB}_A1P1.M${mmatchnum}.fq not found!"
    exit 1
fi

if [ ! -e "${DIR}/${JB}/${JB}_A1P2.M${mmatchnum}.fq" ]; then
    echo "File ${DIR}/${JB}/${JB}_A1P2.M${mmatchnum}.fq not found!"
    exit 1
fi

# Check if subset has been indicated and check for mapping file either way.
if [ -n "$OPTION_O" ]; then
    if [ -e "${DIR}/${JB}/${JB}_map.${OPTION_O}.txt" ]; then
        MAPFILE="${DIR}/${JB}/${JB}_map.${OPTION_O}.txt"
        ODIR="${DIR}/${JB}_${OPTION_O}_output"
        JB2=${JB}_${OPTION_O}
    else
        echo "File ${DIR}/${JB}/${JB}_map.${OPTION_O}.txt not found!"
        exit 1
    fi
else
    if [ -e "${DIR}/${JB}/${JB}_map.txt" ]; then
        MAPFILE="${DIR}/${JB}/${JB}_map.txt"
        ODIR="${DIR}/${JB}_output"
        JB2=${JB}
    else
        echo "File ${DIR}/${JB}/${JB}_map.txt not found!"
        exit 1
    fi
fi
echo "File ${MAPFILE} is the mapping file you have selected."

echo "Do you want to continue processing for $JB? (yes/no)"
read response
while [[ "$response" != "yes" && "$response" != "no" ]]; do
    echo "Invalid response. Please enter 'yes' or 'no'."
    read response
done

if [ "$response" == "no" ]; then
    echo "Skipping $JB."
    continue
fi

# initiate log
timestamp="$(date +"%Y%m%d_%H:%M:%S")"
output_file="${DIR}/part2_${timestamp}.log"
exec > "$output_file" 2>&1

# create header for log file
echo " - -- --- ---- ---- --- -- -"
echo "Log file for Part 2 of the Microbiome Pipeline. Processing the following arguments:
Working Directory: ${DIR}
Data ID: ${JB}
Subset if specified: ${OPTION_O}
Mismatches if specified: ${mmatchnum}
 - -- --- ---- ---- --- -- -"

# make the output directory
mkdir -vp ${ODIR}

echo
echo " - -- --- ---- ---- --- -- -"
echo "Generating fastq information"
echo " - -- --- ---- ---- --- -- -"

#This code helps you to understand the make-up of your reads
# Extract fastq info and create individual directories for each JB inside fastq_info
mkdir -vp "${ODIR}/fastq_info"
usearch -fastx_info "${DIR}/${JB}/${JB}_A1P1.M${mmatchnum}.fq" -quiet -output "${ODIR}/fastq_info/${JB2}_A1P1.M${mmatchnum}.txt" || (echo "Error generating fastq info for ${JB}." && exit 1)

# Create barcodes.fa file for each JB
grep -P "^[A-Z]" "${MAPFILE}" | awk '{print ">"$1"\n"$2}' > "${ODIR}/${JB2}_barcodes.fa"
[[ ! -e "${ODIR}/${JB2}_barcodes.fa" ]] && echo "File ${ODIR}/${JB2}_barcodes.fa was not generated!" && exit 1

# Demultiplexing
usearch -fastx_demux "${DIR}/${JB}/${JB}_A1P1.M${mmatchnum}.fq" -quiet -index "${DIR}/${JB}/${JB}_A1P2.M${mmatchnum}.fq" -barcodes "${ODIR}/${JB2}_barcodes.fa" -fastqout "${ODIR}/${JB2}_A1P1.M${mmatchnum}.fq" && echo "File ${ODIR}/${JB2}_A1P1.M${mmatchnum}.fq has been generated." || (echo "Error: File ${ODIR}/${JB2}_A1P1.M${mmatchnum}.fq was not generated." && exit 1)

cp "${DIR}/${JB}/${JB}_A1P2.M${mmatchnum}.fq" "${ODIR}/${JB2}_A1P2.M${mmatchnum}.fq"
[[ ! -e "${ODIR}/${JB2}_A1P2.M${mmatchnum}.fq" ]] && echo "File "${ODIR}/${JB2}_A1P2.M${mmatchnum}.fq" was not generated!" && exit 1

echo
echo "General information about your data has been saved in ${ODIR}/fastq_info/${JB2}_A1P1.M${mmatchnum}.txt."  | tee /dev/tty

# At this point, you need to decide if you should truncate the reads. 
# This will print out some stats showing how many of the reads will remain 
# at each trunc length once you filter for quality. 
# (The more you trim, the more reads remain after filtering.)

# First, use eestats to determine read qualities.
# You choose a "min" truncated length and assign it to the variable STARTAT. 
# You'll get a printout of all the lengths above that. 
# This will depend on the length of your reads. 220 is a good STARTAT for 301 bp reads and 140 is good for 151 bp reads.

echo
echo " - -- --- ---- ---- --- -- -
Generating stats about the potential effects trimming and filtering will have on your reads
 - -- --- ---- ---- --- -- -"

# Calculate the maximum length of the first 100 reads
MAX_LENGTH=$(head -n 400 ${ODIR}/${JB2}_A1P1.M${mmatchnum}.fq | awk '{if(NR%4==2) print length($1)}' | sort -nr | head -n 1)

# Set STARTAT based on the calculated max length
if (( MAX_LENGTH > 250 )); then
    STARTAT=220; # For 301 bp reads
else
    STARTAT=140; # For 151 bp reads
fi

INC=1; # Increment value

[[ -e ${ODIR}/all_eestats${STARTAT},${INC}.txt ]] && rm ${ODIR}/all_eestats${STARTAT},${INC}.txt

usearch -fastq_eestats2 ${ODIR}/${JB2}_A1P1.M${mmatchnum}.fq -quiet -output "${ODIR}/${JB2}.M${mmatchnum}_eestats.start_${STARTAT}.inc_${INC}.txt" -length_cutoffs ${STARTAT},*,${INC}

# final message - what is next
echo
echo " - -- --- ---- ---- --- -- -
Stats ready. Please view and select your trim length accordingly.
Stats have been saved as: ${ODIR}/${JB2}.M${mmatchnum}_eestats.start_${STARTAT}.inc_${INC}.txt

Part 3 will combine as many of your data files as you specify into a single ASV table. 
You will indicate your trim length under -l and an output file for your results under -o.

-j will be all the samples you want to include by output directory name 
("${ODIR}" is the one you selected for this run, but you may have others.)

For example, next you might run: 
part3.sh -d ${DIR} -j "${ODIR}",another_ID,another_ID2 -l <trim length> -o final_output_directory
 - -- --- ---- ---- --- -- -" | tee /dev/tty


