#!/bin/bash
# See mbio_tutorial.md for further guidance!

### USAGE ###
# This script accepts two types of input:
# 1. Direct command-line arguments:
#    putative_hovis.sh -i output_dir

## Required Flags
# -i: The output directory generated in part 3 that you also ran part 4 on.
# -n: The number of threads available for BLAST
# -6: A fasta file of sequences for our 6 hyalorbilia fungi

## Optional Flags:
# --nostrategy1: Don't process the output from strategy 1. (This is the default strategy.)
# --strategy2: Process the output from strategy 2.
# --strategy3: Process the output from strategy 3.

# 2. A parameter template file:
#    AmQUB_part4.sh params.csv
#    Where params.csv contains the following rows, comma delimited, no white space between.
#    The labels at the beginning of each row should be the same as below.
#    The following shows what all possible rows:

#        Part 3 Output Folder To Process,output_dir
#        Number of Threads Available,256
#        Fasta of Our 6,~/shared/mbio_pipeline_files/our_6_db
#        Skip Strategy 1,true
#        Process Strategy 2,true
#        Process Strategy 3,true
#
#Any optional line can be left out of the file if you with to use default settings.

# Set error handling
set -e  # Exit on any error

# Function to parse parameter file
parse_parameter_file() {
    local param_file="$1"
    while IFS= read -r line; do
        case "$line" in
            "Part 3 Output Folder To Process,"*) output_dir="${line#*,}" ;;
            "Number of Threads Available,"*) NUMTHREADS="${line#*,}" ;;
            "Fasta of Our 6,"*) our_6="${line#*,}" ;;
            "Skip Strategy 1,"*) skSTR1="${line#*,}" ;;
            "Process Strategy 2,"*) STR2_simp="${line#*,}" ;;
            "Process Strategy 3,"*) STR3_simp="${line#*,}" ;;
        esac
    done < "$param_file"
}

# Initialize default values
skSTR1=false
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

Putative Hovis Add-On

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
            -n|--numthreads) NUMTHREADS="$2"; shift 2 ;;
            -6|--our_6) our_6="$2"; shift 2 ;;
            --nostrategy1) skSTR1=true; shift ;;
            --strategy2) STR2_simp=true; shift ;;
            --strategy3) STR3_simp=true; shift ;;
            *) echo "Unknown option: $1" >&2; exit 1 ;;
        esac
    done
fi

# Check for mandatory arguments
if [ -z "$output_dir" ] || [ -z "$NUMTHREADS" ] || [ -z "$our_6" ]  ; then
    echo "Usage: $0 -i <output directory from part 3> -n <number of threads for blast> -6 <our 6 fungi blast database> [--nostrategy1 --strategy2 --strategy3]"
    exit 1
fi

# Make sure output directory and blast database exists
if [[ ! -d "$output_dir" ]]; then
    echo "Error: Output directory '$output_dir' does not exist."
    exit 1
fi

if [[ ! -e "$our_6" ]]; then
    echo "Error: Fasta file '$our_6' not found."
    exit 1
fi

# Make sure numthreads is a valid positive integer.
if ! [[ "$NUMTHREADS" =~ ^[0-9]+$ ]]; then
    echo "Error: Number of threads must be a positive integer."
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
    DIRS+=("${output_dir}/${typ}s/STRATEGY3/${typ}s")
fi

if [ "$skSTR1" == true ]; then
    if [[ "$STR2_simp" == false && "$STR3_simp" == false ]]; then
        echo "Error: You chose to skip Strategy 1 but did not specify Strategy 2 or 3 processing."
        exit 1
    fi
fi

# initiate log
timestamp="$(date +"%y%m%d_%H:%M")"
output_file="${output_dir}/putative_hovis_${timestamp}.log"
exec > >(tee -ia "$output_file") 2>&1

# log header
echo "Log file for Putative Hovis Finder. Processed the following arguments:
Output directory to be processed: ${output_dir}
Threads to use: ${NUMTHREADS}
Location of fasta file for our 6 fungi: ${our_6}" | tee /dev/tty
if [ "$skSTR1" != false ]; then 
    echo "Strategy 1 (default) results are not being processed."
fi

if [ "$STR2_simp" != false ]; then 
    echo "Strategy 2 results will be processed."
fi

if [ "$STR3_simp" != false ]; then 
    echo "Strategy 3 results will be processed."
fi

echo " - -- --- ---- ---- --- -- -"

# detect if classifier assignments are present and use them in final summary file if generated
if [ -d "${output_dir}/${typ}s/classifier_output" ]; then
    CP=true 
fi

base_name="${our_6%.*}"  # This removes the last extension
our_6_db="$output_dir/${base_name}_db/${base_name}"
mkdir -p "$output_dir/${base_name}_db"
makeblastdb -in "$our_6" -dbtype nucl -out "${our_6_db}" || { 
    echo "Failed to create BLAST database."; 
    exit 1; 
}

# run putative hyalorbilia steps on each directory
OPTS="qseqid sseqid pident length mismatch evalue bitscore staxids stitle qcovs"
for WDIR in "${DIRS[@]}"; do
    # create directory
    mkdir -p "${WDIR}/finding_more_hovis"
    # align taxonomic units to our fungi
    INFASTA="${WDIR}/${typ}s.fa"
    OUT="${WDIR}/finding_more_hovis/our_6.blastout"
    blastn -task blastn -db $our_6_db -query $INFASTA -max_target_seqs 100 -evalue 0.001 -num_threads $NUMTHREADS -outfmt "6 $OPTS" > ${OUT}
    
    # generate "chunks" of the matching reads to make trees from
    matching_fas="${WDIR}/finding_more_hovis/6_matching_otus.fa"
    awk '{print $1}' "$OUT" | sort | uniq > "${WDIR}/finding_more_hovis/ids.txt"
    # pull sequences
    seqkit grep -f "${WDIR}/finding_more_hovis/ids.txt" "${INFASTA}" > "${matching_fas}"
    rm -rf "${WDIR}/finding_more_hovis/6_matching_otus_chunks"
    # split into 50-sequence chunks
    seqkit split2 -s 50 -O "${WDIR}/finding_more_hovis/6_matching_otus_chunks" "${matching_fas}"

    PH_via_trees.py "${WDIR}/finding_more_hovis/6_matching_otus_chunks" \
        --ref "/rhome/bpeacock/shared/mbio_pipeline_files/PH_control_seqs.fa" 

    # create final_putative_hovis.fa
    cat "${WDIR}/finding_more_hovis/6_matching_otus_chunks/"*_PH.fa "/rhome/bpeacock/shared/mbio_pipeline_files/PH_control_seqs.fa" > "${WDIR}/finding_more_hovis/PH_with_controls.fa"


    {
        first=1
        for f in "${WDIR}/finding_more_hovis/6_matching_otus_chunks/"*_PH.txt; do
            (( first )) || printf '\n'
            cat "$f"
            first=0
        done
    } | sort > "${WDIR}/finding_more_hovis/PH.txt"
    awk 'NF' "${WDIR}/finding_more_hovis/PH.txt" > tmp && mv tmp "${WDIR}/finding_more_hovis/PH.txt"

    {
        first=1
        for f in "${WDIR}/finding_more_hovis/6_matching_otus_chunks/"*_PPH.txt; do
            (( first )) || printf '\n'
            cat "$f"
            first=0
        done
    } | sort > "${WDIR}/finding_more_hovis/PPH.txt"
    awk 'NF' "${WDIR}/finding_more_hovis/PPH.txt" > tmp && mv tmp "${WDIR}/finding_more_hovis/PPH.txt"

    # create new output files with putative hovis annotated
    rm -rf "${WDIR}/finding_more_hovis/new_output_files"
    mkdir -vp "${WDIR}/finding_more_hovis/new_output_files"
    
    add_ph_to_strings() {
        local list_file="$1"
        local output_dir="$2"
        local string="$3"
        shift 3
    
        local target_files=("$@")
    
        mkdir -p "$output_dir"
    
        for file in "${target_files[@]}"; do
            local out_file="${output_dir}/$(basename "$file")"
    
            [[ -f "$out_file" ]] || cp "$file" "$out_file"
    
            awk -F'\t' -v OFS='\t' -v list="$list_file" -v rep="$string" '
            BEGIN {
                while ((getline otu < list) > 0) {
                    if (otu != "") {
                        map[otu] = 1
                    }
                }
            }
            {
                for (i = 1; i <= NF; i++) {
                    if ($i in map) {
                        $i = rep "-" $i
                    }
                }
                print
            }
            ' "$out_file" > "${out_file}.tmp" && mv "${out_file}.tmp" "$out_file"
    
        done
    }
        
    strings=(PH PPH)
    
    for S in "${strings[@]}"; do
        add_ph_to_strings \
            "${WDIR}/finding_more_hovis/${S}.txt" \
            "${WDIR}/finding_more_hovis/new_output_files" \
            "$S" \
            ${WDIR}/*.txt
    
        add_ph_to_strings \
            "${WDIR}/finding_more_hovis/${S}.txt" \
            "${WDIR}/finding_more_hovis/new_output_files" \
            "$S" \
            "${WDIR}/${typ}s.fa"
    
        add_ph_to_strings \
            "${WDIR}/finding_more_hovis/${S}.txt" \
            "${WDIR}/finding_more_hovis/new_output_files" \
            "$S" \
            "${WDIR}/Detailed_Informational_otu_Table.tsv"
    done
    
    #source for bash helper functions
    source "qiime_shell_helper_functions.sh"
            
    #convert all to biom and qza
    rm -rf "${WDIR}/finding_more_hovis/new_output_files"/*00* "${WDIR}/finding_more_hovis/new_output_files"/*01* "${WDIR}/finding_more_hovis/new_output_files"/*02*
    for F in "${WDIR}/finding_more_hovis/new_output_files"/*03*.txt; do    
        biom convert -i "$F" -o "${F%.txt}.biom" --table-type="OTU table" --to-hdf5
        qiime tools import \
        --input-path "${F%.txt}.biom" \
        --type 'FeatureTable[Frequency]' \
        --output-path "${F%.txt}.qza" 
    done
    
    # Created detailed informational otu table xlsx version
    input_tsv="${WDIR}/finding_more_hovis/new_output_files/Detailed_Informational_otu_Table.tsv"
    output_xlsx="${WDIR}/finding_more_hovis/new_output_files/Detailed_Informational_otu_Table.xlsx"
    
    python3 -c "
import pandas as pd
df = pd.read_csv('${input_tsv}', sep='\t')
df.to_excel('${output_xlsx}', index=False)
print('Wrote ${output_xlsx}')
"

    echo "Putative Hovi Pipeline Complete"
done
