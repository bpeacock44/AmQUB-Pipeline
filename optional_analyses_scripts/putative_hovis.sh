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
exec > "$output_file"
exec 2> >(tee -a "$output_file" >&2)

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
    # create "top" versions - only the best scoring match for each taxonomic unit
    awk '!seen[$1]++' "${OUT}" > "${WDIR}/finding_more_hovis/our_6.top.blastout"
    # filter ASVs for those with taxonomic predictions indicative of possible hovis (i.e. they fit various criteria)
    filtering_tus_for_hovis.py "${WDIR}/finding_more_hovis/our_6.top.blastout" \
        "${WDIR}/classifier_output/taxonomy.tsv" \
        "${WDIR}/blast/tax_assignments.txt" \
        "${WDIR}/finding_more_hovis/solid_hovi_${typ}s.txt" \
        "${WDIR}/finding_more_hovis/putative_hovi_${typ}s.txt"
        # solid hovis are those that were specifically hyalorbilia oviparasitica/etc.
        # putative hovis match various criteria making them more likely to be hovi.
    # Extract all taxonomic units that matched up with our 6 fungi
    cut -f1 "${WDIR}/finding_more_hovis/our_6.top.blastout" > "${WDIR}/finding_more_hovis/our_6.ids.txt"
    # clean up the blast output so sizes aren't appended to the headers
    awk -F'\t' '{sub(/_[^_]*$/, "", $1)}1' OFS='\t' "${WDIR}/blast/filtered.blastout" > "${WDIR}/blast/filtered.clean.blastout"
    # Filter BLAST results to only include taxonomic units that matched to our 6 fungi 
    awk 'NR==FNR {ids[$1]; next} $1 in ids || $0 ~ /^#/ {print}' "${WDIR}/finding_more_hovis/our_6.ids.txt" "${WDIR}/blast/filtered.clean.blastout" > "${WDIR}/blast/filtered.clean.our6.blastout"
    # Filter blast to find all taxonomic units that meet our 2nd level of criteria (3 conditions)
    final_hovi_filter.py "${WDIR}/blast/filtered.clean.our6.blastout" "${WDIR}/finding_more_hovis/blast_summary_output.tsv" 
    # find all putative hovis that also meet the critera for the final hovi filter and extract seqs
    grep -Fw -f "${WDIR}/finding_more_hovis/putative_hovi_${typ}s.txt" "${WDIR}/finding_more_hovis/blast_summary_output.tsv" > "${WDIR}/finding_more_hovis/final_putative_hovis.txt"
    cat "${WDIR}/finding_more_hovis/solid_hovi_${typ}s.txt" >> "${WDIR}/finding_more_hovis/final_putative_hovis.txt"
    sort "${WDIR}/finding_more_hovis/final_putative_hovis.txt" | uniq > "${WDIR}/finding_more_hovis/temp.txt" && mv "${WDIR}/finding_more_hovis/temp.txt" "${WDIR}/finding_more_hovis/final_putative_hovis.txt"
    seqkit grep -f "${WDIR}/finding_more_hovis/final_putative_hovis.txt" $INFASTA > "${WDIR}/finding_more_hovis/final_putative_hovis.fa"
    cat ~/shared/mbio_pipeline_files/control_hya_seqs_for_tree.fa >> "${WDIR}/finding_more_hovis/final_putative_hovis.fa"
    mafft --phylipout --thread ${NUMTHREADS} --auto "${WDIR}/finding_more_hovis/final_putative_hovis.fa" > "${WDIR}/finding_more_hovis/final_putative_hovis.phy"

    # create new output files with putative hovis annotated
    mkdir -vp "${WDIR}/finding_more_hovis/new_output_files"
    
    # Function to add '-PH' to each string in the list, unconditionally
    add_ph_to_strings() {
        local list_file="$1"    # Input file with strings (e.g., final_putative_hovis.txt)
        local target="$2"       # A single file or a glob pattern (e.g., *.txt or a specific file)
        local output_dir="$3"   # Directory where modified files will be saved
    
        # Expand the glob pattern to an array of files
        eval "files=($target)"
    
        # Loop through all the files in the target
        for file in "${files[@]}"; do
            # Copy the original file to the temporary file
            cp "$file" "${output_dir}/$(basename "$file")"
    
            # Loop over each string from the list file
            while IFS= read -r otu; do
                sed -i "s/\b$otu\b/&-PH/g" "${output_dir}/$(basename "$file")"
            done < "$list_file"
        done
    }
    
    # Process all *.txt files
    add_ph_to_strings "${WDIR}/finding_more_hovis/final_putative_hovis.txt" \
                        "${WDIR}/*.txt" \
                        "${WDIR}/finding_more_hovis/new_output_files"
    
    # Process a specific pattern *.s.fa files
    add_ph_to_strings "${WDIR}/finding_more_hovis/final_putative_hovis.txt" \
                        "${WDIR}/${typ}s.fa" \
                        "${WDIR}/finding_more_hovis/new_output_files"
    
    add_ph_to_strings "${WDIR}/finding_more_hovis/final_putative_hovis.txt" \
                        "${WDIR}/Detailed_Informational_otu_Table.tsv" \
                        "${WDIR}/finding_more_hovis/new_output_files"
    
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
    
    # create new version of Detailed file
    # Define the file paths
    input_file="${WDIR}/finding_more_hovis/new_output_files/Detailed_Informational_otu_Table.tsv"  # Input TSV file
    putative_file="${WDIR}/finding_more_hovis/final_putative_hovis.txt"  # List of putative hovis IDs
    
    # Step 1: Add a new column "putative_hovi" with "no" for all rows except the first two
    awk '
    BEGIN {FS="\t"; OFS="\t"}
    
    # First row: Add "putative_hovi" as the second column
    NR==1 {
        print $1, "putative_hovi", $0;
    }
    
    # Second row: Add "nd" as the second column
    NR==2 {
        print $1, "nd", $0;
    }
    
    # All other rows: Add "no" as the second column
    NR > 2 {
        print $1, "no", $0;
    }
    ' "$input_file" > temp_file.tsv
    
    # Step 2: Read the putative hovis IDs into a hash map and replace "no" with "yes" where appropriate
    declare -A putative_hovis_ids
    while read -r id; do
        id_with_suffix="${id}-PH"  # Append "-PH" to the ID
        putative_hovis_ids["$id_with_suffix"]=1  # Store ID in the hash map
    done < "$putative_file"
    
    # Step 3: Update the rows in temp_file.tsv based on the hash map
    # Convert the associative array keys into a space-separated string
    putative_hovis_list=$(printf "%s " "${!putative_hovis_ids[@]}")
    
    # Use awk with a direct variable
    awk -v ids="$putative_hovis_list" '
    BEGIN {
        FS = OFS = "\t";
        split(ids, putative_hovis_ids, " ");
        for (i in putative_hovis_ids) {
            lookup[putative_hovis_ids[i]] = 1;
        }
    }
    NR == 1 { $2 = "putative_hovi"; print; next }  # First row header change
    NR == 2 { $2 = "nd"; print; next }             # Second row placeholder
    NR > 2 { 
        if ($3 in lookup) {
            $2 = "yes";
        } else {
            $2 = "no";
        }
        print;
    }
    ' temp_file.tsv > "$input_file"
    
    # Cleanup
    rm temp_file.tsv
done
