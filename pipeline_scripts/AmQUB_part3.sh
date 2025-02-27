#!/bin/bash
# See mbio_tutorial.md for further guidance!

### USAGE ###
# This script accepts two types of input:
# 1. Direct command-line arguments:
#    AmQUB_part3.sh --in ID1_raw_output+ID2_raw_output --len 200 --out output_dir --al UPARSE
#    AmQUB_part3.sh --in ID1_raw_output --len 200 --out output_dir --al UNOISE3 \
#        --min 8 --ag yes --alpha 1.5
#    AmQUB_part3.sh --in ID1_raw_output+ID2_raw_output --len 200 --out output_dir --al UPARSE \
#        --fr yes --fe yes --fm 0.5 --fd 2 --ug yes --us 2 \
#        --min 8 --ag yes --pid 0.98 --map merged_map.txt \
#        --col Soil_Type+Tissue --pre IDx_OTUs.fa \
#        --tblid 0.98 --un yes --rmf yes              

## Required Flags
# --in Part 2 Output Folders To Process Together (list with "+" between each folder)
# --len Trim Length
# --out Output Folder
# --al UPARSE or UNOISE3

## Optional Flags:
# --fr FILTERING-Generation of Removed Reads File
# --fe FILTERING-EE Value Appended to IDs
# --fm FILTERING-MaxEE Value
# --fd FILTERING-Discard if there are > n Ns in the read
# --ug UNIQUES-Generation of Sequence Binning Report
# --us UNIQUES-Set minuniquesize
# --min ALGORITHM-Set Min Size
# --ag ALGORITHM-Generation of Processing Report
# --pid UPARSE-Set Percent Identity Threshold
# --alpha UNOISE-Set Alpha Parameter
# --map STR2-Mapping File 
# --col STR2-Treatment Column (can be +-delimited list if multiple columns to be processed)
# --pre STR3-Pre-existing OTUsASVs 
# --tblid TABLE-Minimum fractional ID
# --un-TABLE-Generation of Fastq of Unassigned Seqs
# --rmf-TABLE-Generation of Read Mapping File

# 2. A parameter template file:
#    AmQUB_part3.sh params.csv
#    Where params.csv contains the following rows, comma delimited, no white space between.
#    The labels at the beginning of each row should be the same as below (i.e. the part before the first comma)
#    The following shows what might be used for the 3rd example in the command-line example above.

#        Part 2 Output Folders To Process Together,ID1_raw_output+ID2_raw_output
#        Trim Length,200
#        Output Folder,output_dir
#        UPARSE or UNOISE3,UPARSE
#        FILTERING-Generation of Removed Reads File,true
#        FILTERING-EE Value Appended to IDs,true
#        FILTERING-MaxEE Value,0.5
#        FILTERING-Discard if there are > n Ns in the read,2
#        UNIQUES-Generation of Sequence Binning Report,true
#        UNIQUES-Set minuniquesize,2
#        ALGORITHM-Set Min Size,8
#        ALGORITHM-Generation of Processing Report,true
#        UPARSE-Set Percent Identity Threshold,0.98
#        UNOISE-Set Alpha Parameter,NA (this line wouldn't be included)
#        STR2-Mapping File,merged_map.txt
#        STR2-Treatment Column,Soil_Type+Tissue
#        STR3-Pre-existing OTUsASVs,IDx_OTUs.fa
#        TABLE-Minimum fractional ID,0.98
#        TABLE-Generation of Fastq of Unassigned Seqs,true
#        TABLE-Generation of Read Mapping File,true
#
#All optional parameters can be left out of the param file.

set -e  # Exit on error

# Function to check if a required label exists in the file
check_label_present() {
    local label="$1"
    if ! grep -iq "^$label," "$2"; then
        echo "Error: Missing expected label '$label' in parameter file '$2'." | tee /dev/tty
        exit 1
    fi
}

# Function to parse parameter file and store inputs
parse_file_input() {
    for label in "Part 2 Output Folders To Process Together" "Trim Length" "Output Folder" "UPARSE or UNOISE3"; do
        check_label_present "$label" "$1"
    done

    while IFS= read -r line; do
        case "$line" in
            "Part 2 Output Folders To Process Together,"*) IFS='+' read -ra INP <<< "${line#*,}" ;;
            "Trim Length,"*) LEN="${line#*,}" ;;
            "Output Folder,"*) OUTDIR="${line#*,}" ;;
            "UPARSE or UNOISE3,"*) ALGORITHM="${line#*,}" ;;
            "FILTERING-Generation of Removed Reads File,"*) FR="${line#*,}" ;;
            "FILTERING-EE Value Appended to IDs,"*) FE="${line#*,}" ;;
            "FILTERING-MaxEE Value,"*) FM="${line#*,}" ;;
            "FILTERING-Discard if there are > n Ns in the read,"*) FD="${line#*,}" ;;
            "UNIQUES-Generation of Sequence Binning Report,"*) UG="${line#*,}" ;;
            "UNIQUES-Set minuniquesize,"*) US="${line#*,}" ;;
            "ALGORITHM-Set Min Size,"*) MIN="${line#*,}" ;;
            "ALGORITHM-Generation of Processing Report,"*) AG="${line#*,}" ;;
            "UPARSE-Set Percent Identity Threshold,"*) UP="${line#*,}" ;;
            "UNOISE-Set Alpha Parameter,"*) ALPHA="${line#*,}" ;;
            "STR2-Mapping File,"*) MAPF="${line#*,}" ;;
            "STR2-Treatment Column,"*) IFS='+' read -ra COL <<< "${line#*,}" ;;
            "STR3-Pre-existing OTUsASVs,"*) PRE="${line#*,}" ;;
            "TABLE-Minimum fractional ID,"*) TBLID="${line#*,}" ;;
            "TABLE-Generation of Fastq of Unassigned Seqs,"*) UN="${line#*,}" ;;
            "TABLE-Generation of Read Mapping File,"*) RMF="${line#*,}" ;;
        esac
    done < "$1"
}

# Ensure at least one dataset is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 --in <comma-separated output directories> --len <trim length> --out <output dir> --al <UPARSE|UNOISE3> [optional flags]"
    exit 1
fi

    echo " - -- --- ---- ---- --- -- -
 █████╗             ██████╗ ██╗   ██╗██████╗ 
██╔══██╗████╗ ████╗██║   ██╗██║   ██║██╔══██╗
███████║██╔████╔██║██║   ██║██║   ██║██████╔╝
██╔══██║██║╚██╔╝██║██║   ██║██║   ██║██╔══██╗
██║  ██║██║ ╚═╝ ██║╚██████╔╝╚██████╔╝██████╔╝
╚═╝  ╚═╝╚═╝     ╚═╝ ╚════██╗ ╚═════╝ ╚═════╝  
PART 3:OTU/ASV Generation╚═╝
 - -- --- ---- ---- --- -- -"

# Set default values for optional parameters if not set
MIN="${MIN:-false}"
FR="${FR:-false}"
FE="${FE:-false}"
FM="${FM:-1.0}"
FD="${FD:-false}"
UG="${UG:-false}"
US="${US:-false}"
AG="${AG:-false}"
UP="${UP:-0.97}"
ALPHA="${ALPHA:-false}"
MAPF="${MAPF:-false}"
COL="${COL:-false}"
PRE="${PRE:-false}"
TBLID="${TBLID:-false}"
UN="${UN:-false}"
RMF="${RMF:-false}"

# Check if the first argument is a parameter file
if [[ -f "$1" ]]; then
    echo "Reading parameters from file: $1
 - -- --- ---- ---- --- -- -"
    parse_file_input "$1"
else
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --in) IFS='+' read -ra INP <<< "$2"; shift ;;
            --len) LEN="$2"; shift ;;
            --out) OUTDIR="$2"; shift ;;
            --al) ALGORITHM="$2"; shift ;;
            --fr) FR="$2"; shift ;;
            --fe) FE="$2"; shift ;;
            --fm) FM="$2"; shift ;;
            --fd) FD="$2"; shift ;;
            --ug) UG="$2"; shift ;;
            --us) US="$2"; shift ;;
            --min) MIN="$2"; shift ;;
            --ag) AG="$2"; shift ;;
            --pid) UP="$2"; shift ;;
            --alpha) ALPHA="$2"; shift ;;
            --map) MAPF="$2"; shift ;;
            --col) IFS='+' read -ra COL <<< "$2"; shift ;;
            --pre) PRE="$2"; shift ;;
            --tblid) TBLID="$2"; shift ;;
            --un) UN="$2"; shift ;;
            --rmf) RMF="$2"; shift ;;
            *) echo "Unknown option $1"; exit 1 ;;
        esac
        shift
    done
fi

# Check for mandatory arguments
if [ -z "$INP" ] || [ -z "$LEN" ] || [ -z "$OUTDIR" ] || [ -z "$ALGORITHM" ]; then
    echo "Usage: $0 --in <in directories> --len <trim length> --out <output dir> --al <UPARSE|UNOISE3> [optional flags]"
    exit 1
fi

# Check if $ALGORITHM is defined as either UNOISE3 or UPARSE
if [[ "$ALGORITHM" != "UNOISE3" && "$ALGORITHM" != "UPARSE" ]]; then
    echo "Error: ALGORITHM must be defined as either 'UNOISE3' or 'UPARSE'."
    exit 1
fi

if { [[ -z "$MAPF" ]] && [[ -n "$COL" ]]; } || { [[ -n "$MAPF" ]] && [[ -z "$COL" ]]; }; then
    echo "Error: If you want to use Strategy 2 to create taxonomic units, you must indicate BOTH a mapping file AND a column name."
    exit 1
fi

# show required input files
for IN in "${INP[@]}"; do
    matching_files=($(find -L "$IN" -type f -name "*demux*" -not -name "*bp*"))
    # Check if more than one file was found
    if [[ ${#matching_files[@]} -gt 1 ]]; then
        echo "Error: More than one demux file found in $IN. You may have run parts 1 and 2 with two different mismatch levels. This is fine, but you need to choose one. Please delete the demux.fq file in ${IN} you don't want to use. The number of mismatches is indicated after the M. (So for example, ID.M2.demux.fq would have 2 mismatches and ID.M0.demux.fq would have 0.)"
        exit 1
    fi
done

######### BEGIN PROCESSING
# make the output directory
if ! mkdir -p "${OUTDIR}"; then
    echo "Error: Failed to create directory '$OUTDIR'."
    exit 1
fi

# initiate log
timestamp="$(date +"%y%m%d_%H:%M")"
output_file="${OUTDIR}/part3_${timestamp}.log"
exec > "$output_file"
exec 2> >(tee -a "$output_file" >&2)

# log header
echo " - -- --- ---- ---- --- -- -"
echo "Log file for Part 3 of the Microbiome Pipeline. Processing the following arguments:" | tee /dev/tty

# Function to conditionally print variables
print_arg() {
  local label="$1"
  local value="$2"
  if [[ "$value" != "false" ]]; then
    echo "${label} ---> ${value}" | tee /dev/tty
  fi
}

# Print each variable only if it’s not "no"
print_arg "Input Directories" "${INP[*]}"
print_arg "Trim Length" "$LEN"
print_arg "Output Directory" "$OUTDIR"
print_arg "Algorithm" "$ALGORITHM"
print_arg "FILTERING - Removed Reads File" "$FR"
print_arg "FILTERING - EE Value in IDs" "$FE"
print_arg "FILTERING - MaxEE Value" "$FM"
print_arg "FILTERING - Max Ns" "$FD"
print_arg "UNIQUES - Sequence Binning Report" "$UG"
print_arg "UNIQUES - Min Unique Size" "$US"
print_arg "ALGORITHM - Min Size" "$MIN"
print_arg "ALGORITHM - Processing Report" "$AG"
# Conditionally print UPARSE only if the algorithm is NOT UNOISE3
if [[ "$ALGORITHM" != "UNOISE3" ]]; then
  print_arg "UPARSE - Percent Identity Threshold" "$UP"
fi
print_arg "UNOISE - Alpha Parameter" "$ALPHA"
print_arg "STR2 - Mapping File" "${MAPF[*]}"
print_arg "STR2 - Treatment Column" "${COL[*]}"
print_arg "STR3 - Pre-existing OTUs/ASVs" "$PRE"
print_arg "TABLE - Minimum fractional ID" "$TBLID"
print_arg "TABLE - Generation of Fastq of Unassigned Seqs" "$UN"
print_arg "TABLE - Generation of Read Mapping File" "$RMF"

rm -f "${OUTDIR}/combined.fq"
rm -f "${OUTDIR}/filtered.fa" 

# Trimming and Filtering Reads
for IN in "${INP[@]}"; do
    matching_files=($(find -L "$IN" -type f -name "*demux*" -not -name "*bp*"))
    # Assign the single matching file to FQ if exactly one file is found
    FQ="${matching_files[0]}"
    FQ_base="${FQ%.fq}"
    # truncate reads at LEN
    echo
    echo " - -- --- ---- ---- --- -- -"
    echo " - RUNNING: usearch -fastx_truncate ${FQ} -quiet -trunclen ${LEN} -fastqout ${FQ_base}_${LEN}bp.fq"
    usearch -fastx_truncate "${FQ}" -quiet -trunclen "${LEN}" -fastqout "${FQ_base}_${LEN}bp.fq"
    # create combined read file
    cat "${FQ_base}_${LEN}bp.fq" >> "${OUTDIR}/combined.fq"
    #maxee quality filtering of demultiplexed/truncated fq files (*** keep THREADS=1 for repeatability ***)
    CMD=("usearch" "-threads" "1" "-fastq_filter" "${FQ_base}_${LEN}bp.fq" "-quiet" "-fastq_maxee" "${FM}" "-fastaout" "${FQ_base}_${LEN}bp.filtered.fa")
    [[ "$FR" == "true" ]] && CMD+=("-fastqout_discarded" "${FQ_base}_${LEN}bp_discarded.fq")
    [[ "$FE" == "true" ]] && CMD+=("-fastq_eeout")
    [[ "$FD" == "true" ]] && CMD+=("-fastq_maxns" "${FD}")
    echo
    echo " - -- --- ---- ---- --- -- -"
    echo " - RUNNING: ${CMD[@]}"
    "${CMD[@]}"

    # Create filtered.fa read file
    cat "${FQ_base}_${LEN}bp.filtered.fa" >> "${OUTDIR}/filtered.fa"
done

#find unique sequences
CMD=("usearch" "-fastx_uniques" "${OUTDIR}/filtered.fa" "-quiet" "-fastaout" "${OUTDIR}/uniques.fa" "-sizeout" "-relabel" "Uniq")
[[ "$UG" == "true" ]] && CMD+=("-tabbedout" "${OUTDIR}/uniques_binning_report.txt")
[[ "$US" == "true" ]] && CMD+=("-minuniquesize" "${US}")
echo " "
echo " - -- --- ---- ---- --- -- -"
echo " - RUNNING: ${CMD[@]}"
"${CMD[@]}"

# sort by size in preparation for algorithm
echo " "
echo " - -- --- ---- ---- --- -- -"
echo " - RUNNING: usearch -sortbysize ${OUTDIR}/uniques.fa -quiet -fastaout ${OUTDIR}/temp.fa -minsize 1"
usearch -sortbysize "${OUTDIR}/uniques.fa" -quiet -fastaout "${OUTDIR}/temp.fa" -minsize 1
mv "${OUTDIR}/temp.fa" "${OUTDIR}/uniques.fa"

# create otus/asvs using algorithm of choice
if [[ "$ALGORITHM" == "UNOISE3" ]]; then
    mkdir -vp "${OUTDIR}/asvs"
    CMD=("usearch" "-unoise3" "${OUTDIR}/uniques.fa" "-quiet" "-zotus" "${OUTDIR}/asvs/asvs.fa")
    [[ "$MIN" == "true" ]] && CMD+=("-minsize" "${MIN}")
    [[ "$AG" == "true" ]] && CMD+=("-tabbedout" "${OUTDIR}/unoise3_processing_report.txt")
    [[ "$ALPHA" == "true" ]] && CMD+=("-unoise_alpha" "${ALPHA}")
    echo " - -- --- ---- ---- --- -- -"
    echo " - RUNNING: ${CMD[@]}"
    "${CMD[@]}"
    sed -i 's/>Zotu/>Asv/g' "${OUTDIR}/asvs/asvs.fa"
    typ="asv"
else 
    mkdir -vp "${OUTDIR}/otus"
    CMD=("usearch" "-cluster_smallmem" "${OUTDIR}/uniques.fa" "-quiet" "-id" "${UP}" "-sortedby" "size" "-relabel" "Otu" "-centroids" "${OUTDIR}/otus/otus.fa")
    [[ "$MIN" == "true" ]] && CMD+=("-minsize" "${MIN}")
    [[ "$AG" == "true" ]] && CMD+=("-uc" "${OUTDIR}/uparse_processing_report.txt")
    echo " - -- --- ---- ---- --- -- -"
    echo " - RUNNING: ${CMD[@]}"
    "${CMD[@]}"
    typ="otu"
fi

#create an ASV table ("Input should be reads before quality filtering and before discarding low-abundance unique sequences, e.g. singletons")
CMD=("usearch" "--otutab" "${OUTDIR}/combined.fq" "-quiet" "-zotus" "${OUTDIR}/${typ}s/${typ}s.fa" "-otutabout" "${OUTDIR}/${typ}s/${typ}_table_00.txt")
[[ "$TBLID" == "true" ]] && CMD+=("-id" "${TBLID}")
[[ "$UN" == "true" ]] && CMD+=("-notmatchedfq" "${OUTDIR}/unbinned_reads.fq")
[[ "$RMF" == "true" ]] && CMD+=("-mapout" "${OUTDIR}/read_binning_report.txt")
echo " "
echo " - -- --- ---- ---- --- -- -"
echo " - RUNNING: ${CMD[@]}"
"${CMD[@]}"

if [[ "$ALGORITHM" == "UNOISE3" ]]; then
    sed -i 's/#OTU/#ASV/g' "${OUTDIR}/${typ}s/${typ}_table_00.txt"
fi

# Run python script to sort qiime table
qiime_table_sorter.py "${OUTDIR}/${typ}s/${typ}_table_00.txt" "${OUTDIR}/${typ}s/${typ}_table_01.txt"

#source for bash helper functions
source "qiime_shell_helper_functions.sh"

#convert to biom
OTBL="${typ}_table_01"
txt2biom_notax "${OUTDIR}/${typ}s/${OTBL}.txt" "${OUTDIR}/${typ}s/${OTBL}.biom"

./add_counts_to_fasta.py "${OUTDIR}/${typ}s/${typ}_table_01.txt" "${OUTDIR}/${typ}s/${typ}s.fa" "${OUTDIR}/${typ}s/${typ}s_counts.fa"

mkdir -vp "${OUTDIR}/${typ}s/blast"
mv -v "${OUTDIR}/${typ}s/${typ}s_counts.fa" "${OUTDIR}/${typ}s/blast"

# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# OPTIONAL ALTERNATVIE OTU/ASV CALLING STRATEGY 2
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
if [[ "$MAPF" != "false" ]]; then
    echo
    echo " - -- --- ---- ---- --- -- -"
    echo " - Beginning Alternative Strategy 2"
    echo
    mkdir -vp "${OUTDIR}/${typ}s/STRATEGY2"
    mkdir -vp "${OUTDIR}/${typ}s/STRATEGY2/subset_filtered_fastas"
    
    # Create an associative array to store unique values from $COL and their associated sample IDs
    declare -A value_to_sampleIDs
    
    # Extract the index of the column headers for $COL and #SampleID
    col_index=$(head -1 "$MAPF" | tr '\t' '\n' | grep -n "^$COL$" | cut -d: -f1)
    sampleID_index=$(head -1 "$MAPF" | tr '\t' '\n' | grep -n "^#SampleID$" | cut -d: -f1)
    
    # Ensure the column was found
    if [[ -z "$col_index" || -z "$sampleID_index" ]]; then
    echo "Error: Column $COL or #SampleID not found in $MAPF"
    exit 1
    fi
    
    # Build the associative array with unique values from $COL and associated sample IDs
    awk -v col="$col_index" -v sid="$sampleID_index" 'NR > 1 { value[$col] = value[$col] ? value[$col] "," $sid : $sid } END { for (v in value) print v "\t" value[v] }' FS='\t' "$MAPF" > tmp_value_to_sampleIDs.txt
    
    # Iterate through each unique value in the tmp file
    while IFS=$'\t' read -r value sampleIDs; do
    # Prepare a list of sample IDs for the current value
    IFS=',' read -r -a sampleID_array <<< "$sampleIDs"
    output_fasta="${OUTDIR}/${typ}s/STRATEGY2/subset_filtered_fastas/${value}_subset.fa"
      
    echo "Subsetting sequences for value: $value -> ${output_fasta}"
      
    # Loop through each sample ID and extract matching sequences
    for sampleID in "${sampleID_array[@]}"; do
        # Use awk to match the sample ID in the header and collect the entire sequence
        awk -v sampleID=";sample=${sampleID}" '
        BEGIN {in_sample=0}
        /^>/ { 
            if (index($0, sampleID) > 0) {
            in_sample=1; 
            print $0;  
            } else {
            in_sample=0;  
            }
        }
        in_sample && !/^>/ { 
            print $0
        }
        ' "${OUTDIR}/filtered.fa" # Process the FASTA file
    done | seqkit seq - > "$output_fasta" # Use seqkit to format the sequences
    done < tmp_value_to_sampleIDs.txt
    
    
    # Clean up temporary file
    rm -rf tmp_value_to_sampleIDs.txt
    
    #find unique sequences
    for F in "${OUTDIR}/${typ}s/STRATEGY2/subset_filtered_fastas/"*subset.fa; do
        BASE=$(echo "$F" | sed 's/\.fa$//')
        CMD=("usearch" "-fastx_uniques" "${F}" "-quiet" "-fastaout" "${BASE}.uniques.fa" "-sizeout" "-relabel" "Uniq")
        [[ "$UG" == "true" ]] && CMD+=("-tabbedout" "${OUTDIR}/uniques_binning_report.txt")
        [[ "$US" == "true" ]] && CMD+=("-minuniquesize" "${US}")
        echo
        echo " - -- --- ---- ---- --- -- -"
        echo " - RUNNING: ${CMD[@]}"
        "${CMD[@]}"

        # sort by size in preparation for algorithm
        echo
        echo " - -- --- ---- ---- --- -- -"
        echo " - RUNNING: usearch -sortbysize ${BASE}.uniques.fa -fastaout ${BASE}.temp.fa -minsize 1"
        usearch -sortbysize "${BASE}.uniques.fa" -quiet -fastaout "${BASE}.temp.fa" -minsize 1
        mv "${BASE}.temp.fa" "${BASE}.uniques.fa"
    
        # create otus/asvs using algorithm of choice
        if [[ "$ALGORITHM" == "UNOISE3" ]]; then
            CMD=("usearch" "-unoise3" "${BASE}.uniques.fa" "-quiet" "-zotus" "${BASE}.asvs.fa")
            [[ "$MIN" == "true" ]] && CMD+=("-minsize" "${MIN}")
            [[ "$AG" == "true" ]] && CMD+=("-tabbedout" "${BASE}.unoise3_processing_report.txt")
            [[ "$ALPHA" == "true" ]] && CMD+=("-unoise_alpha" "${ALPHA}")
            echo
            echo " - -- --- ---- ---- --- -- -"
            echo " - RUNNING: ${CMD[@]}"
            "${CMD[@]}"
            sed -i 's/>Zotu/>Asv/g' "${BASE}.asvs.fa"
            typ="asv"
        else 
            CMD=("usearch" "-cluster_smallmem" "${BASE}.uniques.fa" "-quiet" "-id" "${UP}" "-sortedby" "size" "-relabel" "Otu" "-centroids" "${BASE}.otus.fa")
            [[ "$MIN" == "true" ]] && CMD+=("-minsize" "${MIN}")
            [[ "$AG" == "true" ]] && CMD+=("-uc" "${BASE}.uparse_processing_report.txt")
            echo
            echo " - -- --- ---- ---- --- -- -"
            echo " - RUNNING: ${CMD[@]}"
            "${CMD[@]}"
            typ="otu"
        fi
    done
    
    # MERGE OTUs/ASVs
    for F in "${OUTDIR}/${typ}s/STRATEGY2/subset_filtered_fastas/"*subset.${typ}s.fa; do
        id=$(basename "$F" | sed "s/.subset.${typ}s.fa$//")
        awk -v suffix="_${COL}_${id}" 'BEGIN {OFS = ""} /^>/ {print $1, suffix; next} {print}' "$F" > "${F}.mod" 
    done
    mkdir -vp "${OUTDIR}/${typ}s/STRATEGY2/${typ}s"
    cat "${OUTDIR}/${typ}s/STRATEGY2/subset_filtered_fastas/"*mod > "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}s.fa"
    merge_fastas.py "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}s.fa"
    mv "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}s_unique.fa" "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}s.fa"
    
    #create a table ("Input should be reads before quality filtering and before discarding low-abundance unique sequences, e.g. singletons")
    CMD=("usearch" "--otutab" "${OUTDIR}/combined.fq" "-quiet" "-zotus" "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}s.fa" "-otutabout" "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}_table_00.txt")
    [[ "$TBLID" == "true" ]] && CMD+=("-id" "${TBLID}")
    echo
    echo " - -- --- ---- ---- --- -- -"
    echo " - RUNNING: ${CMD[@]}"
    "${CMD[@]}"


    if [[ "$ALGORITHM" == "UNOISE3" ]]; then
        sed -i 's/#OTU/#ASV/g' "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}_table_00.txt"
    fi
    
    # Run python script to sort qiime table by decreasing abundance of taxonomic unit
    qiime_table_sorter.py "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}_table_00.txt" "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}_table_01.txt"
    
    # Add counts to taxonomic units and sort them
    ./add_counts_to_fasta.py --typ "${typ^}" --modify-headers "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}_table_00.txt" "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}s.fa" "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}s_counts.fa"
    mv "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}s_counts.fa" "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}s.fa"
    #create a table AGAIN with correct headers.
    CMD=("usearch" "--otutab" "${OUTDIR}/combined.fq" "-quiet" "-zotus" "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}s.fa" "-otutabout" "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}_table_00.txt")
    [[ "$TBLID" == "true" ]] && CMD+=("-id" "${TBLID}")
    [[ "$UN" == "true" ]] && CMD+=("-notmatchedfq" "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/unbinned_reads.fq")
    [[ "$RMF" == "true" ]] && CMD+=("-mapout" "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/read_binning_report.txt")
    echo
    echo " - -- --- ---- ---- --- -- -"
    echo " - RUNNING: ${CMD[@]}"
    "${CMD[@]}"

    if [[ "$ALGORITHM" == "UNOISE3" ]]; then
        sed -i 's/#OTU/#ASV/g' "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}_table_00.txt"
    fi
    
    # Run python script to sort qiime table by decreasing abundance of taxonomic unit
    qiime_table_sorter.py "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}_table_00.txt" "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}_table_01.txt"
    
    # Add counts to taxonomic units and sort them
    ./add_counts_to_fasta.py "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}_table_01.txt" "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}s.fa" "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}s_counts.fa"

    #source for bash helper functions
    source "qiime_shell_helper_functions.sh"
    
    #convert to biom
    OTBL="${typ}_table_01"
    txt2biom_notax "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${OTBL}.txt" "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${OTBL}.biom"
    
    # move counts taxonomic unit file to the blast folder in preparation for blast
    mkdir -vp "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/blast"
    mv -v "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/${typ}s_counts.fa" "${OUTDIR}/${typ}s/STRATEGY2/${typ}s/blast"
fi

# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# OPTIONAL ALTERNATVIE OTU/ASV CALLING METHOD 3
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
if [[ "$PRE" != "false" ]]; then
    echo
    echo " - -- --- ---- ---- --- -- -"
    echo "Beginning Alternative Strategy 3"
    mkdir -vp "${OUTDIR}/${typ}s/STRATEGY3"
    mkdir -vp "${OUTDIR}/${typ}s/STRATEGY3/${typ}s"
    
    # make an otu table using the pre-existing OTUs/ASVs and collect remaining unbinned reads.
    CMD=("usearch" "--otutab" "${OUTDIR}/combined.fq" "-quiet" "-zotus" "${PRE}" "-otutabout" "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/pre-existing_${typ}_table_00.txt" "-notmatchedfq" "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/reads_unbinned_in_pre-exisiting.fq")
    [[ "$TBLID" == "true" ]] && CMD+=("-id" "${TBLID}")
    [[ "$RMF" == "true" ]] && CMD+=("-mapout" "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/pre-existing_read_binning_report.txt")
    echo
    echo " - -- --- ---- ---- --- -- -"
    echo " - RUNNING: ${CMD[@]}"
    "${CMD[@]}"
    
    # store unbinned reads as a variable to feed into the OTU/ASV generation pipeline
    F="${OUTDIR}/${typ}s/STRATEGY3/${typ}s/reads_unbinned_in_pre-exisiting.fq"
    FQ_base="${F%.fq}"
    
    # Filter the unbinned reads
    CMD=("usearch" "-threads" "1" "-fastq_filter" "${F}" "-quiet" "-fastq_maxee" "${FM}" "-fastaout" "${FQ_base}.filtered.fa")
    [[ "$FR" == "true" ]] && CMD+=("-fastqout_discarded" "${FQ_base}_discarded.fq")
    [[ "$FE" == "true" ]] && CMD+=("-fastq_eeout")
    [[ "$FD" == "true" ]] && CMD+=("-fastq_maxns" "${FD}")
    echo
    echo " - -- --- ---- ---- --- -- -"
    echo " - RUNNING: ${CMD[@]}"
    "${CMD[@]}"
    
    #find unique sequences
    CMD=("usearch" "-fastx_uniques" "${FQ_base}.filtered.fa" "-quiet" "-fastaout" "${FQ_base}.uniques.fa" "-sizeout" "-relabel" "Uniq")
    [[ "$UG" == "true" ]] && CMD+=("-tabbedout" "${FQ_base}.uniques_binning_report.txt")
    [[ "$US" == "true" ]] && CMD+=("-minuniquesize" "$US")
    echo
    echo " - -- --- ---- ---- --- -- -"
    echo " - RUNNING: ${CMD[@]}"
    "${CMD[@]}"

    # sort by size in preparation for algorithm
    echo
    echo " - -- --- ---- ---- --- -- -"
    echo " - RUNNING: usearch -sortbysize "${FQ_base}.uniques.fa" "-quiet" -fastaout "${OUTDIR}/temp.fa" -minsize 1"
    usearch -sortbysize "${FQ_base}.uniques.fa" -quiet -fastaout "${OUTDIR}/temp.fa" -minsize 1
    mv "${OUTDIR}/temp.fa" "${FQ_base}.uniques.fa"
    
    # create otus/asvs using algorithm of choice
    if [[ "$ALGORITHM" == "UNOISE3" ]]; then
        CMD=("usearch" "-unoise3" "${FQ_base}.uniques.fa" "-quiet" "-zotus" "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/unbinned_${typ}s.fa")
        [[ "$MIN" == "true" ]] && CMD+=("-minsize" "$MIN")
        [[ "$AG" == "true" ]] && CMD+=("-tabbedout" "${OUTDIR}/${typ}s/STRATEGY3/unbinned_unoise3_processing_report.txt")
        [[ "$ALPHA" == "true" ]] && CMD+=("-unoise_alpha" "$ALPHA")
        echo
        echo " - -- --- ---- ---- --- -- -"
        echo " - RUNNING: ${CMD[@]}"
        "${CMD[@]}"

        sed -i 's/>Zotu/>Asv/g' "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/unbinned_${typ}s.fa"
        typ="asv"
    else 
        CMD=("usearch" "-cluster_smallmem" "${FQ_base}.uniques.fa" "-quiet" "-id" "$UP" "-sortedby" "size" "-relabel" "Otu" "-centroids" "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/unbinned_${typ}s.fa")
        [[ "$MIN" == "true" ]] && CMD+=("-minsize" "$MIN")
        [[ "$AG" == "true" ]] && CMD+=("-uc" "${OUTDIR}/${typ}s/STRATEGY3/unbinned_uparse_processing_report.txt")
        echo
        echo " - -- --- ---- ---- --- -- -"
        echo " - RUNNING: ${CMD[@]}"
        "${CMD[@]}"
        typ="otu"
    fi

    if [ -s "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/unbinned_${typ}s.fa" ]; then
        # rename OTUs so we know they are from the "unbinned" OTUs
        sed -i 's/>/>unbinned_/g' "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/unbinned_${typ}s.fa"
    
        #create an ASV table ("Input should be reads before quality filtering and before discarding low-abundance unique sequences, e.g. singletons")
        CMD=("usearch" "--otutab" "$F" "-quiet" "-zotus" "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/unbinned_${typ}s.fa" "-otutabout" "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/unbinned_${typ}_table_00.txt")
        [[ "$TBLID" == "true" ]] && CMD+=("-id" "$TBLID")
        [[ "$UN" == "true" ]] && CMD+=("-notmatchedfq" "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/unbinned_unbinned_reads.fq")
        [[ "$RMF" == "true" ]] && CMD+=("-mapout" "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/unbinned_read_binning_report.txt")
        echo
        echo " - -- --- ---- ---- --- -- -"
        echo " - RUNNING: ${CMD[@]}"
        "${CMD[@]}"
    
        # Merge tables - put "zero" in for samples that aren't there. 
        merge_tables.py "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/pre-existing_${typ}_table_00.txt" "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/unbinned_${typ}_table_00.txt" "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/${typ}_table_00.txt"
        rm -rf 
        cat "${PRE}" "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/unbinned_${typ}s.fa" > "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/${typ}s.fa"
    else
        # File is empty or doesn't exist, exit the loop or script
        echo "No reads found in ${OUTDIR}/${typ}s/STRATEGY3/${typ}s/unbinned_${typ}s.fa, meaning that no ${typ}s were generated from the unbinned reads. This may be due to you setting a large minsize parameter or other similar issue."
        cp "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/pre-existing_${typ}_table_00.txt" "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/${typ}_table_00.txt" 
        cp "${PRE}" "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/${typ}s.fa"
    fi

    if [[ "$ALGORITHM" == "UNOISE3" ]]; then
        sed -i 's/#OTU/#ASV/g' "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/${typ}_table_00.txt"
    fi
    
    # Run python script to sort qiime table by decreasing abundance of taxonomic unit
    qiime_table_sorter.py "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/${typ}_table_00.txt" "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/${typ}_table_01.txt"
        
    # Add counts to taxonomic units and sort them
    ./add_counts_to_fasta.py "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/${typ}_table_00.txt" "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/${typ}s.fa" "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/${typ}s_counts.fa"

    #source for bash helper functions
    source "qiime_shell_helper_functions.sh"
        
    #convert to biom
    txt2biom_notax "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/${typ}_table_01.txt" "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/${typ}_table_01.biom"
    
    # move counts taxonomic unit file to the blast folder in preparation for blast
    mkdir -vp "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/blast"
    mv -v "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/${typ}s_counts.fa" "${OUTDIR}/${typ}s/STRATEGY3/${typ}s/blast"
fi

echo " - -- --- ---- ---- --- -- -
Final Recommendations
 - -- --- ---- ---- --- -- -
It's time to assign taxonomy! You will either do this by BLAST-ing your sequences against a 
local BLAST database and optionally by using Qiime2 to classify them.
Make sure you have an up-to-date NCBI nucleotide database (or other database if preferred) 
available and bound to your singularity container (see tutorial) as well as sufficient 
computational power.
" | tee /dev/tty