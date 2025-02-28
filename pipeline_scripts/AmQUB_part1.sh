#!/bin/bash
# See mbio_tutorial.md for further guidance.

### USAGE ###
# This script accepts input in two ways:
# 1. Direct command-line arguments:
#    AmQUB_part1.sh -f ID1_raw.fq -p ID1_map.txt -m 2
#    AmQUB_part1.sh -f ID2_raw.fq -p ID2_map.txt

## Required Flags
# -f Fastq file for your flowcell
# -p Mapping file 

## Optional Flags:
# -m The number of allowed mismatches. This can be 1-5. Default is 0.

# 2. A parameter template file:
#    AmQUB_part1.sh params.csv
#    Where params.csv contains the following 3 rows, comma delimited.
#    The labels at the beginning of each row should be the same as below.
#        Raw Fastq File,ID1_raw.fq,ID2_raw.fq
#        Mapping File,ID1_map.txt,ID2_map.txt
#        Mismatch Bases,2,DEFAULT (Default is 0)
#
# Note that parameter file can contain multiple flowcells - one per column with 
# associated parameters. If you want to use the default value for all flowcells, you can omit the line.
# If you want to use the default value for only some of them, indicate with "DEFAULT" as above.


set -e  # Exit on error

# Arrays to hold multiple sets of inputs
declare -a FQ_ARRAY MAPF_ARRAY MMATCH_ARRAY

# Function to check if the label exists in the file
check_label_present() {
    local label="$1"
    if ! grep -q "$label" "$2"; then
        echo "Error: Missing expected label '$label' in parameter file." | tee /dev/tty
        exit 1
    fi
}

# Function to parse parameter file and store multiple sets
parse_file_input() {
    # Check if all labels are present in the file
    check_label_present "Raw Fastq File," "$1"
    check_label_present "Mapping File," "$1"
    check_label_present "Mismatch Bases," "$1"

    # Read the file line by line
    while IFS= read -r line; do
        # Check if the line starts with "Raw Fastq File,"
        if [[ "$line" == "Raw Fastq File,"* ]]; then
            values="${line#Raw Fastq File,}"
            IFS=',' read -r -a FQ_ARRAY <<< "$values"  # Populate global array for FQ_ARRAY

        # Check if the line starts with "Mapping File,"
        elif [[ "$line" == "Mapping File,"* ]]; then
            values="${line#Mapping File,}"
            IFS=',' read -r -a MAPF_ARRAY <<< "$values"  # Populate global array for MAPF_ARRAY

        # Check if the line starts with "Mismatch Bases,"
        elif [[ "$line" == "Mismatch Bases,"* ]]; then
            values="${line#Mismatch Bases,}"
            IFS=',' read -r -a MMATCH_ARRAY <<< "$values"  # Populate global array for MMATCH_ARRAY

        fi
    done < "$1"

    # Strip leading/trailing whitespace from all arrays
    for arr in FQ_ARRAY MAPF_ARRAY MMATCH_ARRAY; do
        for i in "${!arr[@]}"; do
            arr[$i]=$(echo "${arr[$i]}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        done
    done
}

# Function to check if all arrays have the same length
check_arrays_length() {
    local fq_len=${#FQ_ARRAY[@]}
    local map_len=${#MAPF_ARRAY[@]}
    local mmatch_len=${#MMATCH_ARRAY[@]}

    if [[ "$fq_len" -ne "$map_len" || "$fq_len" -ne "$mmatch_len" ]]; then
        echo "Error: The arrays have different lengths. Please ensure all arrays are of equal length." | tee /dev/tty
        exit 1
    fi
}


echo "
        ┌── ~<>
 ┌──────┤
 │      └── ~~<>
─┤
 │ ┌── [A m Q U B]  
 └─┤
   └──── ~<>~   

Part 1: Pre-Processing

 - -- --- ---- ---- --- -- -" 
# Check if the first argument is a file (parameter template)
if [[ -f "$1" ]]; then
    echo "Reading parameters from file: $1
 - -- --- ---- ---- --- -- -"
    parse_file_input "$1"
    check_arrays_length
else
    # Process single command-line arguments (overrides parameter file if provided)
    MMATCH_ARRAY=("0")  

    while getopts ":f:p:m:" opt; do
        case $opt in
            f) FQ_ARRAY=("$OPTARG") ;;
            p) MAPF_ARRAY=("$OPTARG") ;;
            m) MMATCH_ARRAY=("$OPTARG") ;;
            \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
        esac
    done

    # Check if the arrays have the same length (after command-line input)
    check_arrays_length
fi


# Ensure at least one dataset is provided
if [[ ${#FQ_ARRAY[@]} -eq 0 || ${#MAPF_ARRAY[@]} -eq 0 ]]; then
    echo "Error: No valid input provided."
    echo "Usage: AmQUB_part1.sh -f <raw fastq file> -p <mapping file> [-m <mismatches>]"
    echo "OR use a parameter file: AmQUB_part1.sh params.csv"
    exit 1
fi

# Iterate over each set of parameters
for i in "${!FQ_ARRAY[@]}"; do
    FQ="${FQ_ARRAY[$i]}"
    MAPF="${MAPF_ARRAY[$i]}"

    # Check if the value in MMATCH_ARRAY is "DEFAULT" and default to 0
    mmatchnum="${MMATCH_ARRAY[$i]}"
    if [[ "$mmatchnum" == "DEFAULT" ]]; then
        mmatchnum="0"
    fi

    # Validate required parameters
    if [[ -z "$FQ" || -z "$MAPF" ]]; then
        echo "Error: Missing required parameters for dataset $((i+1))."
        exit 1
    fi

    # Check if files exist
    [[ -e "$FQ" ]] || { echo "Error: Required file '$FQ' not found!"; exit 1; }
    [[ -e "$MAPF" ]] || { echo "Error: Required file '$MAPF' not found!"; exit 1; }

    # Set mismatch option
    case $mmatchnum in
        0) VAR="-m0" ;;
        1) VAR="-m1" ;;
        2) VAR="-m12" ;;
        3) VAR="-m123" ;;
        4) VAR="-m1234" ;;
        5) VAR="-m12345" ;;
        *) echo "Error: Invalid mismatch value: $mmatchnum"; exit 1 ;;
    esac
done

for i in "${!FQ_ARRAY[@]}"; do
    FQ="${FQ_ARRAY[$i]}"
    MAPF="${MAPF_ARRAY[$i]}"
    
    # Check if the value in MMATCH_ARRAY is "DEFAULT" and default to 0
    mmatchnum="${MMATCH_ARRAY[$i]}"
    if [[ "$mmatchnum" == "DEFAULT" ]]; then
        mmatchnum="0"
    fi

    # Prepare output directory
    BASE=$(basename "$FQ" .fastq)
    BASE=$(basename "$BASE" .fq)
    OUTDIR="part1_${BASE}_output"
    mkdir -p "$OUTDIR"

    # Log initialization
    timestamp="$(date +"%y%m%d_%H:%M")"
    output_file="${OUTDIR}/part1_${timestamp}.log"
    exec > "$output_file"
    exec 2> >(tee -a "$output_file" >&2)
    echo "Processing ${FQ} with mapping file ${MAPF}, allowing ${mmatchnum} mismatches
 - -- --- ---- ---- --- -- -" | tee /dev/tty

    # Run main pipeline commands
    usearch -search_phix "${FQ}" -quiet -notmatchedfq "${OUTDIR}/${BASE}.phiX_clean.fq" -alnout "${OUTDIR}/${BASE}.phiX_clean.alnout"

    if [ "$mmatchnum" -ne 0 ]; then
        echo "Processing barcode mismatches..."

        _BC_=$(grep -cP "^[A-Z]" "${MAPF}")
        check_barcode_collisions.pl -i "${OUTDIR}/${BASE}.phiX_clean.fq" -m "${MAPF}" -M${mmatchnum} -C -o "${OUTDIR}/${BASE}.BC${_BC_}_M${mmatchnum}.collisions.txt"
        filter_barcode_noncollisions.py -k -i "${OUTDIR}/${BASE}.BC${_BC_}_M${mmatchnum}.collisions.txt" $VAR --output_for_fastq_convert > "${OUTDIR}/${BASE}_M${mmatchnum}.fbncs"
        fastq_convert_mm2pm_barcodes.py -t read -i "${OUTDIR}/${BASE}.phiX_clean.fq" -m "${OUTDIR}/${BASE}_M${mmatchnum}.fbncs" -o "${OUTDIR}/${BASE}.M${mmatchnum}.fq"
        extract_barcodes.go -f "${OUTDIR}/${BASE}.M${mmatchnum}.fq" && mv -v "${OUTDIR}/barcodes.fastq" "${OUTDIR}/${BASE}_BC.M${mmatchnum}.fq"

        # Check if files were generated
        [[ -e "${OUTDIR}/${BASE}.M${mmatchnum}.fq" ]] || { echo "Error: File ${OUTDIR}/${BASE}.M${mmatchnum}.fq was not generated!"; exit 1; }
        [[ -e "${OUTDIR}/${BASE}_BC.M${mmatchnum}.fq" ]] || { echo "Error: File ${OUTDIR}/${BASE}_BC.M${mmatchnum}.fq was not generated!"; exit 1; }

    else
        extract_barcodes.go -f "${OUTDIR}/${BASE}.phiX_clean.fq" && mv -v "${OUTDIR}/barcodes.fastq" "${OUTDIR}/${BASE}_BC.M${mmatchnum}.fq"
        ln -sf "${BASE}.phiX_clean.fq" "${OUTDIR}/${BASE}.M${mmatchnum}.fq"
    fi

    echo "Generating fastq info file..."
    usearch -fastx_info "${OUTDIR}/${BASE}.M${mmatchnum}.fq" -quiet -output "${OUTDIR}/${BASE}.M${mmatchnum}.fastq_info.txt"

    echo "Counting reads per barcode..."
    bc_counter.py "${MAPF}" "${OUTDIR}/${BASE}_BC.M${mmatchnum}.fq" "${OUTDIR}/${BASE}.M${mmatchnum}.read_counts.txt"

echo " - -- --- ---- ---- --- -- -
Final Recommendations
 - -- --- ---- ---- --- -- -
Output files are stored in this directory:
${OUTDIR}

The number of reads per sample that resulted from this script for ${BASE} can be found 
in this file: ${OUTDIR}/${BASE}.M${mmatchnum}.read_counts.txt
You should check this file, as it may indicate that you should remove or ignore certain 
samples downstream.

There is also a file describing the makeup of your reads here:
${OUTDIR}/${BASE}.M${mmatchnum}.fastq_info.txt
 - -- --- ---- ---- --- -- -
"  | tee /dev/tty
done