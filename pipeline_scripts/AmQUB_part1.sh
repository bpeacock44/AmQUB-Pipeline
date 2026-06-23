#!/bin/bash
# See mbio_tutorial.md for further guidance.

### USAGE ###
# This script accepts input in two ways:
# 1. Direct command-line arguments:
#    AmQUB_part1.sh -f ID1_raw.fq -p ID1_map.txt -m 2
#    AmQUB_part1.sh -f ID1_raw.fq -p ID1_map.txt

## Required Flags
# -f Fastq file for your flowcell
# -p Mapping file

## Optional Flags:
# -m The number of allowed mismatches. This can be 1-5. Default is 0.

# 2. A parameter template file:
#    AmQUB_part1.sh params.csv
#    Where params.csv contains the following rows, comma delimited.
#    The labels at the beginning of each row should be the same as below.
#        Raw Fastq File,ID1_raw.fq
#        Mapping File,ID1_map.txt
#        Mismatch Bases,2
#
# The "Mismatch Bases" line is optional. Omit it entirely, leave it blank, or set
# it to DEFAULT to use the default of 0.
#
# This script uses the following helper scripts:
# check_barcode_collisions.pl
# filter_barcode_noncollisions.py
# fastq_convert_mm2pm_barcodes.py
# bc_counter.py
# extract_barcodes.go

set -e  # Exit on error

# Single-flowcell inputs
FQ=""
MAPF=""
mmatchnum="0"

# Function to check if a label exists in the parameter file
check_label_present() {
    local label="$1"
    if ! grep -q "$label" "$2"; then
        echo "Error: Missing expected label '$label' in parameter file."
        exit 1
    fi
}

# Strip leading/trailing whitespace from a single value
trim() {
    printf '%s' "$1" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//'
}

# Function to parse the parameter file (one flowcell)
parse_file_input() {
    local file="$1" line value

    # Required labels
    check_label_present "Raw Fastq File," "$file"
    check_label_present "Mapping File," "$file"
    # "Mismatch Bases," is optional and is therefore not required here.

    while IFS= read -r line; do
        if [[ "$line" == "Raw Fastq File,"* ]]; then
            value="${line#Raw Fastq File,}"
            FQ="$(trim "${value%%,*}")"        # first field after the label
        elif [[ "$line" == "Mapping File,"* ]]; then
            value="${line#Mapping File,}"
            MAPF="$(trim "${value%%,*}")"
        elif [[ "$line" == "Mismatch Bases,"* ]]; then
            value="${line#Mismatch Bases,}"
            mmatchnum="$(trim "${value%%,*}")"
        fi
    done < "$file"
}


echo "
        ┌── ===
 ┌──────┤
 │      └── ooo
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
else
    # Process command-line arguments
    while getopts ":f:p:m:" opt; do
        case $opt in
            f) FQ="$OPTARG" ;;
            p) MAPF="$OPTARG" ;;
            m) mmatchnum="$OPTARG" ;;
            :)  echo "Error: Option -$OPTARG requires an argument." >&2; exit 1 ;;
            \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
        esac
    done
fi

# Normalize the mismatch default (handles blank or DEFAULT)
if [[ -z "$mmatchnum" || "$mmatchnum" == "DEFAULT" ]]; then
    mmatchnum="0"
fi

# Ensure the required inputs are present
if [[ -z "$FQ" || -z "$MAPF" ]]; then
    echo "Error: No valid input provided."
    echo "Usage: AmQUB_part1.sh -f <raw fastq file> -p <mapping file> [-m <mismatches>]"
    echo "OR use a parameter file: AmQUB_part1.sh params.csv"
    exit 1
fi

# Check that the input files exist
[[ -e "$FQ" ]]   || { echo "Error: Required file '$FQ' not found!"; exit 1; }
[[ -e "$MAPF" ]] || { echo "Error: Required file '$MAPF' not found!"; exit 1; }

# Validate that the mapping file has the columns the pipeline needs (fail fast,
# before any tool runs, since step 1 and the read counter both rely on these).
maphdr=$(head -n 1 "$MAPF" | tr -d '\r')
for col in "#SampleID" "BarcodeSequence"; do
    printf '%s' "$maphdr" | tr '\t' '\n' | grep -qxF -- "$col" || {
        echo "Error: Mapping file '$MAPF' is missing required column '$col'."
        exit 1
    }
done

# Validate the mismatch value and set the corresponding option
case "$mmatchnum" in
    0) VAR="-m0" ;;
    1) VAR="-m1" ;;
    2) VAR="-m12" ;;
    3) VAR="-m123" ;;
    4) VAR="-m1234" ;;
    5) VAR="-m12345" ;;
    *) echo "Error: Invalid mismatch value: $mmatchnum"; exit 1 ;;
esac

# Prepare output directory
BASE=$(basename "$FQ" .fastq)
BASE=$(basename "$BASE" .fq)
NDIR=$(dirname "$FQ")
OUTDIR="part1_${BASE}_output"
mkdir -p "$OUTDIR"

# Log initialization. Output goes to the terminal AND the log file.
timestamp="$(date +"%y%m%d_%H:%M")"
output_file="${OUTDIR}/part1_${timestamp}.log"
exec > >(tee "$output_file")
exec 2> >(tee -a "$output_file" >&2)

echo "Processing ${FQ} with mapping file ${MAPF}, allowing ${mmatchnum} mismatches
 - -- --- ---- ---- --- -- -"

# Run main pipeline commands
if [ "$mmatchnum" -ne 0 ]; then
    echo "Processing barcode mismatches..."

    # '|| true' keeps a zero-match grep (exit status 1) from tripping 'set -e'
    _BC_=$(grep -cP "^[A-Z]" "${MAPF}" || true)
    check_barcode_collisions.pl -i "${FQ}" -m "${MAPF}" -M${mmatchnum} -C -o "${OUTDIR}/${BASE}.BC${_BC_}_M${mmatchnum}.collisions.txt"
    filter_barcode_noncollisions.py -k -i "${OUTDIR}/${BASE}.BC${_BC_}_M${mmatchnum}.collisions.txt" $VAR --output_for_fastq_convert > "${OUTDIR}/${BASE}_M${mmatchnum}.fbncs"
    fastq_convert_mm2pm_barcodes.py -t read -i "${FQ}" -m "${OUTDIR}/${BASE}_M${mmatchnum}.fbncs" -o "${OUTDIR}/${BASE}.M${mmatchnum}.fq"
    extract_barcodes.go -f "${OUTDIR}/${BASE}.M${mmatchnum}.fq" && mv -v "${OUTDIR}/barcodes.fastq" "${OUTDIR}/${BASE}_BC.M${mmatchnum}.fq"

    # Check if files were generated
    [[ -e "${OUTDIR}/${BASE}.M${mmatchnum}.fq" ]]    || { echo "Error: File ${OUTDIR}/${BASE}.M${mmatchnum}.fq was not generated!"; exit 1; }
    [[ -e "${OUTDIR}/${BASE}_BC.M${mmatchnum}.fq" ]] || { echo "Error: File ${OUTDIR}/${BASE}_BC.M${mmatchnum}.fq was not generated!"; exit 1; }

else
    extract_barcodes.go -f "${FQ}" && mv -v "${NDIR}/barcodes.fastq" "${OUTDIR}/${BASE}_BC.M${mmatchnum}.fq"
    ln -sf "$(realpath "$FQ")" "${OUTDIR}/${BASE}.M${mmatchnum}.fq"
fi

echo "Generating fastq info file..."
usearch -fastx_info "${OUTDIR}/${BASE}.M${mmatchnum}.fq" -quiet -secs 300 -output "${OUTDIR}/${BASE}.M${mmatchnum}.fastq_info.txt"

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
"