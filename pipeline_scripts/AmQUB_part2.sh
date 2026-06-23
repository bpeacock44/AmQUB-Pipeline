#!/bin/bash
# See mbio_tutorial.md for further guidance!

### USAGE ###
# This script accepts input in two ways:
# 1. Direct command-line arguments:
#    AmQUB_part2.sh -f ID1_raw_output -p ID1_map.txt -m 2 -r 200-301 -i 10
#    AmQUB_part2.sh -f ID2_raw_output -p ID2_map.subset.txt

## Required Flags
# -f Part 1 output folder you want to further process
# -p Mapping file (this may be a new subset of the previous mapping file if you want to exclude samples)

## Optional Flags:
# -s Subset ID. This will be appended to the name of the output folder so you can remember which map you used. Helpful if you are going to make different subsets from the same dataset. (Default is none, so output folder will have no subset ID attached.)
# -m The number of allowed mismatches. This should only be used if it was first used in part 1.
# -r Trim length stats range (default is from 1 to the end of the read; e.g. 1-301 for a 301 base dataset)
# -i Trim length stats interval (default is every 25 and the highest 5)

# 2. A parameter template file:
#    AmQUB_part2.sh params.csv
#    Where params.csv contains the following rows, comma delimited, no white space between.
#    The labels at the beginning of each row should be the same as below.
#        Part 1 Output Folder,ID1_raw_output
#        Mapping File,ID1_map.txt
#        Subset ID,sub2
#        Mismatch Bases,2
#        Trim Length Stats Range,200-301
#        Trim Length Stats Interval,10
#
# Only "Part 1 Output Folder" and "Mapping File" are required. Any optional line may be
# omitted (or set to DEFAULT) to use its default value.

set -e  # Exit on error

# Custom error handler function
error_handler() {
    local error_message=$1
    echo "Error: Issue encountered with '$error_message'"
}

# Trap errors and call the error handler
trap 'error_handler "$BASH_COMMAND"' ERR

# Single-flowcell inputs
OUTF=""
MAPF=""
SUB="DEFAULT"
mmatchnum="0"
RANGE="DEFAULT"
INV="DEFAULT"

# Function to check if a required label exists in the file
check_label_present() {
    local label="$1"
    if ! grep -iq "^$label," "$2"; then
        echo "Error: Missing expected label '$label' in parameter file '$2'."
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
    check_label_present "Part 1 Output Folder" "$file"
    check_label_present "Mapping File" "$file"
    # The remaining labels are optional and default if absent.

    while IFS= read -r line; do
        if [[ "$line" == "Part 1 Output Folder,"* ]]; then
            value="${line#Part 1 Output Folder,}"; OUTF="$(trim "${value%%,*}")"
        elif [[ "$line" == "Mapping File,"* ]]; then
            value="${line#Mapping File,}"; MAPF="$(trim "${value%%,*}")"
        elif [[ "$line" == "Subset ID,"* ]]; then
            value="${line#Subset ID,}"; SUB="$(trim "${value%%,*}")"
        elif [[ "$line" == "Mismatch Bases,"* ]]; then
            value="${line#Mismatch Bases,}"; mmatchnum="$(trim "${value%%,*}")"
        elif [[ "$line" == "Trim Length Stats Range,"* ]]; then
            value="${line#Trim Length Stats Range,}"; RANGE="$(trim "${value%%,*}")"
        elif [[ "$line" == "Trim Length Stats Interval,"* ]]; then
            value="${line#Trim Length Stats Interval,}"; INV="$(trim "${value%%,*}")"
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

Part 2: Demultiplexing

 - -- --- ---- ---- --- -- -"

# Check if the first argument is a parameter file
if [[ -f "$1" ]]; then
    echo "Reading parameters from file: $1
 - -- --- ---- ---- --- -- -"
    parse_file_input "$1"
else
    # Process command-line arguments
    while getopts ":f:p:s:m:r:i:" opt; do
        case $opt in
            f) OUTF="$OPTARG" ;;
            p) MAPF="$OPTARG" ;;
            s) SUB="$OPTARG" ;;
            m) mmatchnum="$OPTARG" ;;
            r) RANGE="$OPTARG" ;;
            i) INV="$OPTARG" ;;
            :)  echo "Error: Option -$OPTARG requires an argument." >&2; exit 1 ;;
            \?) echo "Error: Invalid option -$OPTARG" >&2; exit 1 ;;
        esac
    done
fi

# Normalize the mismatch value (handles blank or DEFAULT)
mmatchnum="$(trim "$mmatchnum")"
if [[ -z "$mmatchnum" || "$mmatchnum" == "DEFAULT" ]]; then
    mmatchnum="0"
fi

# Validate the mismatch value (part 2 selects the matching part-1 file by this
# number; it is not re-applied as a flag, so no VAR is needed here)
case "$mmatchnum" in
    0|1|2|3|4|5) ;;
    *) echo "Error: Invalid mismatch value: $mmatchnum"; exit 1 ;;
esac

# Tidy required inputs
OUTF="$(trim "$OUTF")"
OUTF="${OUTF%/}"                 # drop a trailing slash so basename logic is reliable
MAPF="$(trim "$MAPF")"

# Ensure the required inputs are present
if [[ -z "$OUTF" || -z "$MAPF" ]]; then
    echo "Error: No valid input provided."
    echo "Usage: AmQUB_part2.sh -f <part 1 output folder> -p <mapping file> [-s <subset ID> -m <mismatches>] [-r <range>] [-i <interval>]"
    echo "OR use a parameter file: AmQUB_part2.sh params.csv"
    exit 1
fi

# Check that the folder and map file exist
[[ -e "$OUTF" ]]    || { echo "Error: Required folder '$OUTF' not found!"; exit 1; }
[[ -e "$MAPF" ]] || { echo "Error: Required file '$MAPF' not found!"; exit 1; }

# Validate that the mapping file has the columns the pipeline needs (fail fast).
maphdr=$(head -n 1 "$MAPF" | tr -d '\r')
for col in "#SampleID" "BarcodeSequence"; do
    printf '%s' "$maphdr" | tr '\t' '\n' | grep -qxF -- "$col" || {
        echo "Error: Mapping file '$MAPF' is missing required column '$col'."
        exit 1
    }
done

# Derive the dataset base name from the folder's BASENAME (so a path or trailing
# slash on -f doesn't break the file search or the output folder name).
BASE="${OUTF##*/}"          # strip any path
BASE="${BASE#part1_}"       # strip the part1_ prefix
BASE="${BASE%_output}"      # strip the _output suffix

# Find the part 1 fastq (the read file, not the _BC barcode file)
matching_files=($(find -L "$OUTF" -type f -name "${BASE}*M${mmatchnum}.fq" ! -name "*_BC.M*"))
if [[ ${#matching_files[@]} -gt 1 ]]; then
    echo "Error: More than one matching file found."
    exit 1
fi

FQ="${matching_files[0]}"
if [ -z "$FQ" ]; then
    echo "No matching file found for ${BASE}*M${mmatchnum}.fq. Did you run part 1?"
    exit 1
fi

BC="${FQ/.M${mmatchnum}.fq/_BC.M${mmatchnum}.fq}"
if [ ! -e "$BC" ]; then
    echo "No barcode file found at ${BC}. Did you run part 1?"
    exit 1
fi

# Calculate the maximum length of the first 100 reads
MAX_LENGTH=$(head -n 400 "${FQ}" | awk '{if(NR%4==2) print length($1)}' | sort -nr | head -n 1)

# Resolve the range (start-end) and interval, applying defaults
RANGE="$(trim "$RANGE")"
if [[ -z "$RANGE" || "$RANGE" == "DEFAULT" ]]; then
    BEG="1"
    FIN=$((MAX_LENGTH - 5))
else
    IFS='-' read -ra SPLIT_RANGE <<< "$RANGE"
    if [[ ${#SPLIT_RANGE[@]} -ne 2 ]]; then
        echo "Error: Invalid range format '$RANGE'. Expected format: start-end (e.g., 200-301)"
        exit 1
    fi
    BEG="$(trim "${SPLIT_RANGE[0]}")"
    FIN="$(trim "${SPLIT_RANGE[1]}")"
fi

INV="$(trim "$INV")"
if [[ -z "$INV" || "$INV" == "DEFAULT" ]]; then
    INV="25"
fi

# Step 1: Validate that BEG, INV, and FIN are valid integers
if ! [[ "$BEG" =~ ^[0-9]+$ ]] || ! [[ "$INV" =~ ^[0-9]+$ ]] || ! [[ "$FIN" =~ ^[0-9]+$ ]]; then
    echo "Error: Range and interval must be made up of valid integers."
    exit 1
fi

# Step 2: Check if BEG is smaller than FIN
if [ "$BEG" -ge "$FIN" ]; then
    echo "Error: The beginning of the range must be smaller than the end of the range."
    exit 1
fi

# Step 3: Check if INV is smaller than the difference between FIN and BEG
RANGE_DIFF=$((FIN - BEG))
if [ "$INV" -ge "$RANGE_DIFF" ]; then
    echo "Error: The interval is too large. It must be smaller than the range."
    exit 1
fi

# Step 4: Ensure MAX_LENGTH is greater than or equal to FIN
if [ "$MAX_LENGTH" -lt "$FIN" ]; then
    echo "Error: The upper limit of trim length stats range is higher than the maximum read length ($MAX_LENGTH)."
    exit 1
fi

######### BEGIN PROCESSING

# Build the output directory name from BASE (path-safe)
if [[ "$SUB" == "DEFAULT" ]]; then
    ODIR="part2_${BASE}_output"
else
    ODIR="part2_${BASE}_${SUB}_subset_output"
fi

# make the output directory
if ! mkdir -p "${ODIR}"; then
    echo "Error: Failed to create directory '$ODIR'."
    exit 1
fi

# initiate log (tee to terminal and log file)
timestamp="$(date +"%y%m%d_%H:%M")"
output_file="${ODIR}/part2_${timestamp}.log"
exec > >(tee "$output_file")
exec 2> >(tee -a "$output_file" >&2)

# create header for log file
echo "Processing ${OUTF} with mapping file ${MAPF}, allowing ${mmatchnum} mismatches.
Range and Interval of Stats: ${BEG}-${MAX_LENGTH}, every ${INV} bases.
 - -- --- ---- ---- --- -- -"

# Create barcodes.fa file from the mapping file
grep -P "^[A-Z]" "${MAPF}" | awk '{print ">"$1"\n"$2}' > "${ODIR}/barcodes.fa"
[[ ! -e "${ODIR}/barcodes.fa" ]] && echo "File ${ODIR}/barcodes.fa was not generated!" && exit 1

# Demultiplexing
if ! usearch -fastx_demux "${FQ}" -quiet -index "${BC}" -barcodes "${ODIR}/barcodes.fa" -fastqout "${ODIR}/${BASE}.M${mmatchnum}.demux.fq"; then
    echo "Error: usearch failed to run."
    exit 1
fi
cp "${BC}" "${ODIR}/${BASE}_BC.M${mmatchnum}.fq"

# Stats
usearch -fastq_eestats2 "${ODIR}/${BASE}.M${mmatchnum}.demux.fq" -quiet -output "${ODIR}/${BASE}.M${mmatchnum}.eestats.txt" -length_cutoffs "${BEG},${FIN},${INV}"
# Generate the final 5 stats
STARTAT=$((MAX_LENGTH - 5))
INC=1; # Increment value
usearch -fastq_eestats2 "${ODIR}/${BASE}.M${mmatchnum}.demux.fq" -quiet -output "${ODIR}/${BASE}.M${mmatchnum}.eestats.temp.txt" -length_cutoffs ${STARTAT},*,${INC}
# Merge the two
tail -n +7 "${ODIR}/${BASE}.M${mmatchnum}.eestats.temp.txt" >> "${ODIR}/${BASE}.M${mmatchnum}.eestats.txt"
rm -rf "${ODIR}/${BASE}.M${mmatchnum}.eestats.temp.txt"

echo " - -- --- ---- ---- --- -- -
Final Recommendations
 - -- --- ---- ---- --- -- -
Stats ready. Please view this file and select your trim length accordingly:
${ODIR}/${BASE}.M${mmatchnum}.eestats.txt
 - -- --- ---- ---- --- -- -
"