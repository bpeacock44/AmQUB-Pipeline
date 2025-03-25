#!/bin/bash
# See mbio_tutorial.md for further guidance!

### USAGE ###
# This script accepts two types of input:
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
#    Where params.csv contains the following 3 rows, comma delimited, no white space between.
#    The labels at the beginning of each row should be the same as below.
#        Part 1 Output Folder,ID1_raw_output,ID2_raw_output
#        Mapping File,ID1_map.txt,ID2_map.sub2.txt 
#        Subset ID,DEFAULT,sub2
#        Mismatch Bases,2,DEFAULT (Default is 0)
#        Trim Length Stats Range,200-301,DEFAULT (Default is full length of read)
#        Trim Length Stats Interval,10,DEFAULT (Default is every 25 and each of the the highest 5)
#
# Note that parameter file can contain multiple flowcells - one per column with 
# associated parameters. If you want to use the default value for all flowcells, you can omit the line.
# If you want to use the default value for only some of them, indicate with "DEFAULT" as above.

set -e  # Exit on error

# Custom error handler function
error_handler() {
    local error_message=$1
    echo "Error: Issue encountered with '$error_message'" | tee /dev/tty
}

# Trap errors and call the error handler
trap 'error_handler "$BASH_COMMAND"' ERR

# Arrays to hold multiple sets of inputs
declare -a OUTF_ARRAY MAPF_ARRAY SUB_ARRAY MMATCH_ARRAY RANGE_ARRAY INTERVAL_ARRAY

# Function to check if a required label exists in the file
check_label_present() {
    local label="$1"
    if ! grep -iq "^$label," "$2"; then
        echo "Error: Missing expected label '$label' in parameter file '$2'." | tee /dev/tty
        exit 1
    fi
}

# Function to parse parameter file and store multiple sets
parse_file_input() {
    # Validate required labels
    for label in "Part 1 Output Folder" "Mapping File" "Subset ID" "Mismatch Bases" "Trim Length Stats Range" "Trim Length Stats Interval"; do
        check_label_present "$label" "$1"
    done

    # Initialize BEG_ARRAY and FIN_ARRAY
    BEG_ARRAY=()
    FIN_ARRAY=()

    # Read the file line by line
    while IFS= read -r line; do
        if [[ "$line" == "Part 1 Output Folder,"* ]]; then
            IFS=',' read -r -a OUTF_ARRAY <<< "${line#Part 1 Output Folder,}"

        elif [[ "$line" == "Mapping File,"* ]]; then
            IFS=',' read -r -a MAPF_ARRAY <<< "${line#Mapping File,}"

        elif [[ "$line" == "Subset ID,"* ]]; then
            IFS=',' read -r -a SUB_ARRAY <<< "${line#Subset ID,}"

        elif [[ "$line" == "Mismatch Bases,"* ]]; then
            IFS=',' read -r -a MMATCH_ARRAY <<< "${line#Mismatch Bases,}"

        elif [[ "$line" == "Trim Length Stats Range,"* ]]; then
            IFS=',' read -r -a RANGE_ARRAY <<< "${line#Trim Length Stats Range,}"
            
            # Check each range entry in RANGE_ARRAY
            for range in "${RANGE_ARRAY[@]}"; do
                # Strip whitespace
                range=$(echo "$range" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
                if [[ "$range" == "DEFAULT" ]]; then
                    # If "DEFAULT", add "DEFAULT" to both BEG_ARRAY and FIN_ARRAY
                    BEG_ARRAY+=("DEFAULT")
                    FIN_ARRAY+=("DEFAULT")
                else
                    # Otherwise, process the range normally
                    IFS='-' read -ra SPLIT_RANGE <<< "$range"
                    if [[ ${#SPLIT_RANGE[@]} -ne 2 ]]; then
                        echo "Error: Invalid range format '$range'. Expected format: start-end (e.g., 200-301)"
                        exit 1
                    fi
                    BEG_ARRAY+=("${SPLIT_RANGE[0]}")
                    FIN_ARRAY+=("${SPLIT_RANGE[1]}")
                fi
            done

        elif [[ "$line" == "Trim Length Stats Interval,"* ]]; then
            IFS=',' read -r -a INTERVAL_ARRAY <<< "${line#Trim Length Stats Interval,}"
        fi
    done < "$1"
}

# Function to check if all arrays have the same length
check_arrays_length() {
    local outf_len=${#OUTF_ARRAY[@]}
    local map_len=${#MAPF_ARRAY[@]}
    local sub_len=${#SUB_ARRAY[@]}
    local mmatch_len=${#MMATCH_ARRAY[@]}
    local beg_len=${#BEG_ARRAY[@]}
    local fin_len=${#FIN_ARRAY[@]}
    local interval_len=${#INTERVAL_ARRAY[@]}

    if [[ "$outf_len" -ne "$map_len" || "$outf_len" -ne "$mmatch_len" || "$outf_len" -ne "$beg_len" || "$outf_len" -ne "$fin_len" || "$outf_len" -ne "$interval_len" || "$outf_len" -ne "$sub_len" ]]; then
        echo "Error: The arrays have different lengths. Please ensure all arrays are of equal length." | tee /dev/tty
        exit 1
    fi
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
    check_arrays_length
else
    # Process single command-line arguments
    OUTF_ARRAY=()  
    MAPF_ARRAY=() 
    SUB_ARRAY=() 
    MMATCH_ARRAY=("0")  
    INTERVAL_ARRAY=()  
    BEG_ARRAY=()
    FIN_ARRAY=()

    while getopts ":f:p:s:m:r:i:" opt; do
        case $opt in
            f) OUTF_ARRAY=("$OPTARG") ;;
            p) MAPF_ARRAY=("$OPTARG") ;;
            s) SUB_ARRAY=("$OPTARG") ;;
            m) MMATCH_ARRAY=("$OPTARG") ;;
            r) 
                if [[ -n "$OPTARG" ]]; then  # Check if the argument for -r is provided
                    IFS='-' read -ra RANGE <<< "$OPTARG"
                    if [[ ${#RANGE[@]} -ne 2 ]]; then
                        echo "Error: Invalid range format '$OPTARG'. Expected format: start-end (e.g., 200-301)"
                        exit 1
                    fi
                    BEG_ARRAY=("${RANGE[0]}")  # First part (start)
                    FIN_ARRAY=("${RANGE[1]}")  # Second part (end)
                else
                    echo "Warning: Range (-r) not provided. Default values will be used."
                    BEG_ARRAY=("DEFAULT")
                    FIN_ARRAY=("DEFAULT")
                fi
                ;;
            i) INTERVAL_ARRAY=("$OPTARG") ;;
            :) 
                echo "Error: Option -$OPTARG requires an argument." >&2
                exit 1
                ;;
            \?) 
                echo "Error: Invalid option -$OPTARG" >&2
                exit 1
                ;;
        esac
    done

    # Check if INTERVAL_ARRAY is empty, if so set to "DEFAULT"
    if [ ${#INTERVAL_ARRAY[@]} -eq 0 ]; then
        INTERVAL_ARRAY=("DEFAULT")
    fi

    if [ ${#SUB_ARRAY[@]} -eq 0 ]; then
        SUB_ARRAY=("DEFAULT")
    fi

    # Check if BEG_ARRAY is empty, if so set to "DEFAULT"
    if [ ${#BEG_ARRAY[@]} -eq 0 ]; then
        BEG_ARRAY=("DEFAULT")
    fi

    # Check if FIN_ARRAY is empty, if so set to "DEFAULT"
    if [ ${#FIN_ARRAY[@]} -eq 0 ]; then
        FIN_ARRAY=("DEFAULT")
    fi
fi


# Ensure at least one dataset is provided
if [[ ${#OUTF_ARRAY[@]} -eq 0 || ${#MAPF_ARRAY[@]} -eq 0 ]]; then
    echo "Error: No valid input provided."
    echo "Usage: AmQUB_part2.sh -f <part 1 output folder> -p <mapping file> [-s <subset ID> -m <mismatches>] [-r <range>] [-i <interval>]"
    echo "OR use a parameter file: AmQUB_part2.sh params.csv"
    exit 1
fi

# Iterate over each set of parameters
for i in "${!OUTF_ARRAY[@]}"; do

    # Set up mismatches if indicated.
    mmatchnum="${MMATCH_ARRAY[$i]}"
    mmatchnum=$(echo "$mmatchnum" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    if [[ "$mmatchnum" == "DEFAULT" ]]; then
        mmatchnum="0"
    fi

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

    OUTF="${OUTF_ARRAY[$i]}"
    OUTF=$(echo "$OUTF" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

    MAPFILE="${MAPF_ARRAY[$i]}"
    MAPFILE=$(echo "$MAPFILE" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

    # Validate required parameters
    if [[ -z "$OUTF" || -z "$MAPFILE" ]]; then
        echo "Error: Missing required parameters for dataset $((i+1))."
        exit 1
    fi

    # Check if the folder and map file indicated exist
    [[ -e "$OUTF" ]] || { echo "Error: Required folder '$OUTF' not found!"; exit 1; }
    [[ -e "$MAPFILE" ]] || { echo "Error: Required file '$MAPFILE' not found!"; exit 1; }

    # Check if folder contains output of part 1.
    BASE1="${OUTF##*/}"  # This will remove the path, keeping just the folder or file name.
    BASE1="${BASE1%_output}"
    BASE1="${BASE1#part1_}"

    BASE2="${OUTF%_output}"
    BASE2="${BASE2#part1_}"

    # Show the output files from part 1. 
    matching_files=($(find -L "$OUTF" -type f -name "${BASE2}*M${mmatchnum}.fq" ! -name "*_BC.M*"))
    # Check if more than one file was found
    if [[ ${#matching_files[@]} -gt 1 ]]; then
        echo "Error: More than one matching file found."
        exit 1
    fi
    
    # Assign the single matching file to FQ if exactly one file is found
    FQ="${matching_files[0]}"
    BC="${FQ/.M${mmatchnum}.fq/_BC.M${mmatchnum}.fq}"

    if [ -z "$FQ" ]; then
        echo "No matching file found for ${BASE2}*M${mmatchnum}.fq. Did you run part 1?"
        exit 1
    fi

    if [ -z "$BC" ]; then
        echo "No matching file found for ${BASE2}*_BC.M${mmatchnum}.fq. Did you run part 1?"
        exit 1
    fi

    # Check if Trim Length Stats Range and Interval are valid
    # Calculate the maximum length of the first 100 reads
    MAX_LENGTH=$(head -n 400 "${FQ}" | awk '{if(NR%4==2) print length($1)}' | sort -nr | head -n 1)
    
    # Step 1: Validate that BEG, INV, and FIN are valid integers
    BEG="${BEG_ARRAY[$i]}"
    BEG=$(echo "$BEG" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    if [[ "$BEG" == "DEFAULT" ]]; then
        BEG="1"
    fi

    FIN="${FIN_ARRAY[$i]}"
    FIN=$(echo "$FIN" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    if [[ "$FIN" == "DEFAULT" ]]; then
        FIN=$((MAX_LENGTH - 5))
    fi

    INV="${INTERVAL_ARRAY[$i]}"
    INV=$(echo "$INV" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    if [[ "$INV" == "DEFAULT" ]]; then
        INV="25"
    fi


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
    ODIR=${OUTF#part1_}
    ODIR=${ODIR%_output}

    SUB="${SUB_ARRAY[$i]}"
    if [[ "$SUB" == "DEFAULT" ]]; then
        ODIR=part2_${ODIR}_output
    else
        ODIR=part2_${ODIR}_${SUB}_subset_output
    fi

    # make the output directory
    if ! mkdir -p "${ODIR}"; then
        echo "Error: Failed to create directory '$ODIR'."
        exit 1
    fi

    # initiate log
    timestamp="$(date +"%y%m%d_%H:%M")"
    output_file="${ODIR}/part2_${timestamp}.log"
    exec > "$output_file"
    exec 2> >(tee -a "$output_file" >&2)

    # create header for log file
    echo "Processing ${OUTF} with mapping file ${MAPFILE}, allowing ${mmatchnum} mismatches.
Range and Interval of Stats: ${BEG}-${MAX_LENGTH}, every ${INV} bases.
 - -- --- ---- ---- --- -- -" | tee /dev/tty

    # Create barcodes.fa file for each JB
    grep -P "^[A-Z]" "${MAPFILE}" | awk '{print ">"$1"\n"$2}' > "${ODIR}/barcodes.fa"
    [[ ! -e "${ODIR}/barcodes.fa" ]] && echo "File ${ODIR}/barcodes.fa was not generated!" && exit 1

    # Demultiplexing
    if ! usearch -fastx_demux "${FQ}" -quiet -index "${BC}" -barcodes "${ODIR}/barcodes.fa" -fastqout "${ODIR}/${BASE1}.M${mmatchnum}.demux.fq"; then
        echo "Error: usearch failed to run."
        exit 1
    fi
    cp "${BC}" "${ODIR}/${BASE1}_BC.M${mmatchnum}.fq"

    # Stats
    usearch -fastq_eestats2 "${ODIR}/${BASE1}.M${mmatchnum}.demux.fq" -quiet -output "${ODIR}/${BASE1}.M${mmatchnum}.eestats.txt" -length_cutoffs "${BEG},${FIN},${INV}"
    # Generate the final 5 stats
    STARTAT=$((MAX_LENGTH - 5))
    INC=1; # Increment value
    usearch -fastq_eestats2 "${ODIR}/${BASE1}.M${mmatchnum}.demux.fq" -quiet -output "${ODIR}/${BASE1}.M${mmatchnum}.eestats.temp.txt" -length_cutoffs ${STARTAT},*,${INC}
    # Merge the two
    tail -n +7 "${ODIR}/${BASE1}.M${mmatchnum}.eestats.temp.txt" >> "${ODIR}/${BASE1}.M${mmatchnum}.eestats.txt"
    rm -rf "${ODIR}/${BASE1}.M${mmatchnum}.eestats.temp.txt"

    echo " - -- --- ---- ---- --- -- -
Final Recommendations
 - -- --- ---- ---- --- -- -
Stats ready. Please view this file and select your trim length accordingly:
${ODIR}/${BASE1}.M${mmatchnum}.eestats.txt
 - -- --- ---- ---- --- -- -
" | tee /dev/tty

done

