#!/bin/bash 

# This script will run BLAST with 1000 max seqs and then check to see if the results are satisfcatory to determine LCA
# for all ASVs that are greater than 1% the size of the largest ASV. 
# In this context, that means that the results reach a point where the bitscore begins to go down so you know you've 
# collected all the highest bitscore results available. If you have not reached this point, BLAST will run again with 
# 30000 max seqs for those ASVs only. 
# If LCA still cannot be determined, then the script will proceed regardless of whether the results are satisfactory
# and all BLAST results will be used.

set -e

# error handling function
error_handler() {
    local error_message=$1
    echo "Error on line $error_message" | tee /dev/tty
}

# Trap errors and call the error handler
trap 'error_handler "$BASH_COMMAND"' ERR

# Define initial values
reblast_iteration="rb0"
maxseqs=1000

DIR=$1
BLAST_FILE=$2
RUN_TYPE=$3
TAX_TYPE=$4
NEW_ONLY_FASTA=${5:-}
REUSED_BLAST=${6:-}

timestamp="$(date +"%y%m%d_%H:%M")"
output_file="${DIR}/blast.${timestamp}.log"
exec > "$output_file" 2>&1

first_run=true  # Add a flag for the first run

criteria_met() {
    if [ "$first_run" == true ]; then
        echo "Skipping criteria check for the first run."
        return 0
    fi

    echo "Checking criteria..."
    reblast_file="${DIR}/${TAX_TYPE}s/blast/${reblast_iteration}.fasta"

    # Check if the reblast file has any lines
    if [ ! -s "$reblast_file" ] || [ "$(wc -l < "$reblast_file")" -eq 0 ]; then
        echo "${reblast_file} is empty. Ending the loop."
        return 1 # End if the reblast file is empty
    else
        echo "${reblast_file} has content. Continuing..."
    fi

    return 0 # Continue if reblast file has content
}

while criteria_met; do
    job_ids=()
    echo "Starting first batch of jobs..."

    if "$first_run"; then
        input_file="${NEW_ONLY_FASTA:-${DIR}/${TAX_TYPE}s/blast/${TAX_TYPE}s_counts.fa}"
    else
        input_file="${DIR}/${TAX_TYPE}s/blast/${reblast_iteration}.fasta"
    fi

    if [[ "$RUN_TYPE" == local ]]; then
        # Run locally
        if [[ -n "$NEW_ONLY_FASTA" && ! -s "$input_file" ]]; then
            : > "${DIR}/${TAX_TYPE}s/blast/${maxseqs}.${reblast_iteration}.blastout"
        else
            ${BLAST_FILE} "${input_file}" "${maxseqs}" > "${DIR}/${TAX_TYPE}s/blast/${maxseqs}.${reblast_iteration}.blastout"
        fi
        echo "BLAST completed locally for $input_file"
    else
        # Submit as sbatch job
        if [[ -n "$NEW_ONLY_FASTA" && ! -s "$input_file" ]]; then
            : > "${DIR}/${TAX_TYPE}s/blast/${maxseqs}.${reblast_iteration}.blastout"
        else
        echo sbatch -o "${DIR}/${TAX_TYPE}s/blast/${maxseqs}.${reblast_iteration}.blastout" ${BLAST_FILE} "${input_file}" "${maxseqs}" 
        job_id=$(sbatch -o "${DIR}/${TAX_TYPE}s/blast/${maxseqs}.${reblast_iteration}.blastout" ${BLAST_FILE} "${input_file}" "${maxseqs}" | awk '{print $NF}')
        echo "Blast job submitted with ID: $job_id for $input_file"
        job_ids+=("$job_id")

        # Wait for completion; use sacct for final exit status when available.
        while [[ -n "$(squeue -j "$job_id" -h -o '%i' 2>/dev/null)" ]]; do
            sleep 10
        done
        if command -v sacct >/dev/null 2>&1; then
            state=$(sacct -j "$job_id" -n -X -o State | head -1 | xargs)
            [[ "$state" == COMPLETED* ]] || { echo "SLURM BLAST job $job_id failed: $state"; exit 1; }
        fi
        echo "Job $job_id has completed."
        fi
    fi

    next_reblast_value() {
        local current="$1"
        case $current in
            "rb0") echo "rb1" ;;
            "rb1") echo "completed" ;;
            *) echo "error"; exit 1 ;;
        esac
    }

    next_value=$(next_reblast_value "$reblast_iteration")
    echo $next_value
    if [[ "$next_value" == "completed" ]]; then
        echo "All reblast iterations have been completed."
        break
    elif [[ "$next_value" == "error" ]]; then
        echo "Unknown reblast iteration value: $reblast_iteration"
        exit 1
    fi

    bout="${DIR}/${TAX_TYPE}s/blast/${maxseqs}.${reblast_iteration}.blastout"
    outfile="${DIR}/${TAX_TYPE}s/blast/${maxseqs}.${reblast_iteration}.blastout.not_enough_hits.txt"
    if [[ "$first_run" == true && -n "$REUSED_BLAST" ]]; then
        cat "$REUSED_BLAST" "$bout" > "${bout}.combined"
        bout="${bout}.combined"
    fi
    "$(dirname "$(readlink -f "$0")")/reblast_check.pl" "$bout" "$outfile" 1000

    # Get the number of reads contributing to the biggest ASV. 
    # Only ASVs that are at least 1% of the largest ASV will be re-blasted if they need it. 
    total=$(awk '/^>/{n=split($1,a,"_"); if (a[n]+0>max) max=a[n]+0} END{print max+0}' "${DIR}/${TAX_TYPE}s/blast/${TAX_TYPE}s_counts.fa")
    echo "Determining which taxonomic units are worth re-blasting."

    not_enough_hits_file="${DIR}/${TAX_TYPE}s/blast/${maxseqs}.${reblast_iteration}.blastout.not_enough_hits.txt"
    big_ASVs_file="${DIR}/${TAX_TYPE}s/blast/${maxseqs}.${reblast_iteration}.blastout.not_enough_hits.big_TUs.txt"
    reblast_seqs_file="${DIR}/${TAX_TYPE}s/blast/${next_value}.fasta"
    
    # Check if the not_enough_hits_file is empty
    if [ -s "$not_enough_hits_file" ]; then
        # Proceed with grep and subsequent commands
        # Counts are encoded in FASTA IDs as ASV_12345; use the last underscore.
        # Compare only identifiers flagged by reblast_check.pl.
        awk -v total="$total" 'NR==FNR {wanted[$1]=1; next}
            /^>/ {id=substr($1,2); n=split(id,a,"_"); count=a[n]+0;
                   if (id in wanted && total>0 && count/total>=0.01) print id}' \
            "$not_enough_hits_file" "${DIR}/${TAX_TYPE}s/blast/${TAX_TYPE}s_counts.fa" > "$big_ASVs_file"

        echo "Extracting reblast seqs."
        if [[ -s "$big_ASVs_file" ]]; then
            seqkit grep -n -f "$big_ASVs_file" "${DIR}/${TAX_TYPE}s/blast/${TAX_TYPE}s_counts.fa" -o "$reblast_seqs_file"
        else
            : > "$reblast_seqs_file"
        fi
    else
        # If the not_enough_hits_file is empty, create an empty file for reblast_seqs
        : > "$reblast_seqs_file"
    fi

    case $reblast_iteration in
        "rb0") maxseqs=30000; reblast_iteration="rb1" ;;
        "rb1") break ;;
    esac

    first_run=false
    echo "First loop completed."
done

echo "Merging all blastout files."

blast_dir="${DIR}/${TAX_TYPE}s/blast"
initial=("${blast_dir}/1000.rb0.blastout")
if [[ -n "$REUSED_BLAST" ]]; then
    initial=("$REUSED_BLAST" "${initial[@]}")
fi
python3 "$(dirname "$(readlink -f "$0")")/merge_blast.py" \
    --initial "${initial[@]}" \
    --expanded "${blast_dir}/30000.rb1.blastout" \
    --output "${blast_dir}/final.blastout"

rm -f ${DIR}/${TAX_TYPE}s/blast/1000.rb0.blastout* ${DIR}/${TAX_TYPE}s/blast/30000.rb0.blastout* ${DIR}/${TAX_TYPE}s/blast/rb1.fasta  


