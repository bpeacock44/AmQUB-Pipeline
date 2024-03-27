#!/bin/bash 

# This script will run BLAST with 5000 max seqs and then check to see if the results are satisfcatory to determine LCA
# for all ASVs that are greater than 1% the size of the largest ASV. 
# In this context, that means that the results reach a point where the bitscore begins to go down so you know you've 
# collected all the highest bitscore results available. If you have not reached this point, BLAST will run again with 
# 30000 max seqs for those ASVs only. 
# If LCA still cannot be determined, then the script will proceed regardless of whether the results are satisfactory
# and all BLAST results will be used.

set -e

HDIR=/home/bpeacock_ucr_edu/real_projects/PN94_singularity_of_microbiome_pipeline/targeted_microbiome_via_blast/helper_functions

# Set strict mode
set -euo pipefail

# Define initial values
reblast_iteration="rb0"
maxseqs=5000

DIR=$1
BLAST_FILE=$2
RUN_TYPE=$3

timestamp="$(date +"%Y%m%d_%H:%M:%S")"
output_file="${DIR}/blast.${timestamp}.log"
exec > "$output_file" 2>&1

first_run=true  # Add a flag for the first run

criteria_met() {
    if [ "$first_run" == true ]; then
        echo "Skipping criteria check for the first run."
        return 0
    fi

    echo "Checking criteria..."
    reblast_file="${DIR}/asvs/rep_set/${reblast_iteration}.fasta"

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
        input_file="${DIR}/asvs/rep_set/seqs_chimera_filtered_ASVs.fasta"
    else
        input_file="${DIR}/asvs/rep_set/${reblast_iteration}.fastq"
    fi

    if [[ "$RUN_TYPE" == local ]]; then
        # Run locally
        ${BLAST_FILE} "${input_file}" "${maxseqs}" > "${DIR}/asvs/rep_set/${maxseqs}.${reblast_iteration}.blastout"
        echo "BLAST completed locally for $input_file"
    else
        # Submit as sbatch job
        echo sbatch -o "${DIR}/asvs/rep_set/${maxseqs}.${reblast_iteration}.blastout" ${BLAST_FILE} "${input_file}" "${maxseqs}" 
        job_id=$(sbatch -o "${DIR}/asvs/rep_set/${maxseqs}.${reblast_iteration}.blastout" ${BLAST_FILE} "${input_file}" "${maxseqs}" | awk '{print $NF}')
        echo "Blast job submitted with ID: $job_id for $input_file"
        job_ids+=("$job_id")

        # Wait until the job starts running
        while true; do
            # Check the job status
            job_status=$(squeue -j "$job_id" -h -o "%t" 2>/dev/null)

            # Break the loop if the job is in the "running" state
            if [ "$job_status" == "R" ]; then
                echo "Job $job_id is now running."
                break
            fi

            # Sleep for a short interval before checking again
            sleep 10
        done

        # Wait for the submitted job to finish
        while squeue -j "$job_id" -h &>/dev/null; do
            sleep 10
        done

        echo "Job $job_id has completed."
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
        exit 0
    elif [[ "$next_value" == "error" ]]; then
        echo "Unknown reblast iteration value: $reblast_iteration"
        exit 1
    fi

    bout="${DIR}/asvs/rep_set/${maxseqs}.${reblast_iteration}.blastout"
    outfile="${DIR}/asvs/rep_set/${maxseqs}.${reblast_iteration}.blastout.not_enough_hits.txt"
    ${HDIR}/reblast_check.pl ${bout} ${outfile}
    echo $bout 
    echo $outfile

    # Get the number of reads contributing to the biggest ASV
    total=$(awk 'NR==1{print $NF}' ${DIR}/asvs/rep_set/seqs_chimera_filtered_ASVs.fasta)
    echo "Determining which ASVs are worth blasting."
    echo ${total}

    not_enough_hits_file="${DIR}/asvs/rep_set/${maxseqs}.${reblast_iteration}.blastout.not_enough_hits.txt"
    big_ASVs_file="${DIR}/asvs/rep_set/${maxseqs}.${reblast_iteration}.blastout.not_enough_hits.big_ASVs.txt"
    reblast_seqs_file="${DIR}/asvs/rep_set/${next_value}.fasta"
    
    # Check if the not_enough_hits_file is empty
    if [ -s "$not_enough_hits_file" ]; then
        # Proceed with grep and subsequent commands
        grep -w -f "$not_enough_hits_file" ${DIR}/asvs/rep_set/seqs_chimera_filtered_ASVs.fasta | \
        awk -v total="$total" '{
            # Remove ">" from the first column
            gsub(">", "", $1)
    
            # Calculate the ratio and store the result after the white space
            ratio = $2 / total
    
            # Print lines where the ratio is 0.01 or less
            if (ratio <= 0.01) {
                print $1, $2, ratio
            }
        }' | \
        awk '$3 > 0.01' > "$big_ASVs_file"
    
        # Additional commands
        awk '/^>/{sub(/ .*/, "");}1' "${DIR}/asvs/rep_set/seqs_chimera_filtered_ASVs.fasta" > ${DIR}/asvs/rep_set/modified_seqs.fasta
    
        echo "Extracting reblast seqs."
        # Extract the fasta sequences needing a reblast
        seqkit grep -n -f "$big_ASVs_file" ${DIR}/asvs/rep_set/modified_seqs.fasta -o "$reblast_seqs_file"
    else
        # If the not_enough_hits_file is empty, create an empty file for reblast_seqs
        touch "$reblast_seqs_file"
    fi

    case $reblast_iteration in
        "rb0") maxseqs=30000; reblast_iteration="rb1" ;;
        "rb1") echo "Completed all reblast iterations. Any ASVs that may require further reBLASTing are stored in ${DIR}/asvs/rep_set/${next_value}.fasta"; exit 0 ;;
    esac

    first_run=false
    echo "First loop completed."
done

echo "Merging all blastout files."

rm -f ${DIR}/asvs/rep_set/final.blastout
cat ${DIR}/asvs/rep_set/5000.rb0.blastout | grep -v "# BLAST processed" >> ${DIR}/asvs/rep_set/final.blastout 
if [ -e "${DIR}/asvs/rep_set/30000.rb1.blastout" ]; then
    cat "${DIR}/asvs/rep_set/30000.rb1.blastout" | grep -v "# BLAST processed" >> "${DIR}/asvs/rep_set/final.blastout"
fi

echo "# BLAST processed" >> ${DIR}/asvs/rep_set/final.blastout 

