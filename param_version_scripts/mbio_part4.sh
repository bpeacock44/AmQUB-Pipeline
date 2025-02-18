#!/bin/bash
# See mbio_tutorial.md for further guidance!

### USAGE ###
# This script accepts two types of input:
# 1. Direct command-line arguments:
#    ./mbio_part4.sh -i output_dir -t 256 -e email@email.com
#    ./mbio_part4.sh  -i output_dir -t 256 -e email@email.com \
#        -b mega-blast -v 0.005 -c classifier.qza -f 0.8 --strategy2 --strategy3 
#    ./mbio_part4.sh  -o output_dir -t 256 -e email@email.com \
#        -c classifier.qza -f disable --strategy2 --strategy3 -s \
           
## Required Flags
# -i The output directory generated in part 3 that you want to add taxonomic assignments to.
# -t The number of threads you have available for blast (and classifier) to use.
# -e Your email for use while querying NCBI online.

## Optional Flags:
# -b The type of BLAST you want to do (discontiguous megablast or mega-blast; blastn is default)
# -v Set the expect threshold for BLAST (0.001 is default)
# -s Skip blast (useful for troubleshooting the rest of the script to avoid waiting for blast to run if it's already been completed)
# -c Assign taxonomy via a classifier in addition to BLAST 
# -f Set confidence interval for classifier (default is 0.7) (you can put "disable" to ignore confidence)
# --strategy2 Also assign taxonomy to OTUs/ASVs generated via Strategy 2 in part 3. 
# --strategy3 Also assign taxonomy to OTUs/ASVs generated via Strategy 3 in part 3. 

# 2. A parameter template file:
#    ./mbio_part4.sh params.txt
#    Where params.txt contains the following rows, comma delimited, no white space between.
#    The labels at the beginning of each row should be the same as below.
#    The following shows what all possible rows:

#        Part 3 Output Folder To Process,output_dir
#        Number of Threads Available,256
#        Email,email@email.com
#        Type of BLAST,megablast
#        Expect Threshold for BLAST,0.005
#        Skip BLAST,true
#        Classifier,classifier.qza
#        Confidence Interval,disable
#        Assign Taxonomy for Strategy 2,true
#        Assign Taxonomy for Strategy 3,true
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
            "Email,"*) EMAIL="${line#*,}" ;;
            "Type of BLAST,"*) TASK="${line#*,}" ;;
            "Expect Threshold for BLAST,"*) EVAL="${line#*,}" ;;
            "Classifier,"*) CFIER="${line#*,}" ;;
            "Confidence Interval,"*) CON="${line#*,}" ;;
            "Assign Taxonomy for Strategy 2,"*) STR2="${line#*,}" ;;
            "Assign Taxonomy for Strategy 3,"*) STR3="${line#*,}" ;;
            "Skip BLAST,"*) skip_blast="${line#*,}" ;;
        esac
    done < "$param_file"
}

# Initialize default values
TASK="blastn"
EVAL="0.001"
STR2=false
STR3=false
skip_blast=false

    echo " - -- --- ---- ---- --- -- -
 █████╗             ██████╗ ██╗   ██╗██████╗ 
██╔══██╗████╗ ████╗██║   ██╗██║   ██║██╔══██╗
███████║██╔████╔██║██║   ██║██║   ██║██████╔╝
██╔══██║██║╚██╔╝██║██║   ██║██║   ██║██╔══██╗
██║  ██║██║ ╚═╝ ██║╚██████╔╝╚██████╔╝██████╔╝
╚═╝  ╚═╝╚═╝     ╚═╝ ╚════██╗ ╚═════╝ ╚═════╝  
PART 4:Taxonomy          ╚═╝
 - -- --- ---- ---- --- -- -"

# Check if the first argument is a parameter file
if [[ -f "$1" ]]; then
    echo "Reading parameters from file: $1
 - -- --- ---- ---- --- -- -"
    parse_parameter_file "$1"
    shift  # Remove the parameter file argument
else
    # Parse long options
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -i|--output_folder) output_dir="$2"; shift 2 ;;
            -t|--threads) NUMTHREADS="$2"; shift 2 ;;
            -e|--email) EMAIL="$2"; shift 2 ;;
            -b|--blast) TASK="$2"; shift 2 ;;
            -v|--eval) EVAL="$2"; shift 2 ;;
            -c|--classifier) CFIER="$2"; shift 2 ;;
            -f|--confidence) CON="$2"; shift 2 ;;
            --strategy2) STR2=true; shift ;;
            --strategy3) STR3=true; shift ;;
            -s|--skip-blast) skip_blast=true; shift ;;
            *) echo "Unknown option: $1" >&2; exit 1 ;;
        esac
    done
fi

# Check for mandatory arguments
if [ -z "$output_dir" ] || [ -z "$NUMTHREADS" ] || [ -z "$EMAIL" ]; then
    echo "Usage: $0 -i <output directory from part 3> -t <number of threads for BLAST> -e <email@email.com> [-s -j]"
    exit 1
fi

# Check if the classifier is a valid QIIME 2 classifier
if [[ "${CFIER}" ]]; then
    if qiime tools peek "${CFIER}" | grep -q 'TaxonomicClassifier'; then
        echo "Qiime2 classifier appears to be correct."
    else
        echo "The classifier doesn't appear to be a valid QIIME 2 classifier."
        exit 1
    fi
fi

# Check if DIS or CON is set without CFIER
if [[ -z "${CFIER}" && ("${CON}" != "") ]]; then
    echo "WARNING: You have changed classifier behavior without actually specifying a classififer." >&2
fi

# Check if email is in correct format
if ! [[ "${EMAIL}" =~ ^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$ ]]; then
    echo "Email is not valid."
    exit 1
fi

# Check if the output directory contains the directories "otus" and/or "asvs"
if [[ -d "${output_dir}/otus" && -d "${output_dir}/asvs" ]]; then
    echo "Error: Both 'otus' and '${typ}s' directories exist in ${output_dir}. Exiting."
    exit 1
elif [[ -d "${output_dir}/otus" ]]; then
    typ="otu"
elif [[ -d "${output_dir}/asvs" ]]; then
    typ="asv"
else
    echo "Error: Neither 'otus' nor 'asvs' directories exist in ${output_dir}. Exiting."
    exit 1
fi

# Check if STR2 or STR3 is defined, then validate the existence of the corresponding directories
if [ "$STR2" = true ]; then
    if [[ ! -d "${output_dir}/${typ}s/STRATEGY2" ]]; then
        echo "Error: STRATEGY2 directory is missing in ${output_dir}/${typ}s. Exiting."
        exit 1
    fi
fi

if [ "$STR3" = true ]; then
    if [[ ! -d "${output_dir}/${typ}s/STRATEGY3" ]]; then
        echo "Error: STRATEGY3 directory is missing in ${output_dir}/${typ}s. Exiting."
        exit 1
    fi
fi

# Define run_type variable.
run_type=local

if [ ! -e "${output_dir}/${typ}s/blast/${typ}s_counts.fa" ]; then
    echo "${output_dir}/${typ}s/blast/${typ}s_counts.fa not found! Part 3 was not completed."
    exit 1
fi

# initiate log
timestamp="$(date +"%y%m%d_%H:%M")"
output_file="${output_dir}/part4_${timestamp}.log"
exec > "$output_file"
exec 2> >(tee -a "$output_file" >&2)

# log header
echo " - -- --- ---- ---- --- -- -
Log file for Part 4 of the Microbiome Pipeline. Processed the following arguments:
Output directory: ${output_dir}
Email of user: ${EMAIL}
Threads for BLAST: ${NUMTHREADS}
BLAST task: ${TASK}
E-value for BLAST: ${EVAL}" | tee /dev/tty

if [ "$skip_blast" = true ]; then
    echo "BLAST will be skipped." | tee /dev/tty
fi

if [[ "${CFIER}" ]]; then
    echo "${CFIER} will also be used to assign taxonomy. These assignments will be saved in ${output_dir}/${typ}s/classifier_output." | tee /dev/tty
fi

if [[ "${DIS}" ]]; then
    echo "Confidence will not be used to assign taxonomy via the Qiime2 classifier." | tee /dev/tty
fi

if [[ "${CON}" ]]; then
    echo "${CON} will be used as the minimum confidence level for classifier taxonomic assignments." | tee /dev/tty
fi

if [ "$STR2" = true ]; then
    echo "Also assigning taxonomy to the tables created using Strategy 2." | tee /dev/tty
fi

if [ "$STR3" = true ]; then
    echo "Also assigning taxonomy to the tables created using Strategy 3." | tee /dev/tty
fi

echo " - -- --- ---- ---- --- -- -" | tee /dev/tty

echo '#!/bin/bash
DATABASE_PATH=/database/nt
OPTS="qseqid sseqid pident length mismatch evalue bitscore staxids stitle qcovs"
INFASTA=$1
MAXTSEQS=$2
blastn -task $TASK -db $DATABASE_PATH -query $INFASTA -max_target_seqs $MAXTSEQS -evalue '$EVAL' -num_threads '"$NUMTHREADS"' -outfmt "7 $OPTS"' > blast.sh

# Make the script executable
chmod +x blast.sh

if [ "$skip_blast" = false ]; then
    
    # Run BLAST script in the background
    blast_iterator.sh "${output_dir}" blast.sh ${run_type} ${typ} &
    
    # Get the process ID of the last background command
    blast_pid=$!
    
    # Wait for the process to finish
    while ps -p $blast_pid > /dev/null; do
        sleep 1
    done
    
    # Check for final.blastout file
    if [ ! -s "${output_dir}/${typ}s/blast/final.blastout" ]; then
        echo "Error: final blast output either does not exist or is empty. Blast has not been completed."
        exit 1
    else
        echo "Blast successfully completed."
    fi
    if [ "$STR2" = true ]; then
        # Run BLAST script in the background
        blast_iterator.sh "${output_dir}/${typ}s/STRATEGY2" blast.sh ${run_type} ${typ} &
    
        # Get the process ID of the last background command
        blast_pid=$!
    
        # Wait for the process to finish
        while ps -p $blast_pid > /dev/null; do
            sleep 1
        done
    
        # Check for final.blastout file
        if [ ! -s "${output_dir}/${typ}s/STRATEGY2/${typ}s/blast/final.blastout" ]; then
            echo "Error: final blast output for strategy 2 either does not exist or is empty. Blast has not been completed."
            exit 1
        else
            echo "Blast successfully completed for strategy 2."
        fi
    fi
    if [ "$STR3" = true ]; then
        # Run BLAST script in the background
        blast_iterator.sh "${output_dir}/${typ}s/STRATEGY3" blast.sh ${run_type} ${typ} &
    
        # Get the process ID of the last background command
        blast_pid=$!
    
        # Wait for the process to finish
        while ps -p $blast_pid > /dev/null; do
            sleep 1
        done
    
        # Check for final.blastout file
        if [ ! -s "${output_dir}/${typ}s/STRATEGY3/${typ}s/blast/final.blastout" ]; then
            echo "Error: final blast output for strategy 3 either does not exist or is empty. Blast has not been completed."
            exit 1
        else
            echo "Blast successfully completed for strategy 3."
        fi
    fi
else
    echo "Skipping BLAST as per user instructions."
fi

process_output_dir() {
    local output_dir=$1

    # Activate the virtual environment containing necessary pip modules
    source pymods.sh || { echo "Error: Unable to activate python-pip-modules environment"; exit 1; }

    # Step 1: Generate a list of likely environmental taxonomic IDs to exclude from the blast file and filter
    likely_env_remove.py "${EMAIL}" "${output_dir}/${typ}s/blast/final.blastout" "${output_dir}/${typ}s/blast/filtered.blastout"

    # Step 2: Parse the blastout into a summary file, keeping only the top bitscore hits
    blast_top_hit_parser.py -i "${output_dir}/${typ}s/blast/filtered.blastout" -o "${output_dir}/${typ}s/blast/top_hit_summary.txt"

    # Step 3: Remove temporary files
    rm -f *.xml

    # Step 4: Assign LCA via BLAST
    assign_LCA_via_blast.py -i "${output_dir}/${typ}s/blast/top_hit_summary.txt" -m "${EMAIL}" -o "${output_dir}/${typ}s/blast/tax_assignments.txt"

    # Step 5: Clean the tax_assignments file using awk
    awk -F'\t' 'BEGIN {OFS="\t"} {sub(/_[0-9]+$/, "", $1); print}' "${output_dir}/${typ}s/blast/tax_assignments.txt" > "${output_dir}/${typ}s/blast/temp.txt"
    mv "${output_dir}/${typ}s/blast/temp.txt" "${output_dir}/${typ}s/blast/tax_assignments.txt"

    # Step 6: Check if the final tax_assignments file exists
    if [ ! -f "${output_dir}/${typ}s/blast/tax_assignments.txt" ]; then
        echo "Error: Output file not found."
        exit 1
    else
        echo "Blast assign taxonomy completed successfully."
    fi

    rm -rf *.xml

    # Step 7: Deactivate the virtual environment
    deactivate

    # Step 8: Source Qiime shell helper functions
    source qiime_shell_helper_functions.sh || { echo "Error: Unable to source Qiime shell helper functions"; exit 1; }

    # Step 9: Run count_taxa_levels command
    count_taxa_levels "${output_dir}/${typ}s/blast/tax_assignments.txt" > "${output_dir}/${typ}s/blast/taxa_levels.txt" || { echo "Error: count_taxa_levels command failed"; exit 1; }

    # Step 10: Conduct taxonomic classification with Qiime2 if classifier is selected
    if [[ "${CFIER}" ]]; then
        qiime tools import \
          --type 'FeatureData[Sequence]' \
          --input-path "${output_dir}/${typ}s/${typ}s.fa" \
          --output-path "${output_dir}/${typ}s/${typ}s.qza"

        # Build the qiime feature-classifier command as an array
        CMD=("qiime" "feature-classifier" "classify-sklearn"
             "--i-classifier" "${CFIER}"
             "--i-reads" "${output_dir}/${typ}s/${typ}s.qza"
             "--o-classification" "${output_dir}/${typ}s/tax_assign_v_classifier.qza"
        )

        # Conditionally add parameters based on environment variables
        [[ -n "${CON}" ]] && CMD+=("--p-confidence" "${CON}")

        # Execute the final command
        echo "Running: ${CMD[@]}"
        "${CMD[@]}"

        # Export the final classification to a text file
        qiime tools export \
          --input-path "${output_dir}/${typ}s/tax_assign_v_classifier.qza" \
          --output-path "${output_dir}/${typ}s/classifier_output"
    fi

    #make normalized version
    qiime tools import \
        --input-path "${output_dir}/${typ}s/${typ}_table_01.biom" \
        --type 'FeatureTable[Frequency]' \
        --input-format BIOMV210Format \
        --output-path "${output_dir}/${typ}s/${typ}_table_01.qza"
    qiime feature-table relative-frequency \
        --i-table "${output_dir}/${typ}s/${typ}_table_01.qza" \
        --o-relative-frequency-table "${output_dir}/${typ}s/${typ}_table_01.norm.qza"
    qiime tools export \
        --input-path "${output_dir}/${typ}s/${typ}_table_01.norm.qza" \
        --output-path "${output_dir}/${typ}s/tmp_dir"
    mv "${output_dir}/${typ}s/tmp_dir/feature-table.biom" "${output_dir}/${typ}s/${typ}_table_01.norm.biom"
    biom2txt "${output_dir}/${typ}s/${typ}_table_01.norm.biom" "${output_dir}/${typ}s/${typ}_table_01.norm.txt"
    rm -rf "${output_dir}/${typ}s/tmp_dir"

    # Define file paths
    norm_file="${output_dir}/${typ}s/${typ}_table_01.norm.txt"
    tax_file="${output_dir}/${typ}s/blast/tax_assignments.txt"
    tax_c_file="${output_dir}/${typ}s/classifier_output/taxonomy.tsv"
    sum_output="${output_dir}/${typ}s/taxonomy_summary_table.tsv"

    rm -rf "${sum_output}"
    if [[ "${CFIER}" ]]; then
        # Run the R script with all arguments (### REMOVE output_dir here!)
        ./mini_summary_file_generator.py ${norm_file} ${tax_file} ${tax_c_file} ${sum_output}
    else 
        ./mini_summary_file_generator.py ${norm_file} ${tax_file} NULL ${sum_output}
    fi
    if [ ! -f ${sum_output} ]; then
        echo "'taxonomy_summary_table.tsv' was not successfully created."
        #exit 1
    fi

}

process_output_dir "$output_dir"

# Finish processing strategy 2 if indicated.
if [ "$STR2" = true ]; then
    process_output_dir "${output_dir}/${typ}s/STRATEGY2"
fi

# Finish processing strategy 3 if indicated.
if [ "$STR3" = true ]; then
    process_output_dir "${output_dir}/${typ}s/STRATEGY3"
fi

echo
echo " - -- --- ---- ---- --- -- -
Final Notes
 - -- --- ---- ---- --- -- -
Taxonomic assignments ready. A summary file was generated: ${sum_output}. 

Before continuing, look over this file to examine your taxonomic assignments.

Note that it contains assignments from BLAST and from the classifier if you elected to generate them. 
However, the BLAST assignments are by default considered the primary assignments in step 5. You can 
opt to use the classifier assignments instead so keep that in mind as you examine your summary.

Use the summary file to create a file with list of all the taxonomic units you want to KEEP in your final tables. 

See tutorial for more guidance." | tee /dev/tty
