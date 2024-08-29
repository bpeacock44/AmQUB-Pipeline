#!/bin/bash
# See mbio_tutorial.md for further guidance!

### USAGE ###
#-d: a working directory, which contains one folder for each of your fastq files named by ID
#-o: the name of your output directory
#-b: the path to your blast script file (only if you are running blast)
#-e: email of the user for NCBI purposes

#Optional arguments:
#-u: universal assay - causes final ASV tables to be split into taxonomic groups prior to normalizing
#-s: skip the blast - skips the blast portion - useful for troubleshooting or re-running taxonomy assignment 
    #steps etc. Note that if -s is enabled, -b is not required.
#-j: this flag creates a specialized excel summary output with more detail about the taxonomic assignments.
    #Runtime will increase, as it requires an analysis examining the top 10 blast hits for each ASV.
#-c: indicate a Qiime2 classifier file you want to use for generating alternative taxonomy.

# Set error handling
set -e

# Custom error handler function
error_handler() {
    local error_message=$1
    echo "Error on line $error_message" | tee /dev/tty
}

# Trap errors and call the error handler
trap 'error_handler "$BASH_COMMAND"' ERR

# ARGUMENTS
split_asv_table=false
skip_blast=false
final_sum_file_gen=false

while getopts ":d:o:b:e:c:usj" opt; do
  case $opt in
    d) DIR="$OPTARG"
    ;;
    o) OUTDIR="$OPTARG"
    ;;    
    b) blast_file="$OPTARG"
    ;;
    e) EMAIL="$OPTARG"
    ;;  
    c) CFIER="$OPTARG"
    ;; 
    u) split_asv_table=true
    ;;
    s) skip_blast=true
    ;;
    j) final_sum_file_gen=true
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

shift $((OPTIND -1))

# Check for mandatory arguments
if [ -z "$DIR" ] || [ -z "$OUTDIR" ] || [ -z "$EMAIL" ]; then
    echo "Usage: $0 -d <directory_path> -o <desired name of output dir> -e <email@email.com> [-b <blast parameter file> -t <filtertax_file> -u -s -j]"
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

# If blast is not skipped, check for the blast_file arguments
if [ "$skip_blast" = false ]; then
    if [ -z "$blast_file" ]; then
        echo "Usage: $0 -d <directory_path> -o <desired name of output dir> -b <blast parameter file> -e <email@email.com> [-t <filtertax_file> -u -s -j]"
        exit 1
    elif [ ! -f "$blast_file" ]; then
        echo "Error: BLAST file '$blast_file' does not exist."
        exit 1
    fi
fi

# Check if email is in correct format
if ! [[ "${EMAIL}" =~ ^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$ ]]; then
    echo "Email is not valid."
    exit 1
fi

# Define initial paths and run_type variable.
output_dir="${DIR}/${OUTDIR}"
run_type=local

echo " - -- --- ---- ---- --- -- -"
echo "Checking for input files"
echo " - -- --- ---- ---- --- -- -"

if [ ! -e "${output_dir}/asvs/blast/asvs_counts.fa" ]; then
    echo "${output_dir}/asvs/blast/asvs_counts.fa not found! Part 3 was not completed."
    exit 1
fi

# initiate log
timestamp="$(date +"%Y%m%d_%H:%M:%S")"
output_file="${DIR}/part4_${timestamp}.log"
exec > "$output_file" 2>&1

# log header
echo " - -- --- ---- ---- --- -- -"
echo "Log file for Part 4 of the Microbiome Pipeline. Processed the following arguments:
Working directory: ${DIR}
Output directory: ${OUTDIR}
Email of user: ${EMAIL}"

if [ "$skip_blast" = true ]; then
    echo "BLAST was skipped."
else
    echo "Blast run file: ${blast_file}"
fi

if [ "$split_asv_table" = true ]; then
    echo "Final ASV tables were split into three domains of life (for universal assay data)"
fi

if [ "$final_sum_file_gen" = true ]; then
    echo "A Final Summary File was generated."
fi

if [[ "${CFIER}" ]]; then
    echo "${CFIER} will also be used to assign taxonomy. These assignments will be saved in ${output_dir}/asvs/classifier_output. They will also be included in the Final Summary File if generated." 
fi

echo " - -- --- ---- ---- --- -- -"

mkdir -vp ${output_dir}/asvs/blast

if [ "$skip_blast" = false ]; then
    echo
    echo " - -- --- ---- ---- --- -- -"
    echo "Running BLAST on ASVs"
    echo " - -- --- ---- ---- --- -- -"
    
    # Run BLAST script in the background
    blast_iterator.sh "${output_dir}" "${blast_file}" ${run_type} &
    
    # Get the process ID of the last background command
    blast_pid=$!
    
    # Wait for the process to finish
    while ps -p $blast_pid > /dev/null; do
        sleep 1
    done
    
    # Check for final.blastout file
    if [ ! -s "${output_dir}/asvs/blast/final.blastout" ]; then
    echo "Error: final blast output either does not exist or is empty. Blast has not been completed."
    exit 1
else
    echo "Skipping BLAST as per user instructions."
fi

echo
echo " - -- --- ---- ---- --- -- -"
echo "Determining Likely Taxonomy of ASVs"
echo " - -- --- ---- ---- --- -- -"

#Mario: Created helper script that activates python virtual environment containing the necessary pip modules
#for the next steps.
source pymods.sh || { echo "Error: Unable to activate python-pip-modules environment"; exit 1; }

# generate a list of likely environmental taxonomic IDs to exclude from the blast file.
retrieve_taxonomy() {
    local query="$1"
    local output_file="$2"

    local max_attempts=3
    local attempt=1
    local success=false

    while [ $attempt -le $max_attempts ]; do
        echo "Attempt $attempt: Retrieving taxonomy for $output_file"
        esearch -db taxonomy -query "$query" | efetch -format uid > "$output_file"
        if [ $? -eq 0 ]; then
            success=true
            break
        else
            echo "Attempt $attempt failed. Retrying in 5 seconds..."
            sleep 5
            ((attempt++))
        fi
    done

    if ! $success; then
        echo "Failed to retrieve taxonomy for $output_file after $max_attempts attempts. This may be because you're requesting too many in quick succession. Run part 4 again with the s flag to skip the BLAST."
        exit 1
    fi
}

FILE="${output_dir}/likely_env_taxids_removed.txt"
retrieve_taxonomy "\"environmental samples\"[subtree] OR \"Environmental Samples\"[subtree] OR \"unclassified\"[subtree] OR \"Unclassified\"[subtree] OR \"uncultured\"[subtree] OR \"Uncultured\"[subtree]" "$FILE"

# filter out likely environmental taxonomic IDs
likely_env_remove.py "${output_dir}/likely_env_taxids_removed.txt" "${output_dir}/asvs/blast/filtered.blastout" "${output_dir}/asvs/blast/filtered.final.blastout"

# this step parses the blastout into a summary file, keeping only the top bitscore hits. 
blast_top_hit_parser.py -i "${output_dir}/asvs/blast/filtered.final.blastout" -o "${output_dir}/asvs/blast/top_hit_summary.txt"

rm -f *.xml

assign_LCA_via_blast.py -i "${output_dir}/asvs/blast/top_hit_summary.txt" -m ${EMAIL} -o "${output_dir}/asvs/blast/tax_assignments.txt"

if [ ! -f "${output_dir}/asvs/blast/tax_assignments.txt" ]; then
    echo "Error: Output file not found."
    exit 1
else
    echo "Blast assign taxonomy completed successfully."
fi

rm -rf *.xml

deactivate

echo
echo " - -- --- ---- ---- --- -- -"
echo "Printing Taxa Levels"
echo " - -- --- ---- ---- --- -- -"

# Source Qiime shell helper functions
source qiime_shell_helper_functions.sh || { echo "Error: Unable to source Qiime shell helper functions"; exit 1; }

# Run count_taxa_levels command
count_taxa_levels "${output_dir}/asvs/blast/tax_assignments.txt" > "${output_dir}/asvs/blast/taxa_levels.txt" || { echo "Error: count_taxa_levels command failed"; }

echo
echo " - -- --- ---- ---- --- -- -"
echo "Adding Taxa and Sequences to ASV Tables"
echo " - -- --- ---- ---- --- -- -"

#add taxa to ASV table
OTBL=asv_table_01
biomAddObservations "${output_dir}/asvs/${OTBL}.biom" "${output_dir}/asvs/asv_table_02_add_taxa.biom" "${output_dir}/asvs/blast/tax_assignments.txt"

# create three additional taxonomic levels of ASV tables
OTBL="asv_table_02_add_taxa"

# first, convert tables to qza format
biom convert \
  -i "${output_dir}/asvs/${OTBL}.biom" \
  --to-json \
  -o "${output_dir}/asvs/${OTBL}.v100.biom"
mv "${output_dir}/asvs/${OTBL}.v100.biom" "${output_dir}/asvs/${OTBL}.biom"

qiime tools import \
  --input-path "${output_dir}/asvs/${OTBL}.biom" \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV100Format \
  --output-path "${output_dir}/asvs/${OTBL}.qza"

tail -n +2 "${output_dir}/asvs/blast/tax_assignments.txt" | cut -f1,2 | sed '1s/^/#ASVID\tTaxon\n/' > temp.txt 

# convert taxonomy assignment file to qza format
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path temp.txt \
  --output-path "${output_dir}/asvs/blast/taxonomy.qza"
rm -rf temp.txt

# generate levels
to_process=(
    "${output_dir}/asvs/asv_table_02_add_taxa_L2"
    "${output_dir}/asvs/asv_table_02_add_taxa_L6"
    "${output_dir}/asvs/asv_table_02_add_taxa_L7"
    "${output_dir}/asvs/asv_table_02_add_taxa"
)

for F in ${to_process[@]}; do
    if [[ "$F" =~ _L([0-9]+) ]]; then
        num="${BASH_REMATCH[1]}"
        qiime taxa collapse \
          --i-table "${output_dir}/asvs/asv_table_02_add_taxa.qza" \
          --i-taxonomy "${output_dir}/asvs/blast/taxonomy.qza" \
          --p-level ${num} \
          --o-collapsed-table "${F}.qza"
    fi
    # generate normalized version
    qiime feature-table relative-frequency \
      --i-table "${F}.qza" \
      --o-relative-frequency-table "${F}.norm.qza"
    # convert raw qza files to biom
    qiime tools export \
      --input-path "${F}.qza" \
      --output-path "${output_dir}/asvs/tmp_dir"
    mv "${output_dir}/asvs/tmp_dir/feature-table.biom" "${F}.biom"
    # repeat with norm
    qiime tools export \
      --input-path "${F}.norm.qza" \
      --output-path "${output_dir}/asvs/tmp_dir"
    mv "${output_dir}/asvs/tmp_dir/feature-table.biom" "${F}.norm.biom"
    rm -rf "${output_dir}/asvs/tmp_dir"
done

# re-add observations if necessary
biomAddObservations "${output_dir}/asvs/asv_table_02_add_taxa.biom" "${output_dir}/asvs/temp.tmp" "${output_dir}/asvs/blast/tax_assignments.txt"
mv "${output_dir}/asvs/temp.tmp" "${output_dir}/asvs/asv_table_02_add_taxa.biom"

biomAddObservations "${output_dir}/asvs/asv_table_02_add_taxa.norm.biom" "${output_dir}/asvs/temp.tmp" "${output_dir}/asvs/blast/tax_assignments.txt"
mv "${output_dir}/asvs/temp.tmp" "${output_dir}/asvs/asv_table_02_add_taxa.norm.biom"

# split into 3 domains if indicated and normalize resulting tables
if [[ "$split_asv_table" == true ]]; then
    for F in "${to_process[@]}"; do
        biom2txt "${F}.biom" "${F}.txt"
        KDOMS=("k__Archaea" "k__Bacteria" "k__Eukaryota")
        for K in "${KDOMS[@]}"; do
            grep -P "(#|$K)" "${F}.txt" > "${F}.${K}.txt"
            if [[ "$F" =~ _L([0-9]+) ]]; then
                # Process files with _L in the name (no taxonomy metadata)
                txt2biom_notax "${F}.${K}.txt" "${F}.${K}.biom"
            else
                # Process files without _L in the name (include taxonomy metadata)
                txt2biom "${F}.${K}.txt" "${F}.${K}.biom"
            fi
            qiime tools import \
              --type 'FeatureTable[Frequency]' \
              --input-path "${F}.${K}.biom" \
              --output-path "${F}.${K}.qza"
            qiime feature-table relative-frequency \
              --i-table "${F}.${K}.qza" \
              --o-relative-frequency-table "${F}.${K}.norm.qza"
            qiime tools export \
              --input-path "${F}.${K}.norm.qza" \
              --output-path "${output_dir}/asvs/tmp_dir"
            mv "${output_dir}/asvs/tmp_dir/feature-table.biom" "${F}.${K}.norm.biom"
            rm -rf "${output_dir}/asvs/tmp_dir"
        done  
    done
fi

for F in "${output_dir}/asvs/"*.biom; do
    txt_file="${F%.biom}.txt"
    if [ ! -f "$txt_file" ]; then
        biom2txt "$F" "$txt_file"
    fi
done

# add seqs to L8 
otblfp="${output_dir}/asvs/asv_table_02_add_taxa.txt"
outfp="${output_dir}/asvs/asv_table_03_add_seqs.txt"

Rscript -e "source('${HDIR}/pipeline_helper_functions.R'); add_sequences_to_asv_table('$otblfp', '${output_dir}/asvs/blast/asvs_counts.fa', '$outfp')"

otblfp="${output_dir}/asvs/asv_table_02_add_taxa.norm.txt"
outfp="${output_dir}/asvs/asv_table_03_add_seqs.norm.txt"

Rscript -e "source('${HDIR}/pipeline_helper_functions.R'); add_sequences_to_asv_table('$otblfp', '${output_dir}/asvs/blast/asvs_counts.fa', '$outfp')"

to_process2=($(find "${output_dir}/asvs" -maxdepth 1 -type f -name "*taxa.k*txt"))

for F in "${to_process2[@]}"; do
    if [ ! -f "$F" ]; then
        echo "Error: File $F not found."
        exit 1
    fi
    FNAME=$(basename "$F" | sed 's|^./asv_table_02_add_taxa||')
    otblfp="${F}"
    outfp="${output_dir}/asvs/asv_table_03_add_seqs${FNAME}"
    Rscript -e "source('${HDIR}/pipeline_helper_functions.R'); add_sequences_to_asv_table('$otblfp', '${output_dir}/asvs/blast/asvs_counts.fa', '$outfp')"
done
echo

# if classifier selected, conduct taxonomic classification of asvs using Qiime2
if [[ "${CFIER}" ]]; then
    qiime tools import \
      --type 'FeatureData[Sequence]' \
      --input-path "${output_dir}/asvs/asvs.fa" \
      --output-path "${output_dir}/asvs/asvs.qza"
    qiime feature-classifier classify-sklearn \
      --i-classifier ${CFIER} \
      --i-reads "${output_dir}/asvs/asvs.qza" \
      --o-classification "${output_dir}/asvs/tax_assign_v_classifier.qza" 
    # Export the final classification to a text file
    qiime tools export \
      --input-path "${output_dir}/asvs/tax_assign_v_classifier.qza" \
      --output-path "${output_dir}/asvs/classifier_output"
fi

if [ "$final_sum_file_gen" = true ]; then
    echo " - -- --- ---- ---- --- -- -"
    echo "Creating Summary File"
    echo " - -- --- ---- ---- --- -- -"

    source pymods.sh || { echo "Error: Unable to activate python-pip-modules environment"; exit 1; }

    rm -rf "${output_dir}/asvs/blast/mixed_family_checker_out.txt"
    # Run the mixed_family_checker.py script with error checking
    mixed_family_checker.py "${output_dir}/asvs/blast/final.blastout" --email "${EMAIL}" || { echo "Error: mixed_family_checker.py failed"; exit 1; }

    # Move the output file with error checking
    mv mixed_family_checker_out.txt "${output_dir}/asvs/blast/" || { echo "Error: Unable to move mixed_family_checker_out.txt"; exit 1; }

    # Define file paths
    norm_file="${output_dir}/asvs/asv_table_03_add_seqs.norm.txt"
    raw_file="${output_dir}/asvs/asv_table_03_add_seqs.txt"
    tax_file="${output_dir}/asvs/blast/tax_assignments.txt"
    tax_c_file="${output_dir}/asvs/classifier_output/taxonomy.tsv"
    mixed_file="${output_dir}/asvs/blast/mixed_family_checker_out.txt"

    rm -rf ASV_summary_table.tsv "${output_dir}/asvs/ASV_summary_table.tsv"
    if [[ "${CFIER}" ]]; then
        # Run the R script with all arguments
        Rscript -e "source('${HDIR}/pipeline_helper_functions.R'); process_data_and_write_summary('${norm_file}', '${raw_file}', '${tax_file}', '${tax_c_file}', '${mixed_file}')"
    else 
        Rscript -e "source('${HDIR}/pipeline_helper_functions.R'); process_data_and_write_summary('${norm_file}', '${raw_file}', '${tax_file}', NULL, '${mixed_file}')"
    fi
    if [ ! -f "ASV_summary_table.tsv" ]; then
        echo "'ASV_summary_table.tsv' was not successfully created."
        #exit 1
    fi

    # Move the generated ASV summary file with error checking
    mv ASV_summary_table.tsv "${output_dir}/asvs" 
else
    echo "No Final Summary File generated, as requested."
fi

echo
echo "All ASV tables have been generated." | tee /dev/tty
echo " - -- --- ---- ---- --- -- -"  | tee /dev/tty
echo "Final Recommendations"  | tee /dev/tty
echo " - -- --- ---- ---- --- -- -"  | tee /dev/tty
echo "Sometimes ASVs can be primarily associated with your controls (likely a source of contamination,
which may be from other samples or even from the individual building the library in the first place). Make sure 
to check your ASV tables for these ASVs - just look at your control columns and see if any of the ASVs are 
relatively high in them and not present in the regular samples. These ASVs should be removed before further 
analysis.

Additionally, there may be ambiguous taxonomic assignments in your ASV table. This may not be problematic if the 
ambigulously assigned ASV isn't abundant, but if it is then you should manually blast the sequence and determine
taxonomy for yourself on NCBI prior to further analyses. If you make changes to taxonomy, note that your L2, L6, 
and L7 ASV tables will likely change as well and you may want to regenerate them. 

Here are some examples of ambiguous taxa:

k__Bacteria;Other
k__Fungi;Other
k__Eukaryota;Other
k__Bacteria;p__unclassified_Bacteria
k__Fungi;p__unclassified_Fungi
k__Eukaryota;p__unclassified_Eukaryota
k__Bacteria_OR_k__unclassified_;Other
k__Fungi_OR_k__unclassified_;Other
k__Eukaryota_OR_k__unclassified_;Other
k__Unassigned;Other" | tee /dev/tty