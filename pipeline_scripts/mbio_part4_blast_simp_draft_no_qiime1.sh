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
#-j: this flag creates a specialized excel summary output that Dr. Borneman specifically requested. 
    #Runtime will increase, as it requires an analysis examining the top 10 blast hits for each ASV.

# CODE FOLLOWS HERE #

#!/bin/bash

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
james_sum_file_gen=false

while getopts ":d:o:b:e:t:usj" opt; do
  case $opt in
    d) DIR="$OPTARG"
    ;;
    o) OUTDIR="$OPTARG"
    ;;    
    b) blast_file="$OPTARG"
    ;;
    e) EMAIL="$OPTARG"
    ;;  
    u) split_asv_table=true
    ;;
    s) skip_blast=true
    ;;
    j) james_sum_file_gen=true
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

# If blast is not skipped, check for the blast_file arguments
if [ "$skip_blast" = false ]; then
    if [ -z "$blast_file" ]; then
        echo "Usage: $0 -d <directory_path> -o <desired name of output dir> -b <blast parameter file> -e <email@email.com> [-t <filtertax_file> -u -s -j]"
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

if [ ! -e "${output_dir}/asvs/asvs_counts.fa" ]; then
    echo "${output_dir}/asvs/asvs_counts.fa not found!"
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

if [ "$james_sum_file_gen" = true ]; then
    echo "A James Summary File was generated."
fi

echo " - -- --- ---- ---- --- -- -"

mkdir -vp ${output_dir}/asvs/blast

if [ "$skip_blast" = false ]; then
    echo
    echo " - -- --- ---- ---- --- -- -"
    echo "Running BLAST on ASVs"
    echo " - -- --- ---- ---- --- -- -"
    
    # Run BLAST script in the background
    blast_iterator_v2.sh "${output_dir}" "${blast_file}" ${run_type} &
    
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
fi
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

#TODO: ADD ACCESSIONS TO CONTAINER AND ADD STEP TO ADD ANY NEW ONES?
# filter out any lines in the blast file that match the strings in the final.blastout file
grep -v -F -f "${DIR}/likely_environmental_accessions.txt" "${output_dir}/asvs/blast/final.blastout" > "${output_dir}/asvs/blast/filtered.blastout"

# this step parses the blastout into a summary file, keeping only the top bitscore hits. 
./blast_top_hit_parser.py -i "${output_dir}/asvs/blast/filtered.blastout" -o "${output_dir}/asvs/blast/top_hit_summary.txt"
#blast_top_hit_parser.py -i "${output_dir}/asvs/blast/final.blastout" -o "${output_dir}/asvs/blast/top_hit_summary.txt"

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
 
#Mario: Added -o flag which specifies the output directory location
#http://qiime.org/scripts/summarize_taxa.html
summarize_taxa.py -i "${output_dir}/asvs/${OTBL}.biom" -L 2,6,7 -o "${output_dir}/asvs" -a

# will need to chnage if we add more levels
to_process=(
    "${output_dir}/asvs/asv_table_02_add_taxa_L2.biom"
    "${output_dir}/asvs/asv_table_02_add_taxa_L6.biom"
    "${output_dir}/asvs/asv_table_02_add_taxa_L7.biom"
    "${output_dir}/asvs/asv_table_02_add_taxa.biom"
)

for F in "${to_process[@]}"; do
    if [ ! -f "$F" ]; then
        echo "Error: File $F not found."
        exit 1
    fi

    FNAME=$(basename "$F" | sed 's|^asv_table_02_add_taxa||' | sed 's/.biom//')
    ID=$(basename "$F" | sed 's/asv_table_02_add_taxa//' | sed 's/.biom//')

    biom2txt "$F" "${output_dir}/asvs/${OTBL}${ID}.txt"

    if [[ "$split_asv_table" == true ]]; then
        KDOMS=("k__Archaea" "k__Bacteria" "k__Eukaryota")
        for K in "${KDOMS[@]}"; do
            grep -P "(#|$K)" "${output_dir}/asvs/${OTBL}${ID}.txt" > "${output_dir}/asvs/${OTBL}${ID}.${K}.txt"
            NEW_OTBL="${OTBL}${ID}.${K}"
            if ! grep -q "$K" "${output_dir}/asvs/${NEW_OTBL}.txt"; then
                rm "${output_dir}/asvs/${NEW_OTBL}.txt"
            else
                if [[ "${NEW_OTBL}" == *"_L"* ]]; then
                    txt2biom_notax -f "${output_dir}/asvs/${NEW_OTBL}.txt" "${output_dir}/asvs/${NEW_OTBL}.biom"
                else 
                    txt2biom -f "${output_dir}/asvs/${NEW_OTBL}.txt" "${output_dir}/asvs/${NEW_OTBL}.biom"
                fi
                biom_table_math_ops.py -i "${output_dir}/asvs/${NEW_OTBL}.biom" -o "${output_dir}/asvs/${NEW_OTBL}_norm.biom" --normalize2unity
                biom2txt "${output_dir}/asvs/${NEW_OTBL}_norm.biom" "${output_dir}/asvs/${NEW_OTBL}_norm.txt"
            fi
        done
    fi

    biom_table_math_ops.py -i "$F" -o "${output_dir}/asvs/${OTBL}${ID}_norm.biom" --normalize2unity
    biom2txt "${output_dir}/asvs/${OTBL}${ID}_norm.biom" "${output_dir}/asvs/${OTBL}${ID}_norm.txt"
done

# add seqs to L8 
otblfp="${output_dir}/asvs/asv_table_02_add_taxa.txt"
outfp="${output_dir}/asvs/asv_table_03_add_seqs.txt"

Rscript -e "source('${HDIR}/pipeline_helper_functions.R'); add_sequences_to_asv_table('$otblfp', '${output_dir}/asvs/asvs_counts.fa', '$outfp')"

otblfp="${output_dir}/asvs/asv_table_02_add_taxa_norm.txt"
outfp="${output_dir}/asvs/asv_table_03_add_seqs_norm.txt"

Rscript -e "source('${HDIR}/pipeline_helper_functions.R'); add_sequences_to_asv_table('$otblfp', '${output_dir}/asvs/asvs_counts.fa', '$outfp')"

to_process2=($(find "${output_dir}/asvs" -maxdepth 1 -type f -name "*taxa.k*txt"))

for F in "${to_process2[@]}"; do
    if [ ! -f "$F" ]; then
        echo "Error: File $F not found."
        exit 1
    fi

    FNAME=$(basename "$F" | sed 's|^./asv_table_02_add_taxa||')
    otblfp="${F}"
    outfp="${output_dir}/asvs/asv_table_03_add_seqs${FNAME}"
    Rscript -e "source('${HDIR}/pipeline_helper_functions.R'); add_sequences_to_asv_table('$otblfp', '${output_dir}/asvs/asvs_counts.fa', '$outfp')"
done
echo

##TODO: FIX JAMES SUM FILE GEN - ONLY ONE SET OF TAX ASSIGNMENTS NOW. Also IMPROVE MIXED_FAMILY_CHECKER.PY
if [ "$james_sum_file_gen" = true ]; then
    echo " - -- --- ---- ---- --- -- -"
    echo "Creating Summary File"
    echo " - -- --- ---- ---- --- -- -"

    source pymods.sh || { echo "Error: Unable to activate python-pip-modules environment"; exit 1; }

    # Error check for mixed_family_checker.py
    mixed_family_checker.py "${output_dir}/asvs/blast/final.blastout" --email "${EMAIL}" || { echo "Error: mixed_family_checker.py failed"; exit 1; }

    # Move output file with error check
    mv mixed_family_checker_out.txt "${output_dir}/asvs/blast/" || { echo "Error: Unable to move mixed_family_checker_out.txt"; exit 1; }

    # Error check for Rscript
    Rscript -e "source('${HDIR}/pipeline_helper_functions.R'); process_data_and_write_excel('${output_dir}/asvs/asv_table_03_add_seqs_norm.txt', '${output_dir}/asvs/blast/tax_assignments.txt', ${output_dir}/asvs/asv_table_03_add_seqs.txt', '${output_dir}/asvs/blast/mixed_family_checker_out.txt')" || { echo "Error: Rscript failed"; exit 1; }
   
    #Move generated ASVSumamry file
    mv ASV_summary* ${output_dir}/asvs


else
    echo "No Borneman summary file generated, as requested."
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
k__Unassigned;Other
" | tee /dev/tty