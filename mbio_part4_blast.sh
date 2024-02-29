#!/bin/bash

### USAGE ###
# This script expects to be given at least 4 aguments:
# -d: a working directory, which contains one folder for each of your fastq files named by ID
# -o: the name of your output directory
# -b: the path to your blast script file
# -r: the type of blast run you want to do (local or slurm)
# -e: email of the user for NCBI purposes

# Optional arguments:
# -t This is a filter file - it will have four tab-delimited columns. 
# Include any taxonomic groups that you want to preferentially keep or reject as well as their taxonomic ID. 
# ALL taxonomies included under these taxonomic IDs will be treated accordingly so check NCBI
# and make sure. If you are doing a universal assay, do not include the -t flag and DO include the -u flag.
# -u: universal assay - causes final ASV tables to be split into taxonomic groups prior to normalizing
# -s: skip the blast - skips the blast portion - useful for troubleshooting or re-running taxonomy assignment steps etc.

# Examples:
# mbio_part4.sh -d /path/to/dir -o test1_out -b /path/to/blast.sh -e email@email.com -r slurm -t ${MDIR}/filterfile.txt 
# mbio_part4.sh -d /path/to/dir -o test2_out -b /path/to/blast.sh -e email@email.com -r slurm -t ${MDIR}/filterfile.txt -m 1 
# mbio_part4.sh -d /path/to/dir -o test3_out -b /path/to/blast.sh -e email@email.com -r local -s
# mbio_part4.sh -d /path/to/dir -o test4_out -b /path/to/blast.sh -e email@email.com -r local -u

### INPUT ###
# This script follows part 3, which must be completed first. 
# The output file will have already been generated in part 3.



# ########## SLURM BLAST FILE EXAMPLE (-b) ########## 
##!/bin/bash 
##SBATCH -p i128
##SBATCH -c 128
## any other parameters or modules needed
#module load blast-plus

##<>#<>#<>#<>#<>
## YOU MUST SET THESE:
##<>#<>#<>#<>#<>
#DATABASE_PATH=/sw/dbs/blast_db_download/nt
#NUMTHREADS=128

##<>#<>#<>#<>#<>
## GENERALLY DON'T CHANGE THESE:
##<>#<>#<>#<>#<>
#OPTS="qseqid sseqid pident length mismatch evalue bitscore staxids stitle qcovs"
#TASK=blastn
#INFASTA=$1
#MAXTSEQS=$2  
#EVAL=0.001
# blastn -task $TASK -db $DATABASE_PATH -query $INFASTA -max_target_seqs $MAXTSEQS -evalue $EVAL -num_threads $NUMTHREADS -outfmt "7 $OPTS" 
# ########## SLURM BLAST FILE EXAMPLE (-b) ########## 



# If running a local blast, use the same format but don't include the SBATCH lines.



# ########## FILTER FILE EXAMPLE (-t) ########## 
# If I am doing fungal ITS taken from a plant sample, then the file might include:
# Name    ID    Rank    Action
# Fungi    4751    k    Keep
# Viridiplantae    33090    k    Reject
# This way, I am giving preference to fungal taxonomies and rejecting any plant ones. Note that RANK is LOWERCASE!!
# Also, rank doesn't particularly matter - it just keeps the file naming convention clean. (So for example, even though 
# bacteria is listed as a superkingdom on NCBI, I just put "k" because it's my personal preference and the retrieval still works.)
# ########## FILTER FILE EXAMPLE (-t) ########## 

# <> # TO DO:
# <> # conda qiime1 activation??

# CODE FOLLOWS HERE #

set -e

# Custom error handler function
error_handler() {
    local line_number=$1
    local error_message=$2
    echo "Error on line $line_number: $error_message"
}

# Trap errors and call the error handler
trap 'error_handler ${BASH_LINENO[0]} "$BASH_COMMAND"' ERR

# ARGUMENTS
split_asv_table=false
skip_blast=false

while getopts ":d:o:b:r:e:t:us" opt; do
  case $opt in
    d) DIR="$OPTARG"
    ;;
    o) OUTDIR="$OPTARG"
    ;;    
    b) blast_file="$OPTARG"
    ;;
    r) 
      case "$OPTARG" in
        local|slurm)
          run_type="$OPTARG"
          ;;
        *)
          echo "Invalid value for -r option. Allowed values are 'local' or 'slurm'."
          exit 1
          ;;
      esac
      ;;
    e) EMAIL="$OPTARG"
    ;;
    t) FILTERFILE="$OPTARG"
    ;;    
    u) split_asv_table=true
    ;;
    s) skip_blast=true
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

shift $((OPTIND -1))

# Check for mandatory arguments
if [ -z "$DIR" ] || [ -z "$OUTDIR" ] || [ -z "$blast_file" ] || [ -z "$EMAIL" ] || [ -z "$run_type" ]; then
    echo "Usage: $0 -d <directory_path> -o <desired name of output dir> -b <blast parameter file> -r <local|slurm> -e <email@email.com> [-t <filtertax_file> -m <mismatch number>]"
    exit 1
fi

# Check if email is in correct format
if ! [[ $EMAIL =~ ^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$ ]]; then
    echo "Email is not valid."
    exit 1
fi

# Check if the -t option is not provided
if [ -z "$FILTERFILE" ]; then
    # Ask for confirmation
    read -p "You have not indicated a taxon-to-filter file. This means you are either analyzing a universal amplicon dataset OR you do not wish to assign your ASVs with any taxa preferentially. Is this correct? (yes/no): " choice

    # Process the user's choice
    case "$choice" in
        yes|Yes|YES)
            echo "Proceeding without a taxon-to-filter file."
            ;;
        no|No|NO)
            echo "Exiting. Please provide a taxon-to-filter file using the -t option."
            exit 1
            ;;
        *)
            echo "Invalid choice. Please enter 'yes' or 'no'."
            exit 1
            ;;
    esac
fi

output_dir="${DIR}/${OUTDIR}"
TAXDIR="${DIR}/${OUTDIR}/tax_dir"

echo " - -- --- ---- ---- --- -- -"
echo "Checking for input files"
echo " - -- --- ---- ---- --- -- -"

if [ ! -e "${output_dir}/asvs/rep_set/seqs_chimera_filtered_ASVs.fasta" ]; then
    echo "${output_dir}/asvs/rep_set/seqs_chimera_filtered_ASVs.fasta not found!"
    exit 1
fi

if [ ! -e "${output_dir}/asvs/asv_table_01.biom" ]; then
    echo "${output_dir}/asvs/asv_table_01.biom not found!"
    exit 1
fi

# initiate log
timestamp="$(date +"%Y%m%d_%H:%M:%S")"
output_file="${DIR}/part4_${timestamp}.log"
exec > "$output_file" 2>&1

# log header
echo " - -- --- ---- ---- --- -- -"
echo "Log file for Part 2 of the Microbiome Pipeline. Processing the following arguments:
Working directory: ${DIR}
Output directory: ${OUTDIR}
Blast run file: ${blast_file}
Type of blast: ${run_type}
Email of user: ${EMAIL}
Filterfile if specified: ${FILTERFILE}"
if [ "$skip_blast" = true ]; then
    echo "BLAST was skipped."
fi
if [ "$split_asv_table" = true ]; then
    echo "Final ASV tables will be split into three domains of life since this is universal assay data."
fi
echo " - -- --- ---- ---- --- -- -"

if [ "$skip_blast" = false ]; then
    echo
    echo " - -- --- ---- ---- --- -- -"
    echo "Running BLAST on ASVs"
    echo " - -- --- ---- ---- --- -- -"
    
    cd "${output_dir}"
    
    # Run BLAST script in the background
    blast_iterator.sh "${output_dir}" "${blast_file}" ${run_type} &
    
    # Get the process ID of the last background command
    blast_pid=$!
    
    # Wait for the process to finish
    while ps -p $blast_pid > /dev/null; do
        sleep 1
    done
    
    # Check the exit status of the BLAST script
    blast_exit_status=$?
    if [ $blast_exit_status -eq 0 ]; then
        echo "BLAST successfully completed."
    else
        echo "BLAST failed with exit status $blast_exit_status."
        exit 1
    fi
else
    echo "Skipping BLAST as per user instructions."
fi

echo
echo " - -- --- ---- ---- --- -- -"
echo "Determining Likely Taxonomy of ASVs Using Filters"
echo " - -- --- ---- ---- --- -- -"

# Create filter files in the taxonomy directory if needed. There are two sections - one for if the user didn't specify
# a filter file and another for if they did.
# Read the tab-delimited file, skipping the header

used_taxa=()

add_file_names() {
    local fand="$1"
    local fnot="$2"

    used_taxa+=("$fand" "$fnot")
}

# Process based on FILTERFILE presence
if [ -z "$FILTERFILE" ]; then
    tail -n +2 <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | while IFS=$'\t' read -r Name ID Rank Action; do
        # Generate filenames
        FAND="${TAXDIR}/${Name}__${Rank}_txid${ID}_AND_Environmental_Samples.txt"
        FNOT="${TAXDIR}/${Name}__${Rank}_txid${ID}_NOT_Environmental_Samples.txt"
        # Check the action column and create files accordingly unless they already exist
        if [ "$Action" == "Keep" ]; then
            touch "$FAND"
            touch "$FNOT"
            echo "Creating or locating files ${Name}__${Rank}_txid${ID}_AND_Environmental_Samples.txt and ${Name}__${Rank}_txid${ID}_NOT_Environmental_Samples.txt"
            add_file_names "$FAND" "$FNOT"
        elif [ "$Action" == "Reject" ]; then
            touch "$FNOT"
            echo "Creating or locating file ${Name}__${Rank}_txid${ID}_NOT_Environmental_Samples.txt."
            add_file_names "" "$FNOT"
        fi
    done
else
    tail -n +2 ${FILTERFILE} | while IFS=$'\t' read -r Name ID Rank Action; do
        # Generate filenames
        FAND="${TAXDIR}/${Name}__${Rank}_txid${ID}_AND_Environmental_Samples.txt"
        FNOT="${TAXDIR}/${Name}__${Rank}_txid${ID}_NOT_Environmental_Samples.txt"
        # Check the action column and create files accordingly unless they already exist
        if [ "$Action" == "Keep" ]; then
            touch "$FAND"
            touch "$FNOT"
            echo "Creating or locating files ${Name}__${Rank}_txid${ID}_AND_Environmental_Samples.txt and ${Name}__${Rank}_txid${ID}_NOT_Environmental_Samples.txt"
            add_file_names "$FAND" "$FNOT"
        elif [ "$Action" == "Reject" ]; then
            touch "$FNOT"
            echo "Creating or locating file ${Name}__${Rank}_txid${ID}_NOT_Environmental_Samples.txt."
            add_file_names "" "$FNOT"
        fi
    done
fi

### Update files 
TAXONS=()
TAXIDS=()

for file in "${used_taxa[@]}"; do
    if [[ $file == *NOT* ]]; then
        taxon=$(echo "$file" | perl -ne '@A=split/_NOT_/;@B=split/__\w_/,$A[0];print"$A[0]\n"')
        taxid=$(echo "$file" | perl -ne '@A=split/_NOT_/;@B=split/__\w_/,$A[0];print"$B[1]\n"')
        TAXONS+=("$taxon")
        TAXIDS+=("$taxid")
    fi
done

N=$((${#TAXONS[@]}-1))
echo "$N taxa to update"
# create AND/NOT search-terms and filenames for esearch/efetch
# Function to retrieve taxonomy IDs
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

# Main loop
for ((i = 0; i <= $N; i++)); do
    FAND="${TAXONS[$i]}_AND_Environmental_Samples.txt"
    FNOT="${TAXONS[$i]}_NOT_Environmental_Samples.txt"

    retrieve_taxonomy "${TAXIDS[$i]}[subtree] AND \"Environmental Samples\"[subtree]" "$FAND"
    retrieve_taxonomy "${TAXIDS[$i]}[subtree] NOT \"Environmental Samples\"[subtree]" "$FNOT"
done

# download and unzip updated merged.dmp file
wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip -P ${TAXDIR}
if [ $? -ne 0 ]; then
    echo "Error: Failed to download taxdmp.zip."
    exit 1
fi
unzip -o ${TAXDIR}/taxdmp.zip -d ${TAXDIR}
if [ $? -ne 0 ]; then
    echo "Error: Failed to unzip taxdmp.zip."
    exit 1
fi
rm -rf ${TAXDIR}/taxdmp.zip 

touch "${TAXDIR}/Placeholder__k_txid0_NOT_Environmental_Samples.txt"

# find the user's usearch path, which is where the AccnsWithDubiousTaxAssigns.txt file will be stored.
userpath=$(which usearch)
userdir=$(dirname "$userpath")
tax_files_dir="${userdir}/mbio_taxa_fz"
sudo mkdir -vp "$tax_files_dir" # TODO: will this be problematic, Mario?
sudo touch "${tax_files_dir}/AccnsWithDubiousTaxAssigns.txt" # TODO: will this be problematic, Mario?
sudo chmod 777 "${tax_files_dir}/AccnsWithDubiousTaxAssigns.txt" # TODO: will this be problematic, Mario?
cd ${output_dir}

# This section will create the ASVs2filter.log, which will be used to assign taxonomy. 
# Again, there are two sections - one for if the user didn't specify a filter file and another for if they did.
if [ -z "$FILTERFILE" ]; then
  blast_taxa_categorizer.py \
    -i "${output_dir}/asvs/rep_set/final.blastout" \
    -k $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -e $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_AND_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -r $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -t "$tax_files_dir" \
    -m "${TAXDIR}/merged.dmp"
else
  blast_taxa_categorizer.py \
      -i "${output_dir}/asvs/rep_set/final.blastout" \
      -k $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' "$FILTERFILE" |paste -sd, -) \
      -e $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_AND_Environmental_Samples.txt"}' "$FILTERFILE" |paste -sd, -) \
      -r $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' "$FILTERFILE" |paste -sd, -)  \
      -t "$tax_files_dir" \
      -m "${TAXDIR}/merged.dmp" #-f
fi

bad_accns="${output_dir}/bad_accns.txt"
dubious_accns="${tax_files_dir}/AccnsWithDubiousTaxAssigns.txt"

# Function to check if a number exists in the dubious_accns file
check_existence() {
    local number="$1"
    if grep -q "^$number$" "$dubious_accns"; then
        return 0 # Number exists in the file
    else
        return 1 # Number does not exist in the file
    fi
}

# Filter out lines starting with # or empty lines
filtered_lines=$(grep -vE '^#|^$' "$bad_accns")

while IFS= read -r line; do
    number=$(echo "$line" | cut -f1)
    
    if [ -z "$number" ]; then
        continue
    fi

    if ! check_existence "$number"; then
        echo "$number" >> "$dubious_accns" # Append the number if it doesn't exist
        echo "Adding $number to AccnsWithDubiousTaxAssigns.txt"
        new_addition=true
    fi
done <<< "$filtered_lines"

# If new additions were made, re-run the loop
if [ "$new_addition" = true ]; then
    # Run the blast_taxa_categorizer.py script
    if [ -z "$FILTERFILE" ]; then
        blast_taxa_categorizer.py \
            -i "${output_dir}/asvs/rep_set/final.blastout" \
            -k $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
            -e $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_AND_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
            -r $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
            -t "$tax_files_dir" \
            -m "${TAXDIR}/merged.dmp"
    else
        blast_taxa_categorizer.py \
            -i "${output_dir}/asvs/rep_set/final.blastout" \
            -k $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' "$FILTERFILE" |paste -sd, -) \
            -e $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_AND_Environmental_Samples.txt"}' "$FILTERFILE" |paste -sd, -) \
            -r $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' "$FILTERFILE" |paste -sd, -)  \
            -t "$tax_files_dir" \
            -m "${TAXDIR}/merged.dmp" #-f
    fi
fi

new_addition=false

# Move the generated ASVs files to the rep_set folder
mv "${output_dir}/ASVs2*" "${output_dir}/asvs/rep_set"

mkdir -vp "${output_dir}/asvs/rep_set/assgntax"

module load py-docopt
module load py-biopython
module load py-xmltodict
echo
echo " - -- --- ---- ---- --- -- -"
echo "Assigning Taxonomy With Filters"
echo " - -- --- ---- ---- --- -- -"
rm -rf "${output_dir}/asvs/rep_set/assgntax/taxonomyDB.json"
rm -f *.xml
blast_assign_taxonomy.py -i "${output_dir}/asvs/rep_set/ASVs2filter.log" \
    --db "${output_dir}/asvs/rep_set/assgntax/taxonomyDB.json" --assign_all --add_sizes \
    -m beth.b.peacock@gmail.com \
    -o "${output_dir}/asvs/rep_set/assgntax/seqs_chimera_filtered_tax_assignments.txt"
rm *.xml
module purge
echo
echo " - -- --- ---- ---- --- -- -"
echo "Printing Taxa Levels With Filters"
echo " - -- --- ---- ---- --- -- -"
# count taxa levels
source /sw/miniconda3/bin/activate qiime2
source qiime_shell_helper_functions.sh
count_taxa_levels ${output_dir}/asvs/rep_set/assgntax/seqs_chimera_filtered_tax_assignments.txt > ${output_dir}/asvs/rep_set/assgntax/taxa_levels.txt
echo
echo " - -- --- ---- ---- --- -- -"
echo "Determining Likely Taxonomy of ASVs Without Filters"
echo " - -- --- ---- ---- --- -- -"
blast_taxa_categorizer.py \
    -i "${output_dir}/asvs/rep_set/final.blastout" \
    -k $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -e $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -r $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -t "$tax_files_dir" \
    -m "${TAXDIR}/merged.dmp" #-f

# now we have a new bad_accns file to process:
# Filter out lines starting with # or empty lines
filtered_lines=$(grep -vE '^#|^$' "$bad_accns")

while IFS= read -r line; do
    number=$(echo "$line" | cut -f1)
    
    if [ -z "$number" ]; then
        continue
    fi

    if ! check_existence "$number"; then
        echo "$number" >> "$dubious_accns" # Append the number if it doesn't exist
        echo "Adding $number to AccnsWithDubiousTaxAssigns.txt"
        new_addition=true
    fi
done <<< "$filtered_lines"

# If new additions were made, re-run the loop
if [ "$new_addition" = true ]; then
    # Run the blast_taxa_categorizer.py script
    blast_taxa_categorizer.py \
        -i "${output_dir}/asvs/rep_set/final.blastout" \
        -k $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
        -e $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
        -r $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
        -t "$tax_files_dir" \
        -m "${TAXDIR}/merged.dmp" #-f
fi

# Rename the generated ASVs files and move to the rep_set folder
mv ./ASVs2filter.log ./nf_ASVs2filter.log
mv ./ASVs2summary.txt ./nf_ASVs2summary.txt
mv ./ASVs2reject.txt ./nf_ASVs2reject.txt
mv ./ASVs2keep.txt ./nf_ASVs2keep.txt
mv ./nf_* ${output_dir}/asvs/rep_set

module load py-docopt
module load py-biopython
module load py-xmltodict
echo
echo " - -- --- ---- ---- --- -- -"
echo "Assigning Taxonomy Without Filters"
echo " - -- --- ---- ---- --- -- -"
rm -rf ${output_dir}/asvs/rep_set/assgntax/nf_taxonomyDB.json
rm -f *.xml
blast_assign_taxonomy.py -i ${output_dir}/asvs/rep_set/nf_ASVs2filter.log \
  --db ${output_dir}/asvs/rep_set/assgntax/nf_taxonomyDB.json --assign_all --add_sizes \
  -m beth.b.peacock@gmail.com \
  -o ${output_dir}/asvs/rep_set/assgntax/nf_seqs_chimera_filtered_tax_assignments.txt
rm *.xml
module purge
echo
echo " - -- --- ---- ---- --- -- -"
echo "Printing Taxa Levels Without Filters"
echo " - -- --- ---- ---- --- -- -"
# count taxa levels
source qiime_shell_helper_functions.sh
count_taxa_levels ${output_dir}/asvs/rep_set/assgntax/nf_seqs_chimera_filtered_tax_assignments.txt > ${output_dir}/asvs/rep_set/assgntax/nf_taxa_levels.txt

echo
echo " - -- --- ---- ---- --- -- -"
echo "Adding Taxa and Sequences to ASV Tables"
echo " - -- --- ---- ---- --- -- -"

#add taxa to ASV table
cd ${output_dir}/asvs

OTBL=asv_table_01
biomAddObservations ${OTBL}.biom asv_table_02_add_taxa.biom rep_set/assgntax/seqs_chimera_filtered_tax_assignments.txt

# create three additional taxonomic levels of ASV tables
OTBL=asv_table_02_add_taxa
summarize_taxa.py -i ${OTBL}.biom -L 2,6,7;
to_process=($(find . -maxdepth 1 -type f -name 'asv_table_02_add_taxa*.biom'))

for F in "${to_process[@]}"; do
    FNAME=$(echo "$F" | sed 's|^./asv_table_02_add_taxa||')
    ID=$(echo "$FNAME" | sed 's/.biom//')
    biom2txt $F "asv_table_02_add_taxa${ID}.txt"
    if [[ $split_asv_table ]]; then
        KDOMS=(k__Archaea k__Bacteria k__Eukaryota)
        for K in "${KDOMS[@]}"; do
            grep -P "(#|$K)" "asv_table_02_add_taxa${ID}.txt" > "asv_table_02_add_taxa${ID}.${K}.txt"
            OTBL="asv_table_02_add_taxa${ID}.${K}"
            if ! grep -q "$K" "${OTBL}.txt"; then
                rm "${OTBL}.txt"
            else
                if [[ "${OTBL}.txt" == *"_L"* ]]; then
                    txt2biom_notax "${OTBL}.txt" "${OTBL}.biom"
                    biom_table_math_ops.py -i "${OTBL}.biom" -o "${OTBL}_norm.biom" --normalize2unity
                    biom2txt_notax "${OTBL}_norm.biom" "${OTBL}_norm.txt"
                else 
                    txt2biom "${OTBL}.txt" "${OTBL}.biom"
                    biom_table_math_ops.py -i "${OTBL}.biom" -o "${OTBL}_norm.biom" --normalize2unity
                    biom2txt "${OTBL}_norm.biom" "${OTBL}_norm.txt"
                fi
            fi
        done
    fi
    biom_table_math_ops.py -i ${F} -o "asv_table_02_add_taxa${ID}_norm.biom" --normalize2unity
    biom2txt "asv_table_02_add_taxa${ID}_norm.biom" "asv_table_02_add_taxa${ID}_norm.txt"
done

export MODULEPATH=$MODULEPATH:/sw/spack/share/spack/modules/linux-centos7-cascadelake/
module load r

# add seqs to L8 
otblfp="asv_table_02_add_taxa.txt"
outfp="asv_table_03_add_seqs.txt"

Rscript -e "source('pipeline_helper_functions.R'); add_sequences_to_asv_table('$otblfp', 'rep_set/seqs_chimera_filtered_ASVs.fasta', '$outfp')"

otblfp="asv_table_02_add_taxa_norm.txt"
outfp="asv_table_03_add_seqs_norm.txt"

Rscript -e "source('pipeline_helper_functions.R'); add_sequences_to_asv_table('$otblfp', 'rep_set/seqs_chimera_filtered_ASVs.fasta', '$outfp')"

to_process2=($(find . -maxdepth 1 -type f -name '*taxa.k*txt'))

for F in ${to_process2[@]}; do
    FNAME=$(echo "$F" | sed 's|^./asv_table_02_add_taxa||')
    otblfp=${F}
    outfp="asv_table_03_add_seqs${FNAME}"
    Rscript -e "source('pipeline_helper_functions.R'); add_sequences_to_asv_table('$otblfp', 'rep_set/seqs_chimera_filtered_ASVs.fasta', '$outfp')"
done

echo " - -- --- ---- ---- --- -- -"
echo "Creating Summary File"
echo " - -- --- ---- ---- --- -- -"

module load py-biopython
# get "top 10 contain multiple families" ASVs
python top_ten_family_checker.py final.blastout
# generate final summary file
Rscript -e "source('pipeline_helper_functions.R'); process_data_and_write_excel('asv_table_03_add_seqs_norm.txt', 'rep_set/assgntax/nf_seqs_chimera_filtered_tax_assignments.txt', 'rep_set/assgntax/seqs_chimera_filtered_tax_assignments.txt', 'asv_table_03_add_seqs.txt', 'rep_set/top_ten_family_checker_out.txt')"

echo "All ASV tables have been generated. A summary file can be found here:" | tee /dev/tty
echo $summary_file_name | tee /dev/tty
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