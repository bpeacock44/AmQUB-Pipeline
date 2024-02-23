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
# and make sure. If you are doing a universal assay, do not include the -t flag.
# -m: number of mismatches, if using (again, this should have been specified from part1)

# Examples:
# ./mbio_part4.sh -d /path/to/dir -o test1_out -b /path/to/blast.sh -e email@email.com -r slurm -t ${MDIR}/filterfile.txt 
# ./mbio_part4.sh -d /path/to/dir -o test2_out -b /path/to/blast.sh -e email@email.com -r slurm -t ${MDIR}/filterfile.txt -m 1 
# ./mbio_part4.sh -d /path/to/dir -o test3_out -b /path/to/blast.sh -e email@email.com -r local
# ./mbio_part4.sh -d /path/to/dir -o test4_out -b /path/to/blast.sh -e email@email.com -r local -m 1 

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
# ########## FILTER FILE EXAMPLE (-t) ########## 

# <> # TO DO:
# <> # HDIR FILES SHOULD BE IN PATH
# <> # conda qiime1 activation??

# CODE FOLLOWS HERE #

while getopts ":d:o:b:r:e:t:m:" opt; do
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
    m) mmatchnum="$OPTARG"
    ;;    
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

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

# Check if mmatchnum is not defined or is set to an empty string, set it to "0"
if [ -z "${mmatchnum+x}" ] || [ -z "$mmatchnum" ]; then
    mmatchnum="0"
fi

# Check if the -t option is not provided
if [ -z "$FILTERFILE" ]; then
    # Ask for confirmation
    read -p "You have not indicated a taxon-to-filter file. This means you are either analyzing a universal amplicon dataset OR you do not wish to assign your OTUs with any taxa preferentially. Is this correct? (yes/no): " choice

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
HDIR=/sw/paul_helper_scripts

if [ ! -e "${output_dir}/zotus/rep_set/seqs_chimera_filtered_otus.fasta" ]; then
    echo "${output_dir}/zotus/rep_set/seqs_chimera_filtered_otus.fasta not found!"
    exit 1
fi

if [ ! -e "${output_dir}/zotus/otu_table_01.biom" ]; then
    echo "${output_dir}/zotus/otu_table_01.biom not found!"
    exit 1
fi

# initiate log
timestamp="$(date +"%Y%m%d_%H:%M:%S")"
output_file="${DIR}/part3_${timestamp}.log"
exec > "$output_file" 2>&1

# log header
echo " - -- --- ---- ---- --- -- -"
echo "Log file for Part 2 of the Microbiome Pipeline. Processing the following arguments:
Working Directory: ${DIR}
Output Directory: ${OUTDIR}
Blast Run File: ${blast_file}
Type of Blast: ${run_type}
Email of user: ${EMAIL}
Filterfile if specified: ${FILTERFILE}
Mismatches if specified: ${mmatchnum}
 - -- --- ---- ---- --- -- -"
echo
echo " - -- --- ---- ---- --- -- -"
echo "Running BLAST on OTUs"
echo " - -- --- ---- ---- --- -- -"

cd "${output_dir}"

# Run BLAST script in the background
${HDIR}/micro_blast_script.sh "${output_dir}" "${blast_file}" ${run_type} &

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
echo
echo " - -- --- ---- ---- --- -- -"
echo "Determining Likely Taxonomy of OTUs Using Filters"
echo " - -- --- ---- ---- --- -- -"

# Create filter files in the taxonomy directory if needed. There are two sections - one for if the user didn't specify
# a filter file and another for if they did.
# Read the tab-delimited file, skipping the header
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
      elif [ "$Action" == "Reject" ]; then
          touch "$FNOT"
          echo "Creating or locating file ${Name}__${Rank}_txid${ID}_NOT_Environmental_Samples.txt."
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
      elif [ "$Action" == "Reject" ]; then
          touch "$FNOT"
          echo "Creating or locating file ${Name}__${Rank}_txid${ID}_NOT_Environmental_Samples.txt."
      fi
  done
fi

### Update files - note this will update ALL taxonomy files, not just the ones you are using.
TAXONS=($(ls ${TAXDIR}/*NOT* | perl -ne '@A=split/_NOT_/;@B=split/__\w_/,$A[0];print"$A[0]\n";'))
TAXIDS=($(ls ${TAXDIR}/*NOT* | perl -ne '@A=split/_NOT_/;@B=split/__\w_/,$A[0];print"$B[1]\n";'))
N=$((${#TAXONS[@]}-1)); echo "$N taxa to update"

#create AND/NOT search-terms and filenames for esearch/efetch
for I in $(seq 0 $N); do
  FAND=${TAXONS[$I]}_AND_Environmental_Samples.txt
  FNOT=${TAXONS[$I]}_NOT_Environmental_Samples.txt
  echo "Retrieving taxonomy for ${FAND}" 
  esearch -db taxonomy -query "${TAXIDS[$I]}[subtree] AND \"Environmental Samples\"[subtree]" | efetch -format uid > $FAND
  echo "Retrieving taxonomy for ${FNOT}" 
  esearch -db taxonomy -query "${TAXIDS[$I]}[subtree] NOT \"Environmental Samples\"[subtree]" | efetch -format uid > $FNOT
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

# This section will create the otus2filter.log, which will be used to assign taxonomy. 
# Again, there are two sections - one for if the user didn't specify a filter file and another for if they did.
if [ -z "$FILTERFILE" ]; then
  "${HDIR}/filter_contaminating_reads.py" \
    -i "${output_dir}/zotus/rep_set/final.blastout" \
    -k $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -e $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_AND_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -r $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -t "$tax_files_dir" \
    -m "${TAXDIR}/merged.dmp"
else
  "${HDIR}/filter_contaminating_reads.py" \
      -i "${output_dir}/zotus/rep_set/final.blastout" \
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
    # Run the filter_contaminating_reads.py script
    if [ -z "$FILTERFILE" ]; then
        "${HDIR}/filter_contaminating_reads.py" \
            -i "${output_dir}/zotus/rep_set/final.blastout" \
            -k $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
            -e $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_AND_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
            -r $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
            -t "$tax_files_dir" \
            -m "${TAXDIR}/merged.dmp"
    else
        "${HDIR}/filter_contaminating_reads.py" \
            -i "${output_dir}/zotus/rep_set/final.blastout" \
            -k $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' "$FILTERFILE" |paste -sd, -) \
            -e $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_AND_Environmental_Samples.txt"}' "$FILTERFILE" |paste -sd, -) \
            -r $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' "$FILTERFILE" |paste -sd, -)  \
            -t "$tax_files_dir" \
            -m "${TAXDIR}/merged.dmp" #-f
    fi
fi

new_addition=false

# Move the generated otus files to the rep_set folder
mv ./otus2* "${output_dir}/zotus/rep_set"

mkdir -vp "${output_dir}/zotus/rep_set/assgntax"

module load py-docopt
module load py-biopython
module load py-xmltodict
echo
echo " - -- --- ---- ---- --- -- -"
echo "Assigning Taxonomy With Filters"
echo " - -- --- ---- ---- --- -- -"
rm -rf "${output_dir}/zotus/rep_set/assgntax/taxonomyDB.json"
rm -f *.xml
${HDIR}/blast_assign_taxonomy.py -i "${output_dir}/zotus/rep_set/otus2filter.log" \
    --db "${output_dir}/zotus/rep_set/assgntax/taxonomyDB.json" --assign_all --add_sizes \
    -m beth.b.peacock@gmail.com \
    -o "${output_dir}/zotus/rep_set/assgntax/seqs_chimera_filtered_tax_assignments.txt"
rm *.xml
module purge
echo
echo " - -- --- ---- ---- --- -- -"
echo "Printing Taxa Levels With Filters"
echo " - -- --- ---- ---- --- -- -"
# count taxa levels
source /sw/miniconda3/bin/activate qiime2
source ${HDIR}/qiime_shell_helper_functions.sh
count_taxa_levels ${output_dir}/zotus/rep_set/assgntax/seqs_chimera_filtered_tax_assignments.txt > ${output_dir}/zotus/rep_set/assgntax/taxa_levels.txt
echo
echo " - -- --- ---- ---- --- -- -"
echo "Determining Likely Taxonomy of OTUs Without Filters"
echo " - -- --- ---- ---- --- -- -"
"${HDIR}/filter_contaminating_reads.py" \
    -i "${output_dir}/zotus/rep_set/final.blastout" \
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
    # Run the filter_contaminating_reads.py script
    "${HDIR}/filter_contaminating_reads.py" \
        -i "${output_dir}/zotus/rep_set/final.blastout" \
        -k $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
        -e $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
        -r $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
        -t "$tax_files_dir" \
        -m "${TAXDIR}/merged.dmp" #-f
fi

# Rename the generated otus files and move to the rep_set folder
mv ./otus2filter.log ./nf_otus2filter.log
mv ./otus2summary.txt ./nf_otus2summary.txt
mv ./otus2reject.txt ./nf_otus2reject.txt
mv ./otus2keep.txt ./nf_otus2keep.txt
mv ./nf_* ${output_dir}/zotus/rep_set

module load py-docopt
module load py-biopython
module load py-xmltodict
echo
echo " - -- --- ---- ---- --- -- -"
echo "Assigning Taxonomy Without Filters"
echo " - -- --- ---- ---- --- -- -"
rm -rf ${output_dir}/zotus/rep_set/assgntax/nf_taxonomyDB.json
rm -f *.xml
${HDIR}/blast_assign_taxonomy.py -i ${output_dir}/zotus/rep_set/nf_otus2filter.log \
  --db ${output_dir}/zotus/rep_set/assgntax/nf_taxonomyDB.json --assign_all --add_sizes \
  -m beth.b.peacock@gmail.com \
  -o ${output_dir}/zotus/rep_set/assgntax/nf_seqs_chimera_filtered_tax_assignments.txt
rm *.xml
module purge
echo
echo " - -- --- ---- ---- --- -- -"
echo "Printing Taxa Levels Without Filters"
echo " - -- --- ---- ---- --- -- -"
# count taxa levels
source ${HDIR}/qiime_shell_helper_functions.sh
count_taxa_levels ${output_dir}/zotus/rep_set/assgntax/nf_seqs_chimera_filtered_tax_assignments.txt > ${output_dir}/zotus/rep_set/assgntax/nf_taxa_levels.txt

echo
echo " - -- --- ---- ---- --- -- -"
echo "Adding Taxa and Sequences to OTU Tables"
echo " - -- --- ---- ---- --- -- -"

#add taxa to otu table
cd ${output_dir}/zotus

OTBL=otu_table_01
biomAddObservations ${OTBL}.biom otu_table_02_add_taxa.biom rep_set/assgntax/seqs_chimera_filtered_tax_assignments.txt

# create three additional taxonomic levels of OTU tables
OTBL=otu_table_02_add_taxa
summarize_taxa.py -i ${OTBL}.biom -L 2,6,7;

for F in otu_table_02_add_taxa*.biom; do
  FNAME=$(echo "$F" | sed 's/otu_table_02_add_taxa//')
  ID=$(echo "$FNAME" | sed 's/.biom//')
  $HDIR/biom_table_math_ops.py -i ${F} -o "otu_table_02_add_taxa_norm${FNAME}" --normalize2unity
  biom2txt $F "otu_table_02_add_taxa${ID}.txt"
  biom2txt "otu_table_02_add_taxa_norm${FNAME}" "otu_table_02_add_taxa_norm${ID}.txt"
done

export MODULEPATH=$MODULEPATH:/sw/spack/share/spack/modules/linux-centos7-cascadelake/
module load r

# add seqs to L8 (regular)
Rscript "${HDIR}/add_seqs_to_OTU.R" "otu_table_02_add_taxa.txt" "otu_table_03_add_seqs.txt"
Rscript "${HDIR}/add_seqs_to_OTU.R" "otu_table_02_add_taxa_norm.txt" "otu_table_03_add_seqs_norm.txt"

echo " - -- --- ---- ---- --- -- -"
echo "Creating Summary File"
echo " - -- --- ---- ---- --- -- -"

module load py-biopython
# get "top 10 contain multiple families" OTUs
python ${HDIR}/top_10_family_checker.py final.blastout
# generate final summary file
Rscript ${HDIR}/final_summary_table_maker.R

echo "All OTU tables have been generated. A summary file can be found here:" | tee /dev/tty
echo $summary_file_name | tee /dev/tty
echo " - -- --- ---- ---- --- -- -"  | tee /dev/tty
echo "Final Recommendations"  | tee /dev/tty
echo " - -- --- ---- ---- --- -- -"  | tee /dev/tty
echo "Sometimes certain OTUs can be primarily associated with your controls (likely a source of contamination,
which may be from other samples or even from the individual building the library in the first place). Make sure 
to check your OTU tables for these OTUs - just look at your control columns and see if any of the OTUs are 
relatively high in them and not present in the regular samples. These OTUs should be removed before further 
analyses. " | tee /dev/tty