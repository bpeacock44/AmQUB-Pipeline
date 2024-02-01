#!/bin/bash

### USAGE ###
# This script expects to be given at least 4 aguments:
# -d: a working directory, which contains a one folder for each of your fastq files named by ID
# -j: the folders created in the last part that you intend to process in a comma-delimited list (ID1_subset1,ID2,ID3_subset2,etc.)
# -l: the length you want to trim your reads to. Note ALL files will be trimmed to this length.
# -o: the name of your output directory
# -b: the path to your blast script file
# -r: the type of blast run you want to do (local or slurm)

# Optional arguments:
# -t file will have four tab-delimited columns. Include any taxonomic groups that you want to preferentially keep or reject
# as well as their taxonomic ID. ALL taxonomies included under these taxonomic IDs will be treated accordingly so check NCBI
# and make sure. If you are doing a universal assay, do not include the -t flag.
# -m: number of mismatches, if using (again, this should have been specified from part1)

# Examples:
# ./mbio_part3.sh -d /path/to/dir -j "JB141_Nickels01_output,JB143_output" -l 150 -o test1_out -b /path/to/blast.sh -r slurm -t ${MDIR}/filterfile.txt 
# ./mbio_part3.sh -d /path/to/dir -j "JB141_Nickels01_output,JB143_output" -l 150 -o test2_out -b /path/to/blast.sh -r slurm -t ${MDIR}/filterfile.txt -m 1 
# ./mbio_part3.sh -d /path/to/dir -j "JB141_Nickels01_output,JB143_output" -l 150 -o test3_out -b /path/to/blast.sh -r local
# ./mbio_part3.sh -d /path/to/dir -j "JB141_Nickels01_output,JB143_output" -l 150 -o test4_out -b /path/to/blast.sh -r local -m 1 

### INPUT ###
# This script follows part 2, which must be completed first. You will look at your trim stats, 
# determine what length you want to trim to, and run this code to finish the analysis.

# So, as an example, your working directory might now include:
#       Folder JB141_Nickels01_output and directory JB143_output, both containing output of part2.

# When this code is run, a new directory named as you indicated will be created for analysis output. 



# ########## SLURM BLAST FILE EXAMPLE ########## 
##!/bin/bash 
##SBATCH -p i128
##SBATCH -c 128
## any other parameters or modules needed
#module load blast-plus
#
##<>#<>#<>#<>#<>
## YOU MUST SET THESE:
##<>#<>#<>#<>#<>
#DATABASE_PATH=/sw/dbs/blast_db_download/nt
#NUMTHREADS=128
#
##<>#<>#<>#<>#<>
## GENERALLY DON'T CHANGE THESE:
##<>#<>#<>#<>#<>
#OPTS="qseqid sseqid pident length mismatch evalue bitscore staxids stitle"
#TASK=blastn
#INFASTA=$1
#MAXTSEQS=$2  
#EVAL=0.001
#echo blastn -task $TASK -db $DATABASE_PATH -query $INFASTA -max_target_seqs $MAXTSEQS -evalue $EVAL -num_threads $NUMTHREADS -outfmt "7 $OPTS" 
# ########## SLURM BLAST FILE EXAMPLE ########## 



# If running a local blast, use the same format but don't include the SBATCH lines.



# ########## FILTER FILE EXAMPLE ########## 
# If I am doing fungal ITS taken from a plant sample, then the file might include:
# Name    ID    Rank    Action
# Fungi    4751    k    Keep
# Viridiplantae    33090    k    Reject
# This way, I am giving preference to fungal taxonomies and rejecting any plant ones. Note that RANK is LOWERCASE!!
# ########## FILTER FILE EXAMPLE ########## 

# <> # TO DO:
# <> # HDIR FILES SHOULD BE IN PATH
# <> # conda qiime1 activation??

# CODE FOLLOWS HERE #

while getopts ":d:j:l:o:b:r:t:m:" opt; do
  case $opt in
    d) DIR="$OPTARG"
    ;;
    j) IFS=',' read -ra JBS <<< "$OPTARG"
    ;;
    l) LEN="$OPTARG"
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
    t) FILTERFILE="$OPTARG"
    ;;
    m) mmatchnum="$OPTARG"
    ;;    
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

# Check for mandatory arguments
if [ -z "$DIR" ] || [ -z "${JBS[*]}" ] || [ -z "$LEN" ] || [ -z "$OUTDIR" ] || [ -z "$blast_file" ] || [ -z "$run_type" ]; then
    echo "Usage: $0 -d <directory_path> -j <comma-separated output directories> -l <trim length> -o <desired name of output dir> -b <blast parameter file> -r <local|slurm> [-t <filtertax_file> -m <mismatch number>]"
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

# Remove "_output" from each element in the array
for ((i=0; i<${#JBS[@]}; i++)); do
  JBS[$i]=${JBS[$i]%%_output}
done

echo " - -- --- ---- ---- --- -- -"
echo "Checking for input files"
echo " - -- --- ---- ---- --- -- -"
HDIR=/sw/paul_helper_scripts

# show your fastq files and map
for JB in "${JBS[@]}"; do
    if [ ! -e "${DIR}/${JB}_output/${JB}_A1P1.M${mmatchnum}.fq" ]; then
        echo "File ${DIR}/${JB}_output/${JB}_A1P1.M${mmatchnum}.fq not found!"
        exit 1
    fi

    if [ ! -e "${DIR}/${JB}_output/${JB}_A1P2.M${mmatchnum}.fq" ]; then
        echo "File ${DIR}/${JB}_output/${JB}_A1P2.M${mmatchnum}.fq not found!"
        exit 1
    fi
done

# initiate log
timestamp="$(date +"%Y%m%d_%H:%M:%S")"
output_file="${DIR}/part3_${timestamp}.log"
exec > "$output_file" 2>&1


# log header
echo " - -- --- ---- ---- --- -- -"
echo "Log file for Part 2 of the Microbiome Pipeline. Processing the following arguments:
Working Directory: ${DIR}
Data IDs: ${JBS[@]}
Trim Length: ${LEN}
Output Directory: ${OUTDIR}
Filterfile if specified: ${FILTERFILE}
Mismatches if specified: ${mmatchnum}
 - -- --- ---- ---- --- -- -"

echo
echo " - -- --- ---- ---- --- -- -"
echo "Trimming and Filtering Reads"
echo " - -- --- ---- ---- --- -- -"

# set up and go to output directory
output_dir="${DIR}/${OUTDIR}"
mkdir -vp ${output_dir}
cd "${output_dir}"

echo "This folder contains the results of ${JBS[@]} trimmed at ${LEN}." > "${output_dir}/summary.txt"

# make tax directory
mkdir -vp "${DIR}/${OUTDIR}/tax_dir"
TAXDIR="${DIR}/${OUTDIR}/tax_dir"

#truncate reads at LEN
for JB in ${JBS[@]}; do
  usearch -fastx_truncate "${DIR}/${JB}_output/${JB}_A1P1.M${mmatchnum}.fq" -quiet -trunclen ${LEN} -fastqout "${DIR}/${JB}_output/${JB}_A1P1_${LEN}bp.fq"; # &
done

rm -f "${output_dir}/combined.fq"
if [ ${#JBS[@]} -gt 1 ]; then
  #EITHER combine if you have multiple files
  for JB in ${JBS[@]}; do
  cat "${DIR}/${JB}_output/${JB}_A1P1_${LEN}bp.fq" >> "${output_dir}/combined.fq"
  done
else
  #OR create a symlink if you have only one file
  ln -s "${DIR}/${JB}_output/${JB}_A1P1_${LEN}bp.fq" "${output_dir}/combined.fq"
fi

#maxee quality filtering of demultiplexed/truncated fq files (*** keep THREADS=1 for repeatability ***)
for JB in ${JBS[@]}; do
  usearch -threads 1 -fastq_filter "${DIR}/${JB}_output/${JB}_A1P1_${LEN}bp.fq" -quiet -fastq_maxee 1.0 -fastaout "${DIR}/${JB}_output/${JB}.filtered.fa"
done
echo
echo " - -- --- ---- ---- --- -- -"
echo "Pooling Samples and Creating OTUs"
echo " - -- --- ---- ---- --- -- -"

#sample pooling (https://www.drive5.com/usearch/manual/pool_samples.html)
rm -f "${output_dir}/filtered.fa"
if [ ${#JBS[@]} -gt 1 ]; then
  #EITHER combine if you have multiple files
  for JB in ${JBS[@]}; do
  cat "${DIR}/${JB}_output/${JB}_A1P1_${LEN}bp.fq" >> "${output_dir}/filtered.fa"
  done
else
  #OR create a symlink if you have only one file
  ln -s "${DIR}/${JB}_output/${JB}.filtered.fa" "${output_dir}/filtered.fa"
fi

#"The data requirements of cluster_otus and unoise3 are the same: the input should be a set of unique sequences which have
# been trimmed and quality filtered. It is therefore easy to generate both 97% OTUs and ZOTUs in the same pipeline, and
# I therefore suggest building an OTU table using both strategies." - Robert Edgar
#
# We decided to use only ZOTUS with unoise3
#find unique sequences
usearch -fastx_uniques "${output_dir}/filtered.fa" -quiet -fastaout "${output_dir}/uniques.fa" -sizeout -relabel Uniq

#make a subdirectory for the zotus
mkdir -vp "${output_dir}/zotus"

#cluster unique sequences into ZOTUS using the UNOISE3 algorithm (default -minsize=8)
usearch -unoise3 "${output_dir}/uniques.fa" -quiet -zotus "${output_dir}/zotus/zotus.fa"

# Convert '>Zotu' to '>Otu' in the file
sed 's/>Zotu/>Otu/g' "${output_dir}/zotus/zotus.fa" > "${output_dir}/zotus/z.fa"

# Check if the replacement was successful before overwriting
if grep -q '>Otu' "${output_dir}/zotus/z.fa"; then
    # Overwrite the original file if '>Otu' is found
    mv -v "${output_dir}/zotus/z.fa" "${output_dir}/zotus/zotus.fa"
else
    echo "Header replacement failed. Original file not overwritten."
fi
echo
echo " - -- --- ---- ---- --- -- -"
echo "Creating Initial OTU Table"
echo " - -- --- ---- ---- --- -- -"

#create an OTU table ("Input should be reads before quality filtering and before discarding low-abundance unique sequences, e.g. singletons")
usearch --otutab "${output_dir}/combined.fq" -quiet -otus "${output_dir}/zotus/zotus.fa" -otutabout "${output_dir}/zotus/otu_table_00.txt"

cd "${output_dir}"
#use R to sort OTU table
export MODULEPATH=$MODULEPATH:/sw/spack/share/spack/modules/linux-centos7-cascadelake/
module load r
Rscript "${HDIR}/processing_otu.R"

source /sw/miniconda3/bin/activate qiime1

#source for bash helper functions
source "${HDIR}/qiime_shell_helper_functions.sh"

#convert to biom
OTBL=otu_table_01
txt2biom_notax "${output_dir}/zotus/${OTBL}.txt" "${output_dir}/zotus/${OTBL}.biom"
biomsummarize "${output_dir}/zotus/${OTBL}.biom" "${output_dir}/zotus/${OTBL}.biom.sum"

cd "${output_dir}/zotus"
#add counts to otu table
Rscript "${HDIR}/add_counts_to_fasta_seqs.R" 
cd "${output_dir}"

echo
echo " - -- --- ---- ---- --- -- -"
echo "Running BLAST on OTUs"
echo " - -- --- ---- ---- --- -- -"

mkdir -vp "${output_dir}/zotus/rep_set"
mv -v "${output_dir}/zotus/seqs_chimera_filtered_otus.fasta" "${output_dir}/zotus/rep_set"
cd "${output_dir}"

# Run BLAST script in the background
${DIR}/micro_blast_v2.sh "${output_dir}" "${blast_file}" ${run_type} &

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

# This section will create the otus2filter.log, which will be used to assign taxonomy. 
# Again, there are two sections - one for if the user didn't specify a filter file and another for if they did.
if [ -z "$FILTERFILE" ]; then
  "${HDIR}/filter_contaminating_reads.py" \
    -i "${output_dir}/zotus/rep_set/5000.rb0.blastout" \
    -k $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -e $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_AND_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -r $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nEukaryota\t2759\tk\tKeep\nBacteria\t2\tk\tKeep\nArchaea\t2157\tk\tKeep\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -m "${TAXDIR}/merged.dmp"
else
  "${HDIR}/filter_contaminating_reads.py" \
      -i "${output_dir}/zotus/rep_set/5000.rb0.blastout" \
      -k $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' "$FILTERFILE" |paste -sd, -) \
      -e $(awk -F'\t' '$4=="Keep"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_AND_Environmental_Samples.txt"}' "$FILTERFILE" |paste -sd, -) \
      -r $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' "$FILTERFILE" |paste -sd, -)  \
      -m "${TAXDIR}/merged.dmp" #-f
fi

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
    -i "${output_dir}/zotus/rep_set/5000.rb0.blastout" \
    -k $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -e $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -r $(awk -F'\t' '$4=="Reject"{print "'${TAXDIR}'/"$1"__"$3"_txid"$2"_NOT_Environmental_Samples.txt"}' <(echo -e "Name\tID\tRank\tAction\nPlaceholder\t0\tk\tReject") | paste -sd, -) \
    -m "${TAXDIR}/merged.dmp" #-f

# Rename the generated otus files and move to the rep_set folder
mv ./otus2filter.log ./nf_otus2filter.log
mv ./otus2summary.txt ./nf_otus2summary.txt
mv ./otus2reject.txt ./nf_otus2reject.txtqq
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

module load r

for F in otu_table_02_add_taxa*.biom; do
  FNAME=$(echo "$F" | sed 's/otu_table_02_add_taxa//')
  ID=$(echo "$FNAME" | sed 's/.biom//')
  $HDIR/biom_table_math_ops.py -i ${F} -o "otu_table_02_add_taxa_norm${FNAME}" --normalize2unity
  biom2txt $F "otu_table_02_add_taxa${ID}.txt"
  biom2txt "otu_table_02_add_taxa_norm${FNAME}" "otu_table_02_add_taxa_norm${ID}.txt"
done

# add seqs to L8 (regular)
Rscript "${HDIR}/add_seqs_to_OTU.R" "otu_table_02_add_taxa.txt" "otu_table_03_add_seqs.txt"
Rscript "${HDIR}/add_seqs_to_OTU.R" "otu_table_02_add_taxa_norm.txt" "otu_table_03_add_seqs_norm.txt"

echo " - -- --- ---- ---- --- -- -"
echo "Creating Summary File"
echo " - -- --- ---- ---- --- -- -"




