#!/bin/bash

### USAGE ###
# This script expects to be given at least 4 aguments:
# -d: a working directory, which contains one folder for each of your fastq files named by ID
# -o: the name of your output directory
# -b: the path to your blast script file (should include ALL jobs, each beginning with a shebang as below)
# -n: number of different jobs you want to run
# -r: the type of blast run you want to do (local or slurm)

# Examples:
# mbio_just_blast.sh -d /path/to/dir -o test1_out -b /path/to/blast.sh -n 2 -r slurm

### INPUT ###
# This script follows part 3, which must be completed first. 
# The output file will have already been generated in part 3.



## NOTE THAT YOU CANNOT SUBMIT BATCHES FROM WITHIN THE SINGULARITY. If you want to run blast across various computational 
## resources, it is better to do that manually on your system and then resume this pipeline after. 
## See tutorials on github for guidance.



# ########## SLURM BLAST FILE EXAMPLE (-b), if you had 2 jobs to run ########## 
##!/bin/bash 
##SBATCH -p epyc
##SBATCH -c 256
##SBATCH --mem=300G 
## any other parameters or modules needed
#module load ncbi-blast/2.6.0+
#module load db-ncbi

##<>#<>#<>#<>#<>
## YOU MUST SET THESE:
##<>#<>#<>#<>#<>
#DATABASE_PATH=$NCBI_DB/nt
#NUMTHREADS=256

##<>#<>#<>#<>#<>
## GENERALLY DON'T CHANGE THESE:
##<>#<>#<>#<>#<>
#OPTS="qseqid sseqid pident length mismatch evalue bitscore staxids stitle qcovs"
#TASK=blastn
#INFASTA=$1
#MAXTSEQS=$2  
#EVAL=0.001
# blastn -task $TASK -db $DATABASE_PATH -query $INFASTA -max_target_seqs $MAXTSEQS -evalue $EVAL -num_threads $NUMTHREADS -outfmt "7 $OPTS" 

##!/bin/bash 
##SBATCH -c 256
##SBATCH --mem=300G 
## any other parameters or modules needed
#module load ncbi-blast/2.6.0+
#module load db-ncbi

##<>#<>#<>#<>#<>
## YOU MUST SET THESE:
##<>#<>#<>#<>#<>
#DATABASE_PATH=$NCBI_DB/nt
#NUMTHREADS=256

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

# CODE FOLLOWS HERE #

set -e

# Custom error handler function
error_handler() {
    local error_message=$1
    echo "Error on line $error_message" | tee /dev/tty
}

# Trap errors and call the error handler
trap 'error_handler "$BASH_COMMAND"' ERR

while getopts ":d:o:b:n:r:" opt; do
  case $opt in
    d) DIR="$OPTARG"
    ;;
    o) OUTDIR="$OPTARG"
    ;;    
    b) blast_file="$OPTARG"
    ;;
    n) num="$OPTARG"
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
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

shift $((OPTIND -1))

# Check for mandatory arguments
if [ -z "$DIR" ] || [ -z "$OUTDIR" ] || [ -z "$blast_file" ] || [ -z "$num" ] || [ -z "$run_type" ]; then
    echo "Usage: $0 -d <directory_path> -o <desired name of output dir> -b <blast parameter file> -n <numbers of blasts> -r <local|slurm>"
    exit 1
fi

output_dir="${DIR}/${OUTDIR}"

echo " - -- --- ---- ---- --- -- -"
echo "Checking for input files"
echo " - -- --- ---- ---- --- -- -"

if [ ! -e "${output_dir}/asvs/rep_set/seqs_chimera_filtered_ASVs.fasta" ]; then
    echo "${output_dir}/asvs/rep_set/seqs_chimera_filtered_ASVs.fasta not found!"
    exit 1
fi

# initiate log
timestamp="$(date +"%Y%m%d_%H:%M:%S")"
output_file="${DIR}/part4_${timestamp}.log"
exec > "$output_file" 2>&1

# log header
echo " - -- --- ---- ---- --- -- -"
echo "Log file for Part 4 of the Microbiome Pipeline. Processing the following arguments:
Working directory: ${DIR}
Output directory: ${OUTDIR}
Blast run file: ${blast_file}
Type of blast: ${run_type}
Number of blast jobs: ${num}"
echo " - -- --- ---- ---- --- -- -"

echo
echo " - -- --- ---- ---- --- -- -"
echo "Running BLAST on ASVs"
echo " - -- --- ---- ---- --- -- -"
    
cd "${output_dir}"

# split ASV file into multiple files, with prefix
fasta_file="${output_dir}/asvs/rep_set/seqs_chimera_filtered_ASVs.fasta"  # Specify your input FASTA file
seqkit split "$fasta_file" -p ${num} -O "${output_dir}/asvs/rep_set"

filename=seqs_chimera_filtered_ASVs.fasta
# Rename files
files=(${output_dir}/asvs/rep_set/seqs_chimera_filtered_ASVs.part*)
for ((i = 0; i < ${#files[@]}; i++)); do
    new_file_name="${output_dir}/asvs/rep_set/$((i + 1))_${filename}"
    mv "${files[i]}" "$new_file_name"
    echo "Created $new_file_name"
done

# split blast_file into multiple files, with prefix
# Count the number of shebangs in the file
shebang_count=$(grep -c '^#!/' "${blast_file}")

# Check if the number of shebangs matches the expected number
if [ "$shebang_count" -ne "$num" ]; then
    echo "Error: The file does not contain $num scripts."
    exit 1
fi

# Split the file
awk 'BEGIN {file_number=1} /^#!/ { if (file_number <= '"$num"') { close((file_number-1) "_blast.sh"); file_number++; } } { print > ((file_number-1) "_blast.sh") }' "$blast_file"

# Run BLAST script in the background
for N in $(seq 1 "$num"); do
	multi_blast_iterator.sh "${output_dir}" "${run_type}" "${N}" &
done

wait

echo "Merging all blastout files."

rm -f "${output_dir}/asvs/rep_set/final.blastout"
cat "${output_dir}/asvs/rep_set/"*".blastout" | grep -v "# BLAST processed" >> "${output_dir}/asvs/rep_set/final.blastout"

# Extract headers from seqs_chimera_filtered_ASVs.fasta
grep '^>' "${output_dir}/asvs/rep_set/seqs_chimera_filtered_ASVs.fasta" | sed 's/^>//g' > "${output_dir}/asvs/rep_set/headers.txt"

# Check if each header is present in final.blastout
while IFS= read -r header; do
    if ! grep -q "${header}" "${output_dir}/asvs/rep_set/final.blastout"; then
        echo "Header '$header' not found in final.blastout. There was an issue with the blast."
        #exit 1
    fi
done < "${output_dir}/asvs/rep_set/headers.txt"

echo "# BLAST processed" >> "${output_dir}/asvs/rep_set/final.blastout"


