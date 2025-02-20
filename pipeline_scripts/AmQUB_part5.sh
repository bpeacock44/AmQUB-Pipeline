#!/bin/bash
# See mbio_tutorial.md for further guidance!

### USAGE ###
# This script accepts two types of input:
# 1. Direct command-line arguments:
#    AmQUB_part5.sh -i output_dir -k to_keep.txt
#    AmQUB_part5.sh -i output_dir -k to_keep.txt -u -c -d -m  

## Required Flags
# -i: The output directory generated in part 3 that you also ran part 4 on.
# -k: A file of taxonomic units you wish to keep in the final table - one per line.

## Optional Flags:
# -u: Universal assay - causes final ASV tables to be split into taxonomic groups prior to normalizing
# -c: Optionally choose to use the classifier taxonomic assignments instead of BLAST.
# -d: Creates a specialized "Detailed Informational OTU/ASV Table" with more detail about the taxonomic assignments.
# -m: Will run an optional analysis checking whether the top 10 BLAST results for each taxonomic unit 
	# had more than one family present. This can be an indicator that the BLAST results should be examined
	# more closely. Note that runtime will increase as this is a more intensive analysis. 

# 2. A parameter template file:
#    AmQUB_part4.sh params.csv
#    Where params.csv contains the following rows, comma delimited, no white space between.
#    The labels at the beginning of each row should be the same as below.
#    The following shows what all possible rows:

#        Part 3 Output Folder To Process,output_dir
#        Taxonomic Units to Keep File,to_keep.txt
#        Universal Assay,true
#        Classifier Assignments Primary,true
#        Detailed Informational OTU/ASV Table,true
#        Mixed Top 10 Analysis,true
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
            "Taxonomic Units to Keep File,"*) KEEP="${line#*,}" ;;
            "Universal Assay,"*) UNI="${line#*,}" ;;
            "Classifier Assignments Primary,"*) CLA="${line#*,}" ;;
            "Detailed Informational OTU/ASV Table,"*) DET="${line#*,}" ;;
            "Mixed Top 10 Analysis,"*) MIX="${line#*,}" ;;
        esac
    done < "$param_file"
}

# Initialize default values
UNI=false
CLA=false
DET=false
MIX=false

    echo " - -- --- ---- ---- --- -- -
 █████╗             ██████╗ ██╗   ██╗██████╗ 
██╔══██╗████╗ ████╗██║   ██╗██║   ██║██╔══██╗
███████║██╔████╔██║██║   ██║██║   ██║██████╔╝
██╔══██║██║╚██╔╝██║██║   ██║██║   ██║██╔══██╗
██║  ██║██║ ╚═╝ ██║╚██████╔╝╚██████╔╝██████╔╝
╚═╝  ╚═╝╚═╝     ╚═╝ ╚════██╗ ╚═════╝ ╚═════╝  
PART 5:Final Tables      ╚═╝
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
            -i|--output) output_dir="$2"; shift 2 ;;
            -k|--keep_file) KEEP="$2"; shift 2 ;;
            -u|--universal) UNI=true; shift ;;
            -c|--classifier) CLA=true; shift ;;
            -d|--detailed_summary) DET=true; shift ;;
            -m|--mixed_fam_analysis) MIX=true; shift ;;
            *) echo "Unknown option: $1" >&2; exit 1 ;;
        esac
    done
fi

# Check for mandatory arguments
if [ -z "$output_dir" ] || [ -z "$KEEP" ]; then
    echo "Usage: $0 -i <output directory from part 3> -k <list of taxonomic units to keep> [-u -c -d -m]"
    exit 1
fi

# Check if the output directory contains the directories "otus" and/or "asvs"
if [[ -d "${output_dir}/otus" && -d "${output_dir}/asvs" ]]; then
    echo "Error: Both 'otus' and 'asvs' directories exist in ${output_dir}. Exiting."
    exit 1
elif [[ -d "${output_dir}/otus" ]]; then
    typ="otu"
elif [[ -d "${output_dir}/asvs" ]]; then
    typ="asv"
else
    echo "Error: Neither 'otus' nor 'asvs' directories exist in ${output_dir}. Exiting."
    exit 1
fi

# initiate log
timestamp="$(date +"%y%m%d_%H:%M")"
output_file="${output_dir}/part5_${timestamp}.log"
exec > "$output_file"
exec 2> >(tee -a "$output_file" >&2)

# log header
echo "Log file for Part 5 of the Microbiome Pipeline. Processed the following arguments:
Output directory to be processed: ${output_dir}
List of taxonomic units to keep: ${KEEP}" | tee /dev/tty

if [ "$UNI" = true ]; then
    echo "The tables will be processed as universal assays and separated into 3 domains." | tee /dev/tty
fi

if [ "$CLA" = true ]; then
    echo "The classifier-generated taxonomic assignments will be used as primary assignments" | tee /dev/tty
fi

if [ "$DET" = true ]; then
    echo "A Detailed Informational OTU/ASV Table will be generated." | tee /dev/tty
fi

if [ "$MIX" = true ]; then
    echo "The mixed family analysis will be done." | tee /dev/tty
fi

echo " - -- --- ---- ---- --- -- -"

# detect if STRATEGY 3 and STRATEGY 2 were made and process if they are present.
STR2=false
STR3=false

# Check for STRATEGY2
if [ -d "${output_dir}/${typ}s/STRATEGY2" ]; then
  STR2=true
fi

# Check for STRATEGY3
if [ -d "${output_dir}/${typ}s/STRATEGY3" ]; then
  STR3=true
fi

# detect if classifier assignments are present and use them in final summary file if generated
if [ -d "${output_dir}/${typ}s/classifier_output" ]; then
	CP=true 
else
	if [ "$CLA" = true ]; then
		echo "Classifier-generated taxonomic assignments were not generated in Part 4. You cannot use them a primary assignments."
 		exit 1
 	fi
fi

source qiime_shell_helper_functions.sh || { echo "Error: Unable to source Qiime shell helper functions"; exit 1; }

# remove taxonomic units as indicated
OTBL=${typ}_table_01
biom subset-table \
  -i "${output_dir}/${typ}s/${OTBL}.biom" \
  -a observation \
  -s ${KEEP} \
  -o "${output_dir}/${typ}s/${typ}_table_02_TUs_removed.biom"

#add taxa to ASV table
OTBL=${typ}_table_02_TUs_removed
# Use classifier assignments if indicated; otherwise use BLAST.
if [ "$CLA" = true ]; then
	awk 'NR==1 {gsub("Feature ID", "#OTU ID"); gsub("Taxon", "taxonomy"); print; next} {print}' "${output_dir}/${typ}s/classifier_output/taxonomy.tsv" > "${output_dir}/${typ}s/classifier_output/taxonomy.fixed_headers.tsv"
	biomAddObservations "${output_dir}/${typ}s/${OTBL}.biom" "${output_dir}/${typ}s/${typ}_table_03_add_taxa.biom" "${output_dir}/${typ}s/classifier_output/taxonomy.fixed_headers.tsv"
	tail -n +2 "${output_dir}/${typ}s/classifier_output/taxonomy.tsv" | cut -d$'\t' -f1,2 > temp.txt
else
	biomAddObservations "${output_dir}/${typ}s/${OTBL}.biom" "${output_dir}/${typ}s/${typ}_table_03_add_taxa.biom" "${output_dir}/${typ}s/blast/tax_assignments.txt"
	tail -n +2 "${output_dir}/${typ}s/blast/tax_assignments.txt" | cut -d$'\t' -f1,2 > temp.txt
fi

# create three additional taxonomic levels of ASV tables
OTBL="${typ}_table_03_add_taxa"

# first, convert tables to qza format
qiime tools import \
  --input-path "${output_dir}/${typ}s/${OTBL}.biom" \
  --type 'FeatureTable[Frequency]' \
  --output-path "${output_dir}/${typ}s/${OTBL}.qza"

# convert taxonomy assignment file to qza format
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path temp.txt \
  --output-path "${output_dir}/${typ}s/taxonomy.qza"
rm -rf temp.txt

# generate levels
to_process=(
    "${output_dir}/${typ}s/${typ}_table_03_add_taxa_L2"
    "${output_dir}/${typ}s/${typ}_table_03_add_taxa_L6"
    "${output_dir}/${typ}s/${typ}_table_03_add_taxa_L7"
    "${output_dir}/${typ}s/${typ}_table_03_add_taxa"
)

for F in ${to_process[@]}; do
    if [[ "$F" =~ _L([0-9]+) ]]; then
        num="${BASH_REMATCH[1]}"
        qiime taxa collapse \
          --i-table "${output_dir}/${typ}s/${typ}_table_03_add_taxa.qza" \
          --i-taxonomy "${output_dir}/${typ}s/taxonomy.qza" \
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
      --output-path "${output_dir}/${typ}s/tmp_dir"
    mv "${output_dir}/${typ}s/tmp_dir/feature-table.biom" "${F}.biom"
    # repeat with norm
    qiime tools export \
      --input-path "${F}.norm.qza" \
      --output-path "${output_dir}/${typ}s/tmp_dir"
    mv "${output_dir}/${typ}s/tmp_dir/feature-table.biom" "${F}.norm.biom"
    rm -rf "${output_dir}/${typ}s/tmp_dir"
done

if [ "$CLA" = true ]; then
	biomAddObservations "${output_dir}/${typ}s/${typ}_table_03_add_taxa.norm.biom" "temp.txt" "${output_dir}/${typ}s/classifier_output/taxonomy.fixed_headers.tsv"
else
	biomAddObservations "${output_dir}/${typ}s/${typ}_table_03_add_taxa.norm.biom" "temp.txt" "${output_dir}/${typ}s/blast/tax_assignments.txt"
fi
mv temp.txt "${output_dir}/${typ}s/${typ}_table_03_add_taxa.norm.biom"

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
              --output-path "${output_dir}/${typ}s/tmp_dir"
            mv "${output_dir}/${typ}s/tmp_dir/feature-table.biom" "${F}.${K}.norm.biom"
            rm -rf "${output_dir}/${typ}s/tmp_dir"
        done  
    done
fi

for F in "${output_dir}/${typ}s/"*.biom; do
    txt_file="${F%.biom}.txt"
    if [ ! -f "$txt_file" ]; then
        biom2txt "$F" "$txt_file"
    fi
done

# add seqs to L8 
otblfp="${output_dir}/${typ}s/${typ}_table_03_add_taxa.txt"
outfp="${output_dir}/${typ}s/${typ}_table_04_add_seqs.txt"

add_sequences_to_asv.py ${otblfp} ${output_dir}/${typ}s/${typ}s.fa ${outfp}

otblfp="${output_dir}/${typ}s/${typ}_table_03_add_taxa.norm.txt"
outfp="${output_dir}/${typ}s/${typ}_table_04_add_seqs.norm.txt"

add_sequences_to_asv.py ${otblfp} ${output_dir}/${typ}s/${typ}s.fa ${outfp}

to_process2=($(find "${output_dir}/${typ}s" -maxdepth 1 -type f -name "${typ}_table_03*taxa.k*txt"))

for F in "${to_process2[@]}"; do
    if [ ! -f "$F" ]; then
        echo "Error: File $F not found."
        #exit 1
    fi
    FNAME=$(basename "$F" | sed 's|^${typ}_table_03_add_taxa||')
    otblfp="${F}"
    outfp="${output_dir}/${typ}s/${typ}_table_04_add_seqs${FNAME}"
    echo $outfp
    add_sequences_to_asv.py ${otblfp} ${output_dir}/${typ}s/${typ}s.fa ${outfp}
done

if [ "$DET" = true ]; then
    # Set default values for optional arguments
	mixed_arg="NULL"
	cp_arg="NULL"

    source pymods.sh || { echo "Error: Unable to activate python-pip-modules environment"; exit 1; }

    rm -rf "${output_dir}/${typ}s/blast/mixed_family_output.txt"
    # Run the mixed_family_checker.py script with error checking

    if [ "$MIX" = true ]; then
    	mixed_file="${output_dir}/${typ}s/blast/mixed_family_output.txt"
    	mixed_family_checker.py "${output_dir}/${typ}s/blast/final.blastout" --email "${EMAIL}" --output "${mixed_file}" || { echo "Error: mixed_family_checker.py failed"; exit 1; }
    	mixed_arg="${mixed_file}"
    fi

    # Define classifier variable if available 
    if [ "$CP" = true ]; then
  		cp_arg="${output_dir}/${typ}s/classifier_output/taxonomy.tsv"
	fi 

    # Define file paths
    norm_file="${output_dir}/${typ}s/${typ}_table_04_add_seqs.norm.txt"
    raw_file="${output_dir}/${typ}s/${typ}_table_04_add_seqs.txt"
    tax_file="${output_dir}/${typ}s/blast/tax_assignments.txt"
    sum_output="${output_dir}/${typ}s/Detailed_Informational_${typ}_Table.tsv"
    rm -rf "${sum_output}"
	
	# Single call to the script
	summary_file_generator.py ${norm_file} ${raw_file} ${tax_file} ${cp_arg} ${mixed_arg} ${sum_output}

    if [ ! -f ${sum_output} ]; then
        echo "${output_dir}/${typ}s/Detailed_Informational_${typ}_Table.tsv was not successfully created."
        #exit 1
    fi
fi

echo
echo " - -- --- ---- ---- --- -- -
Final Recommendations
 - -- --- ---- ---- --- -- -
The basic pipeline is complete. You now have both normalized and raw count tables for your 
taxonomic units. 

Note: Sometimes taxonomic units can be primarily associated with your controls (likely a 
source of contamination, which may be from other samples or even from the individual building 
the library in the first place). Make sure to check your tables for these units - just look at 
your control columns and see if any of the units are relatively high in them and not present 
in the regular samples. These units should probably be removed before further analysis.
" | tee /dev/tty