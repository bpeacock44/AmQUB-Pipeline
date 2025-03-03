#!/bin/bash
# See mbio_tutorial.md for further guidance!

### USAGE ###
# This script accepts two types of input:
# 1. Direct command-line arguments:
#    AmQUB_part5.sh -i output_dir -k to_keep.txt
#    AmQUB_part5.sh -i output_dir -k to_keep.txt -2 to_keep_S2.txt -3 to_keep_S3.txt -u -c -d -m  

## Required Flags
# -i: The output directory generated in part 3 that you also ran part 4 on.
# -k: A file of taxonomic units you wish to keep in the final table - one per line.

## Optional Flags:
# -2: A file of taxonomic units you wish to keep in the final table for STRATEGY 2 - one per line.
# -3: A file of taxonomic units you wish to keep in the final table for STRATEGY 3 - one per line.
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
#        Taxonomic Units to Keep File STRATEGY2,to_keep_S2.txt
#        Taxonomic Units to Keep File STRATEGY3,to_keep_S3.txt
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
            "Taxonomic Units to Keep File STRATEGY2,"*) STR2="${line#*,}" ;;
            "Taxonomic Units to Keep File STRATEGY3,"*) STR3="${line#*,}" ;;
            "Universal Assay,"*) UNI="${line#*,}" ;;
            "Classifier Assignments Primary,"*) CLA="${line#*,}" ;;
            "Detailed Informational OTU/ASV Table,"*) DET="${line#*,}" ;;
            "Mixed Top 10 Analysis,"*) MIX="${line#*,}" ;;
        esac
    done < "$param_file"
}

# Initialize default values
STR2=false
STR3=false
UNI=false
CLA=false
DET=false
MIX=false

echo "
        ┌── ===
 ┌──────┤
 │      └── ooo
─┤
 │ ┌── [A m Q U B]  
 └─┤
   └──── ~<>~   

Part 5: Post-Processing Tables

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
            -2|--strategy2) STR2=$2; shift 2 ;;
            -3|--strategy3) STR3=$2; shift 2 ;;
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
    echo "Usage: $0 -i <output directory from part 3> -k <list of taxonomic units to keep> [-2 -3 -u -c -d -m]"
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

# detect if STRATEGY 3 and STRATEGY 2 were made and process if they are present.
DIRS=("${output_dir}/${typ}s")

# check if input files exist

if [ "$STR2" != false ]; then 
    if [ -d "${output_dir}/${typ}s/STRATEGY2" ]; then
        DIRS+=("${output_dir}/${typ}s/STRATEGY2/${typ}s")
    else 
        echo "You have provided a taxonomic units to keep file for Strategy 2, but no Strategy 2 output directory (${output_dir}/${typ}s/STRATEGY2) was detected."
    fi
fi

if [ "$STR3" != false ]; then 
    if [ -d "${output_dir}/${typ}s/STRATEGY3" ]; then
        DIRS+=("${output_dir}/${typ}s/STRATEGY3/${typ}s")
    else 
        echo "You have provided a taxonomic units to keep file for Strategy 3, but no Strategy 3 output directory (${output_dir}/${typ}s/STRATEGY3) was detected."
    fi
fi

if [ ! -f "$KEEP" ]; then
    echo "The taxonomic units file $KEEP does not exist! Exiting." >&2
    exit 1
fi
if [ "$STR2" != false ] && [ ! -f "$STR2" ]; then
    echo "The taxonomic units file for Strategy 2 ($STR2) does not exist! Exiting." >&2
    exit 1
fi
if [ "$STR3" != false ] && [ ! -f "$STR3" ]; then
    echo "The taxonomic units file for Strategy 3 ($STR3) does not exist! Exiting." >&2
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
if [ "$STR2" != false ]; then 
    echo "Strategy 2 results will be processed."
fi

if [ "$STR3" != false ]; then 
    echo "Strategy 3 results will be processed."
fi

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

# Print the array to check the result
echo "Array contents: ${array[@]}"

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

for DIR in ${DIRS[@]}; do
    # remove taxonomic units as indicated
    OTBL=${typ}_table_01
    biom subset-table \
      -i "${DIR}/${OTBL}.biom" \
      -a observation \
      -s ${KEEP} \
      -o "${DIR}/${typ}_table_02_TUs_removed.biom"
    
    #add taxa to ASV table
    OTBL=${typ}_table_02_TUs_removed
    # Use classifier assignments if indicated; otherwise use BLAST.
    if [ "$CLA" = true ]; then
    	awk 'NR==1 {gsub("Feature ID", "#OTU ID"); gsub("Taxon", "taxonomy"); print; next} {print}' "${DIR}/classifier_output/taxonomy.tsv" > "${DIR}/classifier_output/taxonomy.fixed_headers.tsv"
    	biomAddObservations "${DIR}/${OTBL}.biom" "${DIR}/${typ}_table_03_add_taxa.biom" "${DIR}/classifier_output/taxonomy.fixed_headers.tsv"
    	tail -n +2 "${DIR}/${typ}s/classifier_output/taxonomy.tsv" | cut -d$'\t' -f1,2 > temp.txt
    else
    	biomAddObservations "${DIR}/${OTBL}.biom" "${DIR}/${typ}_table_03_add_taxa.biom" "${DIR}/blast/tax_assignments.txt"
    	tail -n +2 "${DIR}/blast/tax_assignments.txt" | cut -d$'\t' -f1,2 > temp.txt
    fi
    
    # create three additional taxonomic levels of ASV tables
    OTBL="${typ}_table_03_add_taxa"
    
    # first, convert tables to qza format
    qiime tools import \
      --input-path "${DIR}/${OTBL}.biom" \
      --type 'FeatureTable[Frequency]' \
      --output-path "${DIR}/${OTBL}.qza"
    
    # convert taxonomy assignment file to qza format
    qiime tools import \
      --type 'FeatureData[Taxonomy]' \
      --input-format HeaderlessTSVTaxonomyFormat \
      --input-path temp.txt \
      --output-path "${DIR}/taxonomy.qza"
    rm -rf temp.txt
    
    # generate levels
    to_process=(
        "${DIR}/${typ}_table_03_add_taxa_L2"
        "${DIR}/${typ}_table_03_add_taxa_L6"
        "${DIR}/${typ}_table_03_add_taxa_L7"
        "${DIR}/${typ}_table_03_add_taxa"
    )
    
    for F in ${to_process[@]}; do
        if [[ "$F" =~ _L([0-9]+) ]]; then
            num="${BASH_REMATCH[1]}"
            qiime taxa collapse \
              --i-table "${DIR}/${typ}_table_03_add_taxa.qza" \
              --i-taxonomy "${DIR}/taxonomy.qza" \
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
          --output-path "${DIR}/tmp_dir"
        mv "${DIR}/tmp_dir/feature-table.biom" "${F}.biom"
        # repeat with norm
        qiime tools export \
          --input-path "${F}.norm.qza" \
          --output-path "${DIR}/tmp_dir"
        mv "${DIR}/tmp_dir/feature-table.biom" "${F}.norm.biom"
        rm -rf "${DIR}/tmp_dir"
    done
    
    if [ "$CLA" = true ]; then
    	biomAddObservations "${DIR}/${typ}_table_03_add_taxa.norm.biom" "temp.txt" "${DIR}/classifier_output/taxonomy.fixed_headers.tsv"
    else
    	biomAddObservations "${DIR}/${typ}_table_03_add_taxa.norm.biom" "temp.txt" "${DIR}/blast/tax_assignments.txt"
    fi
    mv temp.txt "${DIR}/${typ}_table_03_add_taxa.norm.biom"
    
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
                  --output-path "${DIR}/tmp_dir"
                mv "${DIR}/tmp_dir/feature-table.biom" "${F}.${K}.norm.biom"
                rm -rf "${DIR}/tmp_dir"
            done  
        done
    fi
    
    for F in "${DIR}/"*.biom; do
        txt_file="${F%.biom}.txt"
        if [ ! -f "$txt_file" ]; then
            biom2txt "$F" "$txt_file"
        fi
    done
    
    # add seqs to L8 
    otblfp="${DIR}/${typ}_table_03_add_taxa.txt"
    outfp="${DIR}/${typ}_table_04_add_seqs.txt"
    
    add_sequences_to_asv.py ${otblfp} ${DIR}/${typ}s.fa ${outfp}
    
    otblfp="${DIR}/${typ}_table_03_add_taxa.norm.txt"
    outfp="${DIR}/${typ}_table_04_add_seqs.norm.txt"
    
    add_sequences_to_asv.py ${otblfp} ${DIR}/${typ}s.fa ${outfp}
    
    to_process2=($(find "${DIR}" -maxdepth 1 -type f -name "${typ}_table_03*taxa.k*txt"))
    
    for F in "${to_process2[@]}"; do
        if [ ! -f "$F" ]; then
            echo "Error: File $F not found."
            #exit 1
        fi
        FNAME=$(basename "$F" | sed 's|^${typ}_table_03_add_taxa||')
        otblfp="${F}"
        outfp="${DIR}/${typ}_table_04_add_seqs${FNAME}"
        echo $outfp
        add_sequences_to_asv.py ${otblfp} ${DIR}/${typ}s.fa ${outfp}
    done
    
    if [ "$DET" = true ]; then
        # Set default values for optional arguments
    	mixed_arg="NULL"
    	cp_arg="NULL"
    
        source pymods.sh || { echo "Error: Unable to activate python-pip-modules environment"; exit 1; }
    
        rm -rf "${DIR}/blast/mixed_family_output.txt"
        # Run the mixed_family_checker.py script with error checking
    
        if [ "$MIX" = true ]; then
        	mixed_file="${DIR}/blast/mixed_family_output.txt"
        	mixed_family_checker.py "${DIR}/blast/final.blastout" --email "${EMAIL}" --output "${mixed_file}" || { echo "Error: mixed_family_checker.py failed"; exit 1; }
        	mixed_arg="${mixed_file}"
        fi
    
        # Define classifier variable if available 
        if [ "$CP" = true ]; then
      		cp_arg="${DIR}/classifier_output/taxonomy.tsv"
    	fi 
    
        # Define file paths
        norm_file="${DIR}/${typ}_table_04_add_seqs.norm.txt"
        raw_file="${DIR}/${typ}_table_04_add_seqs.txt"
        tax_file="${DIR}/blast/tax_assignments.txt"
        sum_output="${DIR}/Detailed_Informational_${typ}_Table.tsv"
        rm -rf "${sum_output}"
	   
    	# Single call to the script
    	summary_file_generator.py ${norm_file} ${raw_file} ${tax_file} ${cp_arg} ${mixed_arg} ${sum_output}
    
        if [ ! -f ${sum_output} ]; then
            echo "${DIR}/Detailed_Informational_${typ}_Table.tsv was not successfully created."
            #exit 1
        fi
    fi
done

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