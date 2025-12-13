```sh
######## ######## ######## ######## ########
### ANALYSES
# make maximum parsimony trees
WDIR=/bigdata/bornemanlab/bpeacock/PN110_param_mbio_pipeline/250730_JB236JB239_test_with_James
ALLDIRS=("part3_table1_2_out/otus"
"part3_table1_2_out/otus/STRATEGY2/otus")
ALLDIRS=("part3_table11_out/otus/STRATEGY3/otus")
NUMTHREADS=60
#module load mafft

for DIR in "${ALLDIRS[@]}"; do
    rm -rf "${WDIR}/${DIR}/finding_more_hovis/final_putative_hovis.aln.fa"
    mafft --thread ${NUMTHREADS} --auto "${WDIR}/${DIR}/finding_more_hovis/final_putative_hovis.fa" > "${WDIR}/${DIR}/finding_more_hovis/final_putative_hovis.aln.fa"
done

module load emboss
for DIR in ${ALLDIRS[@]}; do
    # Extract descriptive components from directory path
    TABLE=$(echo "$DIR" | grep -o 'part3_table[0-9]\+')
    STRATEGY=$(echo "$DIR" | grep -o 'STRATEGY[0-9]\+')
    # If STRATEGY is empty, set it to "STRATEGY1"
    if [[ -z "$STRATEGY" ]]; then
        STRATEGY="STRATEGY1"
    fi
    mkdir -vp "${WDIR}/${DIR}/finding_more_hovis/analyses/tree"
    # Construct a more descriptive output file name
    OUTTREEFILE="${WDIR}/${DIR}/finding_more_hovis/analyses/tree/${TABLE}_${STRATEGY}_final_tree.fix.nwk"
    # Run fdnapars with the new output tree file name
    fdnapars -sequence "${WDIR}/${DIR}/finding_more_hovis/final_putative_hovis.aln.fa"  \
        -outfile "${WDIR}/${DIR}/finding_more_hovis/analyses/tree/phylip" \
        -outtreefile "$OUTTREEFILE" \
        -auto
done

# correlation
WDIR=/bigdata/bornemanlab/bpeacock/PN110_param_mbio_pipeline/250730_JB236JB239_test_with_James
map="${WDIR}/data/merged_map.txt"
run_analysis() {
    local -n DIRS_REF=$1  # Accepts an array of directories
    local -n COLUMNS_REF=$2  # Accepts an array of columns

    for DIR in "${DIRS_REF[@]}"; do
        table="${WDIR}/${DIR}/finding_more_hovis/new_output_files/otu_table_03_add_taxa.norm.txt"
        out="${WDIR}/${DIR}/finding_more_hovis/analyses/corr"
        mkdir -vp "${out}"
        for COL in "${COLUMNS_REF[@]}"; do
            corr_otu_v_col.py --asv_file "${table}" \
               --metadata_file "${map}" \
               --column "${COL}" \
               --output_dir "${out}"
        done
    done
}

# 1. Cyst Taxa vs Nematode Suppression (Corr1)
# 2. Female Taxa vs Nematode Suppression (Corr2)
# 3. Cyst + Female Taxa vs Nematode Suppression (Corr3)
columns=(Corr1 Corr2 Corr3) 
DIRS=("part3_table1_2_out/otus"
"part3_table1_2_out/otus/STRATEGY2/otus"
"part3_table11_out/otus/STRATEGY3/otus")
run_analysis DIRS columns

# II. Correlations of Taxa vs Taxa # UNFINISHED - NEED TO DISCUSS STRATEGY FOR DECREASING SIZE OF OUTPUT
#WDIR=/bigdata/bornemanlab/bpeacock/PN110_param_mbio_pipeline/250730_JB236JB239_test_with_James
WDIR=/bigdata/bornemanlab/bpeacock/PN110_param_mbio_pipeline/250519_JB236JB239_big
map="${WDIR}/data/merged_map.txt"
run_analysis() {
    local -n DIRS_REF=$1  # Accepts an array of directories
    local -n columns_REF=$2  # Accepts an array of columns
    local setA="$3"
    local setB="$4"

    for DIR in "${DIRS_REF[@]}"; do
        table="${WDIR}/${DIR}/finding_more_hovis/new_output_files/otu_table_03_add_taxa.norm.txt"
        out="${WDIR}/${DIR}/finding_more_hovis/analyses/corr"
        mkdir -vp "${out}"
        for COL in "${columns_REF[@]}"; do
            ./otu.py --asv_file ${table} \
                --metadata_file ${map} --output_dir ${out} \
                --column ${COL} --corr_threshold_high "0.3" \
                --corr_threshold_low "-0.3" --pval_threshold "0.05" \
                --setA ${setA} --setB ${setB} \
                #--min_avg_abundance "0.0005"
        done
    done
}

# 8. Cyst + Female Taxa vs Cyst + Female Taxa (Corr8)
columns=(Corr8)  
setA="CFtaxa"
setB="CFtaxa"
DIRS=("part3_table1_2_out/otus")
#"part3_table1_2_out/otus/STRATEGY2/otus")
run_analysis DIRS columns "$setA" "$setB"
# started at 1:30pm on Dec 5 with no min_avg_abundance which is 16998 OTUs

# III. Correlations of Fold Increase of Soil vs Taxa
WDIR=/bigdata/bornemanlab/bpeacock/PN110_param_mbio_pipeline/250730_JB236JB239_test_with_James
map="${WDIR}/data/merged_map.txt"
run_analysis() {
    local -n DIRS_REF=$1  # Accepts an array of directories
    local -n columns_REF=$2  # Accepts an array of columns
    local setA="$3"
    local setB_before="$4"
    local setB_after="$5"

    for DIR in "${DIRS_REF[@]}"; do
        table="${WDIR}/${DIR}/finding_more_hovis/new_output_files/otu_table_03_add_taxa.norm.txt"
        out="${WDIR}/${DIR}/finding_more_hovis/analyses/corr"
        mkdir -vp "${out}"
        for COL in "${columns_REF[@]}"; do
            corr_otu_v_foldincr.py --asv_file ${table} \
                --metadata_file ${map} --output_dir ${out} \
                --column ${COL} --corr_threshold_high "0.3" \
                --corr_threshold_low "-0.3" --pval_threshold "0.05" \
                --treatment_col "SoilNumber" --setA ${setA} --setB_before ${setB_before} \
                --setB_after ${setB_after} 
        done
    done
}

# 18. Fold Increase of Soil vs. Cyst Taxa (Corr18)
columns=(Corr18) 
setA="Ctaxa" 
setB_before="BStaxa"
setB_after="AStaxa"
DIRS=("part3_table11_out/otus/STRATEGY3/otus")
run_analysis DIRS columns ${setA} ${setB_before} ${setB_after}

# IV. Prism Plots
WDIR=/bigdata/bornemanlab/bpeacock/PN110_param_mbio_pipeline/250730_JB236JB239_test_with_James
map="${WDIR}/data/merged_map.txt"
run_analysis() {
    local -n DIRS_REF=$1
    local -n columns_REF=$2
    local -n num_ref=$3

    # Generate PRISM-compatible file
    for DIR in "${DIRS_REF[@]}"; do
        ASV_NORM="${WDIR}/${DIR}/finding_more_hovis/new_output_files/otu_table_04_add_seqs.norm.txt"
        out="${WDIR}/${DIR}/finding_more_hovis/analyses/prism"
        mkdir -vp ${out}
        for col in "${columns_REF[@]}"; do
            OUTPUT_PRISM_FILE="${out}/prism.${col}.tsv"
            OUTPUT_AVG_ABUNDANCE_FILE="${out}/prism.${col}.avgabun"
            mkdir -vp $OUTPUT_AVG_ABUNDANCE_FILE
            ./prism_maker_sep.py \
                --file_path "$ASV_NORM" \
                --map_file "${map}" \
                --column "${col}" \
                --output_avg_abundance_directory "$OUTPUT_AVG_ABUNDANCE_FILE" \
                --output_prism_file "$OUTPUT_PRISM_FILE" \
                --num_taxa ${num_ref}
            rm -rf ${OUTPUT_AVG_ABUNDANCE_FILE}
        done
    done
}

# 1-5 here aren't necessary - they are generated automatically in 6-10 as the "all" row.
# 1. Cyst Taxa All Soils Plots (Taxa1)
num=3
columns=(Taxa1)
DIRS=("part3_table1_2_out/otus"
"part3_table1_2_out/otus/STRATEGY2/otus"
"part3_table11_out/otus/STRATEGY3/otus")
run_analysis DIRS columns num

# 3. Cyst + Female Taxa All Soils Plots (Taxa3)
num=3 # NOT WORKING?
columns=(Taxa3)
DIRS=("part3_table1_2_out/otus"
"part3_table1_2_out/otus/STRATEGY2/otus"
"part3_table11_out/otus/STRATEGY3/otus")
run_analysis DIRS columns num

# 8. Cyst + Female Taxa Separate Soils Plots (Taxa8)
num=15
columns=(Taxa8)
DIRS=("part3_table1_2_out/otus"
"part3_table1_2_out/otus/STRATEGY2/otus"
"part3_table11_out/otus/STRATEGY3/otus")
run_analysis DIRS columns num

# 14. Before Soil Separate Samples Plots (Taxa14)
num=15
columns=(Taxa14)
DIRS=("part3_table11_out/otus/STRATEGY3/otus")
run_analysis DIRS columns num





# collect results to send
ALLDIRS=("part3_table1_2_out/otus"
"part3_table1_2_out/otus/STRATEGY2/otus"
"part3_table11_out/otus/STRATEGY3/otus")
for F in ${ALLDIRS[@]}; do
    mkdir -vp check/${F}
    #cp -r ${F}/finding_more_hovis/phylip check/${F}/phylip
    cp -r ${F}/finding_more_hovis/new_output_files/analyses/prism check/${F}/prism
    cp -r ${F}/finding_more_hovis/new_output_files/Det* check/${F}
done
```