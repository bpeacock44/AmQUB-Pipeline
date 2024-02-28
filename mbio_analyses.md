# Some Basic Post-Run Analyses

These functions are available for some rudimentary analyses on the generated ASV tables. 

To run them, you will need a merged mapping file containing all the samples in your ASV table. You can create this by concatenating sample rows from all the mapping files you want to process. Note that it isn't problematic if the mapping files have MORE samples than the ASV table represents. These will just be ignored. But if you are missing samples from your mapping file, then the analyses will not be able to use the data from those samples. 

All of them except the correlation analysis are from the Qiime1 program, as this is the version of Qiime1 included in our program.
Keep in mind that the ASV tables created here can be used with Qiime2 for further analyes as well!
First, specify your asv tables (both regular and normalized) and your merged mapping file.
```sh 
reg_asv_table=asv_table_02_add_taxa.biom
norm_asv_table=asv_table_02_add_taxa_norm.biom
mapping_file=merged_map.txt
```

## SETUP
Run this section to set up your environment:
```sh 
HDIR=/home/bpeacock_ucr_edu/qiime1_setup/paul_helper_scripts
source /sw/anaconda3/etc/profile.d/conda.sh
conda activate qiime1
source ${HDIR}/qiime_shell_helper_functions.sh
module load r
```
Run the following to remove everything from your regular and normalized ASV tables except for those marked "SAMPLE" in the mapping file.
```sh 
filter_samples_from_OTU_table.py -i ${reg_asv_table} -o nc.${reg_asv_table} -m ${mapping_file} --valid_states "SampleType:SAMPLE"
filter_samples_from_OTU_table.py -i ${norm_asv_table} -o nc.${norm_asv_table} -m ${mapping_file} --valid_states "SampleType:SAMPLE"
```


## DIVERSITY ANALYSES
### BETA DIVERSITY (use raw counts)
First, create a metric parameter file.
```sh 
# Run this to see available metrics.
beta_diversity.py -s 
# Add or removed desired metrics to the variable "MET" below. They must have commas between them. 
B_MET="thellinger,abund_jaccard"
```
Now run the analysis. You can ignore the warning about "rank" being depreciated.
```sh 
echo -e 'beta_diversity:metrics\'${B_MET} > bd_parameters.txt # creates parameter file
beta_diversity_through_plots.py -p bd_parameters.txt -i nc.${reg_asv_table} -m ${mapping_file} -o bdivs.nc.${reg_asv_table} # runs analysis
cd bdivs.nc.${reg_asv_table}
condense_multiple_PCOA_plots
cd ..
```

### ALPHA DIVERSITY (use raw counts)
First, create a metric variable.
```sh
alpha_diversity.py -s # run this to see available metrics
# Add or removed desired metrics to the variable "MET" below. They must have commas between them. 
A_MET="chao1,shannon"
```
Now run the analysis.
```sh
# Now run the rest of this code without making changes. 
alpha_diversity.py -m ${A_MET} -i nc.${reg_asv_table} -o adivs.nc.${reg_asv_table}.txt
```


## OTHER ANALYSES
These next three require you to indicate which column of metadata in your mapping file you want to use to analyze the data. They must appear in the command exactly as they do in the mapping file, which must be called "merged_map.txt." (Capitalization, white space, etc.)
### DIFFERENTIAL ABUNDANCE (use raw counts)
Run as: 
${HDIR}/asv_diff_abun_wrapper.sh [ASV table path] [mapping file] [column of factors] [factor1] [factor2]
- column of factors - column containing traits you want to compare
- factor1 and factor2 - the two traits you want to compare
```sh
${HDIR}/asv_diff_abun_wrapper.sh ${norm_asv_table} ${mapping_file} tissue stem_whole stem_scrapings
```

## CORRELATION (use normalized counts)
Run as: 
${HDIR}/asv_diff_abun_wrapper.sh [ASV table path] [mapping file] [column of measurements to correlate against]
```sh
${HDIR}/asv_pearson_corr_wrapper.sh ${norm_asv_table} ${mapping_file} rating
```

## TAXA PRISM TABLE SUMMARY MAKER: (use normalized counts)
This last script will generate a csv file that can be used as input to the program prism. 

Run as: 
${HDIR}/asv_prism_summary_wrapper.sh [ASV table path] [mapping file] [column name] [comma-delimited treatment list] [number of taxa desired in final output]
- column name - the name of the column the treatments you want to splits samples by is located
- comma-delimited treatment list - a list of the treatments you want included in the summary
- number of taxa desired - the number of most abundant taxa you'd like in the summary (i.e. how many different groups will appear in your stacked bar plot)
Note that any taxa not included in that number will be lumped into an "other" category.
```sh
${HDIR}/asv_prism_summary_wrapper.sh ${norm_asv_table} ${mapping_file} tissue stem_whole,leaf_scrapings,leaf_whole 20
```

