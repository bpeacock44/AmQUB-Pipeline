# Some Basic Post-Run Analyses

These functions are available for some rudimentary analyses on the generated ASV tables. 

To run them, you will need a merged mapping file containing all the samples in your ASV table. You can create this by concatenating sample rows from all the mapping files you want to process. Note that it isn't problematic if the mapping files have MORE samples than the ASV table represents. These will just be ignored. But if you are missing samples from your mapping file, then the analyses will not be able to use the data from those samples. 

All of them except the correlation analysis are from the Qiime1 program, as this is the version of Qiime1 included in our program.
Keep in mind that the ASV tables created here can be used with Qiime2 for further analyes as well!

*Make sure you are within the container before you start.*

Specify your asv tables (both regular and normalized) and your merged mapping file.

```sh 
reg_asv_table=asv_table_02_add_taxa.biom
norm_asv_table=asv_table_02_add_taxa_norm.biom
mapping_file=merged_map.txt
```

Note that when using the mapping file with Qiime1, it has to be formatted correctly. The following command can check for issues/errors with your mapping file before beginning:

```sh
validate_mapping_file.py -m ${mapping_file} -o validation_output/
```

## SETUP
Run the following to remove everything from your regular and normalized ASV tables except for those marked "SAMPLE" in the mapping file.
```sh 
filter_samples_from_otu_table.py -i ${reg_asv_table} -o nc.${reg_asv_table} -m ${mapping_file} --valid_states "SampleType:SAMPLE"
filter_samples_from_otu_table.py -i ${norm_asv_table} -o nc.${norm_asv_table} -m ${mapping_file} --valid_states "SampleType:SAMPLE"
```

## DIVERSITY ANALYSES
### BETA DIVERSITY (use raw counts)
First, create a metric parameter file.
```sh 
# Run this to see available metrics.
beta_diversity.py -s 
# Add or removed desired metrics to the variable "MET" below. They must have commas between them. 
B_MET="thellinger,abund_jaccard"
# Create the file.
echo -e 'beta_diversity:metrics\'${B_MET} > bd_parameters.txt # creates parameter file
```
Now run the analysis. You can ignore the warning about "rank" being depreciated.
```sh 
beta_diversity_through_plots.py -p bd_parameters.txt -i nc.${reg_asv_table} -m ${mapping_file} -o bdivs.nc.${reg_asv_table} # runs analysis
source qiime_shell_helper_functions.sh || { echo "Error: Unable to source Qiime shell helper functions"; exit 1; } # loads some helper functions
cd bdivs.nc.${reg_asv_table}; condense_multiple_PCOA_plots; cd -
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
asv_diff_abun_wrapper.sh -t [ASV table path in biom format] -m [mapping file] -c [column of factors] -v1 [factor1] -v2 [factor2] 
- column of factors - column containing traits you want to compare
- factor1 and factor2 - the two traits you want to compare
```sh
# create a text version of your raw count table
source qiime_shell_helper_functions.sh || { echo "Error: Unable to source Qiime shell helper functions"; exit 1; } # loads some helper functions
biom2txt "nc.${reg_asv_table}" "nc.${reg_asv_table%.biom}.txt"
# run the analysis (edgeR)
asv_diff_abun_wrapper.sh -t "nc.${reg_asv_table}" -m ${mapping_file} -c tissue -v1 stem_whole -v2 stem_scrapings 
```

## CORRELATION (use normalized counts)
Run as: 
asv_diff_abun_wrapper.sh [ASV table path] [mapping file] [column of measurements to correlate against]
```sh
asv_pearson_corr_wrapper.sh "nc.${norm_asv_table}" ${mapping_file} rating
```

## TAXA PRISM TABLE SUMMARY MAKER: (use normalized counts)
This last script will generate a csv file that can be used as input to the program prism. 

Run as: 
asv_prism_summary_wrapper.sh [ASV table path] [mapping file] [column name] [comma-delimited treatment list] [number of taxa desired in final output]
- column name - the name of the column the treatments you want to splits samples by is located
- comma-delimited treatment list - a list of the treatments you want included in the summary
- number of taxa desired - the number of most abundant taxa you'd like in the summary (i.e. how many different groups will appear in your stacked bar plot)
Note that any taxa not included in that number will be lumped into an "other" category.
```sh
asv_prism_summary_wrapper.sh ${norm_asv_table} ${mapping_file} tissue stem_whole,leaf_scrapings,leaf_whole 20
```

