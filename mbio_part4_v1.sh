#mbio_part4_v1.sh

# These functions are available for some rudimentary analyses on the generated OTU tables.
# To run them, you will need a merged mapping file containing all the samples in your OTU table.
# You can create this by concatenating sample rows from all the mapping files you want to process. 
# Note that it isn't problematic if the mapping files have MORE samples than the OTU table represents. 
# These will just be ignored. But if you are missing samples from your mapping file, then the analyses 
# will not be able to use the data from those samples. 

# These analysis will use the normalized and non-normalized versions of otu_table_02_add_taxa, which is 
# the OTU table before more metadata is added besides taxonomy. All but correlation are from the Qiime1 program.
HDIR=paul_helper_scripts
source /sw/anaconda3/etc/profile.d/conda.sh
conda activate qiime1
source ${HDIR}/qiime_shell_helper_functions.sh
module load r

##### ##### ##### ##### ##### 
# First run the following to remove everything from your OTU tables except for those marked "SAMPLE" in 
# the mapping file.
filter_samples_from_otu_table.py -i otu_table_02_add_taxa.biom -o otu_table_02_add_taxa.nc.biom -m merged_map.txt --valid_states "SampleType:SAMPLE"
filter_samples_from_otu_table.py -i otu_table_02_add_taxa_norm.biom -o otu_table_02_add_taxa_norm.nc.biom -m merged_map.txt --valid_states "SampleType:SAMPLE"





##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
# These first two are diversity analyses, and I have left them outside of a typical script in case you want to
# play around with the parameters (the type of analysis, etc.)

##### ##### ##### ##### ##### 
# BETA DIVERSITY (normalization happens internally)
# Create parameter file.
echo -e 'beta_diversity:metrics\thellinger,abund_jaccard\n#,weighted_unifrac,unweighted_unifrac\n' > bd_parameters.txt
beta_diversity_through_plots.py -p bd_parameters.txt -i otu_table_02_add_taxa.nc.biom -m merged_map.txt -o bdivs.otu_table_02_add_taxa.nc.biom
cd bdivs.otu_table_02_add_taxa.nc.biom
condense_multiple_PCOA_plots
cd ..

##### ##### ##### ##### ##### 
# ALPHA DIVERSITY (normalization happens internally)
# TO BE WRITTEN





##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
# These next two require you to indicate which column of metadata in your mapping file you want to use to 
# analyze the data. They must appear in the command exactly as they do in the mapping file 
# (Capitalization, white space, etc.)

##### ##### ##### ##### ##### 
# DIFFERENTIAL ABUNDANCE (normalization happens internally)
# Run as: 
# ${HDIR}/otu_diff_abun_wrapper.sh <column of factors> <factor1> <factor2>

${HDIR}/otu_diff_abun_wrapper.sh Tissue stem_whole stem_scrapings

##### ##### ##### ##### ##### 
# CORRELATION (pre-normalized file required)
# Run as: 
# ${HDIR}/otu_diff_abun_wrapper.sh <column of measurments to correlate against>

${HDIR}/otu_pearson_corr_wrapper.sh Rating
