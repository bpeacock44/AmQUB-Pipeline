#This script loads an ASV table, sorts columns and rows, then saves the result
#Usage:
# Rscript sort_asv_table.R

source('/home/bpeacock_ucr_edu/real_projects/PN94_singularity_of_microbiome_pipeline/targeted_microbiome_via_blast/helper_functions/pipeline_helper_functions.R')

library("optparse")

args <- commandArgs(trailingOnly=TRUE)

parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input_atbl"), type="character", default="asv_table_00.txt",
        help="Input ASV table filepath [required].")
parser <- add_option(parser, c("-o", "--output_atbl"), type="character", default="asv_table_01.txt",
        help="Output ASV table filepath [required].")
parser <- add_option(parser, c("-f", "--force_overwrite"), action="store_true", default=FALSE,
        help="Force overwrite of output file [optional].")

opts <- parse_args(parser, args=args)


#options sanity check
if(!file.exists(opts$input_atbl)){
  stop(paste0('Cannot find [',opts$input_atbl,']\n'))
}
#prevent accidental overwrite
if(file.exists(opts$output_atbl) & !opts$force_overwrite) {
    stop(paste0('** Warning! ** Output file [',opts$output_atbl,'] exists! (Use "-f" to overwrite)'))
}

#load asv table
cat(paste0('Loading [',opts$input_atbl,']\n'))
tbl <- loadQIIMEasvtable(opts$input_atbl)

cat('Sorting...\n')

#sort columns alphabetically
tbl <- sortQIIMEasvtable(tbl, sortby="col", normalize_sort=FALSE)

#sort rows descending by rowSums
tbl <- sortQIIMEasvtable(tbl, sortby="row", normalize_sort=FALSE)

#save
cat('Saving [',opts$output_atbl,']\n')
write.table(tbl, file=opts$output_atbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)