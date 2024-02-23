#This script loads an OTU table, sorts columns and rows, then saves the result
#Usage:
# Rscript sort_otu_table.R

source('/sw/paul_helper_scripts/pipeline_helper_functions.R')

library("optparse")

args <- commandArgs(trailingOnly=TRUE)

parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input_otbl"), type="character", default="otu_table_00.txt",
        help="Input OTU table filepath [required].")
parser <- add_option(parser, c("-o", "--output_otbl"), type="character", default="otu_table_01.txt",
        help="Output OTU table filepath [required].")
parser <- add_option(parser, c("-f", "--force_overwrite"), action="store_true", default=FALSE,
        help="Force overwrite of output file [optional].")

opts <- parse_args(parser, args=args)


#options sanity check
if(!file.exists(opts$input_otbl)){
  stop(paste0('Cannot find [',opts$input_otbl,']\n'))
}
#prevent accidental overwrite
if(file.exists(opts$output_otbl) & !opts$force_overwrite) {
    stop(paste0('** Warning! ** Output file [',opts$output_otbl,'] exists! (Use "-f" to overwrite)'))
}

#load otu table
cat(paste0('Loading [',opts$input_otbl,']\n'))
tbl <- loadQIIMEotutable(opts$input_otbl)

cat('Sorting...\n')

#sort columns alphabetically
tbl <- sortQIIMEotutable(tbl, sortby="col", normalize_sort=FALSE)

#sort rows descending by rowSums
tbl <- sortQIIMEotutable(tbl, sortby="row", normalize_sort=FALSE)

#save
cat('Saving [',opts$output_otbl,']\n')
write.table(tbl, file=opts$output_otbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)