# The original QIIME1 adonis.r code has been modified here
# See '## added code ##' sections to view code changes

# Runs vegan function adonis on QIIME distance matrix
# usage:
# R --slave --args --source_dir $QIIME_HOME/qiime/support_files/R/ -d unifrac.txt -m Fasting_Map.txt -c Treatment -o adonis < adonis.r
#
# print help string:
# R --slave --args -h --source_dir $QIIME_HOME/qiime/support_files/R/ < adonis.r
#
# Requires command-line param --source_dir pointing to QIIME R source dir

# load libraries and source files
args <- commandArgs(trailingOnly=TRUE)

if(!is.element('--source_dir', args)){
    stop("\n\nPlease use '--source_dir' to specify the R source code directory.\n\n")
}
sourcedir <- args[which(args == '--source_dir') + 1]
source(sprintf('%s/loaddata.r',sourcedir))
source(sprintf('%s/util.r',sourcedir))
load.library('optparse')
load.library('vegan')

# make option list and parse command line
option_list <- list(
    make_option(c("--source_dir"), type="character",
        help="Path to R source directory [required]."),
    make_option(c("-d", "--distmat"), type="character",
        help="Input distance matrix [required]."),
    make_option(c("-m", "--mapfile"), type="character",
        help="Input metadata mapping file [required]."),
    make_option(c("-c", "--category"), type="character",
        help="Metadata column header giving cluster IDs [required]."),
    make_option(c("-n", "--num_permutations"), type="integer", default=999,
        help="Number of permutations [default %default]."),
    # make_option(c("-o", "--outdir"), type="character", default='.',
        # help="Output directory [default %default]."),

    make_option(c("-o", "--outdir"), type="character",
        help="Output directory."),
    make_option(c("--output_fp"), type="character",
        help="Output filepath [default %default]."),

    make_option(c("-p", "--do_pairwise"), action="store_true", default=FALSE,
        help="Perform pairwise analyses between categories (instead of all-against-all)")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# error checking
if(is.null(opts$mapfile)) stop('Please supply a mapping file.')
if(is.null(opts$category)) stop('Please supply a mapping file header.')
if(is.null(opts$distmat)) stop('Please supply a distance matrix.')
if(is.null(opts$outdir) & is.null(opts$output_fp)) stop("Please use either '--outdir' or '--output_fp'")
if(!is.null(opts$outdir) & !is.null(opts$output_fp)) stop("Please use either '--outdir' or '--output_fp', but not both.")
#get output directory from output_fp. dirname returns "." if it has no dir
if(is.null(opts$outdir)) opts$outdir <- dirname(opts$output_fp)

# create output directory if needed
if(opts$outdir != ".") dir.create(opts$outdir,showWarnings=FALSE, recursive=TRUE)

# load qiime data
map <- load.qiime.mapping.file(opts$mapfile)

################
## added code ##
#load.qiime.mapping.file() may not work when the first line begins with a '#' (which they normally do),
#so we verify the map has loaded properly or try again
if(!"BarcodeSequence" %in% colnames(map)){
    map2 <- read.delim(opts$mapfile, header=TRUE, as.is=TRUE)
    if("BarcodeSequence" %in% colnames(map2)[2]){
        rownames(map2)=map2[,1]
        map <- map2[-1,]
    } else {
        cat("*** CAUTION *** (from adonis.r) The mapping file may not have loaded properly?\n")
    }
    cat("Mapping file reloaded okay IF this portion looks right:\n")
    cat("> map[1:3,1:2]\n")
    print(map[1:2,1:2])
}

#since we're using a map containing all samples we need to eliminate all rows not part of current category
smpls <- which(map[[opts$category]]!="")
map <- map[smpls,]
## added code ##
################

distmat <- load.qiime.distance.matrix(opts$distmat)

if(! opts$do_pairwise){
    # do qiime's default analysis for the current category #
    
    qiime.data <- remove.nonoverlapping.samples(map=map, distmat=distmat)

    # run adonis
    results <- adonis(
        formula = as.dist(qiime.data$distmat) ~ qiime.data$map[[opts$category]],
        permutations = opts$num_permutations)

    if(is.null(opts$output_fp)) {
        #qiime's default output filename
        filepath <- sprintf('%s/adonis_results.txt',opts$outdir)
    } else {
        #paul's optional output filepath
        filepath <- opts$output_fp
    }
    # write output file
    sink(filepath)
    print(results)
    sink(NULL)
} else {
    #################
    ## added code ##
    # do pairwise analyses for the current category #
    #save an original copy of the current map before pairwise subsetting begins
    map.orig <- map
    #get groups
    groups <- levels(factor(map.orig[[opts$category]]))
    #create pair combinations
    cmb <- c()
    if(length(groups) > 1) {
        cmb <- combn(groups,2)
    } else if(length(groups)==1) {
        stop(paste('\r*** Only 1 category-type [',groups,'] found in map column [',opts$category,'] ***', sep=""))
    }
    #get distance metric from dm filepath
    metric <- gsub('.*\\/(.*)_dm\\.txt','\\1', opts$distmat)
    
    
    for(i in 1:length(cmb[1,])){
        samples <- rownames(map.orig[which(map.orig[[opts$category]] %in% cmb[,i]),])
        pair <- sprintf('%s_v_%s',cmb[1,i],cmb[2,i])
        
        #prevent '__'s (from group names) getting into output filename (as [compare_categories2matrix.py] expects only 1 to parse filename correctly. See "need: <metric>..." below)
        pair <- gsub('__','_', pair)
        
        # define an output filepath
        #need: <metric>_<categ>__<pair1>_v_<pair2>_results.txt
        filepath <- sprintf('%s/%s_%s__%s_results.txt', opts$outdir, metric, opts$category, pair)
        # make output go there
        sink(filepath)
        
        #subset map
        map <- map.orig[samples,]
        #subset distmat using map
        qiime.data <- remove.nonoverlapping.samples(map=map, distmat=distmat)
        
        results <- adonis(
            formula = as.dist(qiime.data$distmat) ~ qiime.data$map[[opts$category]],
            permutations = opts$num_permutations)
        
        print(results)
        sink(NULL)
    }
    ## added code ##
    #################
}
