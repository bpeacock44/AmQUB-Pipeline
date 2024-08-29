
find_last_header_line <- function(fp) {
#finds the first line (N) starting WITHOUT a '#' char and RETURNS N-1
    con <- file(fp, "r")
    N = 0
    while(TRUE){
        line <- readLines(con, n=1)
        X <- grep("^#", line, perl=TRUE, value=FALSE)
        if(length(X) == 0){
            break
        } else {
            N = N + 1
        }
    }
    close(con)
    return(N)
}

loadQIIMEmap <- function(fp) {
#loads a QIIME-happy mapping file
    skip <- find_last_header_line(fp) - 1
    map <- read.table(file=fp, skip=skip, sep="\t", header=TRUE, as.is=TRUE, comment.char="")
    colnames(map)[1] <- "#SampleID"
    rownames(map) <- map[[1]]
    return(map)
}

loadQIIMEasvtable <- function(fp) {
#loads a QIIME-happy text asv table file
    skip <- find_last_header_line(fp) - 1
    tbl <- read.table(file=fp, skip=skip, sep="\t", header=TRUE, fill=TRUE, comment.char="", stringsAsFactors=FALSE)
    colnames(tbl)[1] <- "#ASV ID"
    rownames(tbl) <- tbl[[1]]
    return(tbl)
}

seq_imp_fct <- function(fileloc) {
#Loads a fasta file and formats specifically for including into edgeR results. Eg: >denovo10_3209#ACTGG...
#Modified from: http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/sequenceAnalysis.txt
    my_fasta <- readLines(fileloc) # reads file line-wise into vector
    y <- regexpr("^[^>]", my_fasta, perl=T) # identifies all fields that DO NOT start with a '>' sign
    y <- as.vector(y); y[y==-1] <- 0; #change all -1's to 0's
    # The 'rep' function is used here to generate a grouping vector/column for writing each sequence into a single field (record?)
    index <- which(y==0); #get indexes of fasta name lines
    distance <- data.frame(start=index[1:(length(index)-1)], end=index[2:length(index)]); #create sequence start-end table (seqs can be on multiple lines)
    distance <- rbind(distance, c(distance[length(distance[,1]),2], length(y)+1)) # gets data for last entry
    distance <- data.frame(distance, dist=distance[,2]-distance[,1]); #add calculated "distance" in a 3rd column
    seq_no <- 1:length(y[y==0]); #[1] 1 2 3 4 5 ...
    index <- rep(seq_no, as.vector(distance[,3])); #[1] 1 1 2 2 3 3 4 4 ... (ie., 3 3 says the 3rd seq can be found at my_fasta[5:6])
    my_fasta <- data.frame(index, y, my_fasta)
    my_fasta[my_fasta[,2]==0,1] <- 0; #change vals in col_1 to 0 IF vals in col_2 are 0
    # Use the 'paste' function in a 'tapply' loop to combine each sequence (that may be on multiple lines) into a single field/line
    seq <- tapply(as.vector(my_fasta[,3]), factor(my_fasta[,1]), paste, collapse="", simplify=F)
    seq <- as.vector(seq[2:length(seq)])
    Desc <- as.vector(my_fasta[c(grep("^>", as.character(my_fasta[,3]), perl = TRUE)),3])

    # modified code starts here #
    Desc2 = gsub("\\s","_",Desc)
    fa2 <- data.frame(Desc=Desc2, Seq=seq)
    fa_names = sub(">", "", Desc2, perl=TRUE)
    fa_names2 = sub("_\\d+", "", fa_names, perl=TRUE)
    rownames(fa2) <- fa_names2
    return(fa2)
}

sortQIIMEasvtable <- function(tbl, sortby="row", normalize_sort=FALSE) {
#returns tbl sorted descending by rowSums or ascending by colSums
    #determine which columns have only numbers
    if(typeof(sortby) != "character") {
        stop("Values for 'sortby' must be either 'row' or 'col' (and be in quotes)")
    }
    column_classes <- unlist(lapply(tbl,class))
    numcols <- which(column_classes == "integer" | column_classes == "numeric")
    if (length(numcols) == 0) {
        stop("No 'integer' or 'numeric' columns found in ASV table! Cannot sort!")
    }
    if(sortby=="row") {
        if(normalize_sort == TRUE) {
            nc <- tbl[,numcols]
            nc <- t(t(nc)/colSums(nc))
            x <- rowSums(nc)
        } else {
            x <- rowSums(tbl[,numcols])
        }
        return(tbl[order(-x),])
    }
    else if(sortby=="col") {
        nc <- tbl[,numcols]
        nd <- nc[,order(colnames(nc))]
        tbl[,numcols] <- nd
        names(tbl)[numcols] <- colnames(nd)
        return(tbl)
    } else {
        stop("Values for 'sortby' must be either 'row' or 'col'")
    }
}

normalize_asvtable <- function(tbl) {
    original_class <- class(tbl)
    #ensure we have a dataframe so sapply works properly
    tbl <- as.data.frame(tbl)
    #determine which columns have only numbers
    column_classes <- sapply(tbl,class)
    numcols <- which(column_classes == "integer" | column_classes == "numeric")
    if (length(numcols) == 0) {
        stop("No 'integer' or 'numeric' columns found in ASV table! Cannot normalize!")
    }
    nc <- tbl[,numcols]
    nc <- t(t(nc)/colSums(nc))
    tbl[,numcols] <- nc
    #convert back to matrix?
    if(original_class == "matrix") {
        tbl <- as.matrix(tbl)
    }
    return(tbl)
}

add_counts_to_fasta_sequences <- function(otbl_fp=NULL, fasta_fp=NULL, out_fp=NULL) {
    if(is.null(otbl_fp) | is.null(fasta_fp) | is.null(out_fp)) {
        stop("Usage is:\nadd_counts_to_fasta_sequences(otbl_fp, fasta_fp, out_fp)")
    }
    #load asv table
    tbl<-loadQIIMEasvtable(otbl_fp)
    #load fasta sequences
    fastaSeqs<-seq_imp_fct(fasta_fp)
    dim(fastaSeqs)
    dim(tbl)
    tbl[1:3,1:5]
    #order table descending by rowSums
    x <- c(which(sapply(tbl,class) == "integer"))
    tbl2 <- tbl[order(-rowSums(tbl[,x])),x]
    #put sequences in the same order
    fS<-fastaSeqs[rownames(tbl2),]
    #get rowSums
    rs<-rowSums(tbl2)
    #make asv names with abundances
    rn<-paste0(">",names(rs)," ",rs)
    #check if rownames agree
    if(! all(names(rs)==rownames(fS)) ) {
        stop("all(names(tbl) != rownames(fasta))")
    }
    #replace sequence names
    fS[,1]<-rn
    #join/save
    ASV_seq <- unlist(c(rbind(fS[,1],fS[,2])))
    write.table(ASV_seq, file=out_fp, sep="", row.names=FALSE, col.names=FALSE, quote=FALSE)
}

add_sequences_to_asv_table <- function(tblfp=NULL, fastafp=NULL, outfp=NULL) {
    usage <- "Usage: add_sequences_to_asv_table(<tablefp>, <fastafp>, <outfp>)"
    if(is.null(tblfp) | is.null(fastafp) | is.null(outfp)) {
        stop(usage)
    }
    errmessg <- ''
    if(!file.exists(tblfp)) {
        errmessg <- paste0("File [", tblfp,"] does not exist!\n")
    }
    if(!file.exists(fastafp)) {
        errmessg <- paste0(errmessg, paste0("File [", fastafp,"] does not exist!\n"))
    }
    if(nchar(errmessg) > 0) {
        cat(usage,"\n")
        stop(errmessg)
    }
    #load asv table
    tbl <- loadQIIMEasvtable(tblfp)
    #load fasta sequences
    fastaSeqs <- seq_imp_fct(fastafp)
    #put sequences in the same order
    fS <- fastaSeqs[match(rownames(tbl),rownames(fastaSeqs)),]
    #check if rownames agree
    if(any(rownames(tbl) != rownames(fS))) {
        stop("all(rownames(tbl) != rownames(fS))")
    }
    #replace sequence names (that may contain abundances)
    fS[,1] <- rownames(fS)
    #join/save
    ASV_seq <- unlist(c(rbind(paste0(">",fS[,1],"#",fS[,2]))))
    tbl2 <- cbind(tbl, sequence=ASV_seq)
    write.table(tbl2, file=outfp, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}

# Added by Beth Peacock
process_data_and_write_summary <- function(norm, raw, tax, tax_c = NULL, mixed) {
    # Load required libraries
    library(dplyr)

    # Read the main data file
    df <- read.table(norm, sep = "\t", comment.char = "", header = TRUE)
    colnames(df)[1] <- "ASV_ID"
    df$taxonomy <- NULL
    
    # Read primary taxonomy file
    tax_df <- read.table(tax, sep = "\t", comment.char = "", header = TRUE)
    colnames(tax_df) <- c("ASV_ID", "taxonomy", "bitscore", "per_ID", "per_qcov")
    df <- merge(df, tax_df[c(1:2, 4:5)], by = "ASV_ID", all.x = TRUE)


    # Check if tax_c is provided and process it if available
    if (!is.null(tax_c)) {

        # Read classifier taxonomy file if given
        tax_df_C <- read.table(tax_c, sep = "\t", comment.char = "", header = TRUE)
        colnames(tax_df_C) <- c("ASV_ID", "c_taxonomy", "confidence")
        df <- merge(df, tax_df_C, by = "ASV_ID", all.x = TRUE)
        
        # Include c_taxonomy and confidence in excluded columns
        excluded_columns <- c("ASV_ID", "taxonomy", "per_ID", "per_qcov", "c_taxonomy", "confidence", "sequence")
    } else {

        # Exclude c_taxonomy and confidence if tax_c is not provided
        excluded_columns <- c("ASV_ID", "taxonomy", "per_ID", "per_qcov", "sequence")
    }
    
    # Calculate mean abundance
    df$avg_abun <- rowMeans(df[, !names(df) %in% excluded_columns], na.rm = TRUE)
    
    # Flag ASVs with mixed families in top 10
    top10_df <- read.table(mixed, sep = "\t", comment.char = "", header = TRUE)
    df$mixed <- "no"
    top10_identifiers <- top10_df[[1]]
    matches <- df$ASV_ID %in% top10_identifiers
    df$mixed[matches] <- "yes"
    
    # Reorder columns
    if (!is.null(tax_c)) {
        df <- df[, c(1, (length(df) - 7):length(df), 2:(length(df) - 8))]
    } else {
        df <- df[, c(1, (length(df) - 5):length(df), 2:(length(df) - 6))]
    }
    # Add row with total raw counts
    raw <- read.table(raw, sep = "\t", comment.char = "", header = TRUE)
    raw[,1] <- NULL
    raw$taxonomy <- NULL
    raw$sequence <- NULL
    sums <- colSums(raw)
    if (!is.null(tax_c)) {
        nd_columns <- rep("nd", 9)
    } else {
        nd_columns <- rep("nd", 7)
    }
    new_vector <- c(nd_columns, sums)
    df <- rbind(new_vector, df)
    
    # Write dataframe to an Excel file
    write.table(df, "ASV_summary_table.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}

