source("lib_rnaseq.R")
require(rtracklayer)

if (count_tool == "STAR"){
  ## Load STAR data
    message("loading STAR gene counts ...")
    dirlist <- list.dirs(path=input_path, recursive=FALSE, full.names=TRUE)
    ## remove 'report' folder
    dirlist <- dirlist[which(! basename(dirlist) %in% c("export", "logs", "report"))]
  
    counts.exprs <- lapply(dirlist, function(d){
        exprs.in <- list.files(path=d, pattern="_counts.csv$", full.names=TRUE, recursive=TRUE)
        if (length(exprs.in) == 1){
            counts.exprs <- read.csv(exprs.in, sep="\t", header=FALSE, row.names=1, check.names=FALSE)
            
            ## Get stranded info
            if (stranded == ''){
                strandfile <- list.files(path=d, pattern="strandness.info", full.names=TRUE, recursive=TRUE)
                stopifnot(length(strandfile)==1)
                stranded <- gsub("STRANDED=","",read.table(strandfile)[1,1])
                message("Sample=", basename(d), " - strandness=", stranded)
            }
            ## Get right column
            if (stranded == "reverse"){
                counts.exprs <- counts.exprs[,3, drop=FALSE]
            }else if (stranded == "yes"){
                counts.exprs <- counts.exprs[,2, drop=FALSE]
            }else if (stranded == "no"){
                counts.exprs <- counts.exprs[,1, drop=FALSE]
            }else{
                stop("STRANDED parameter is not defined")
            }           
            return(counts.exprs)
        }
    })

    counts.exprs <- data.frame(counts.exprs)
    colnames(counts.exprs) <- sapply(dirlist, basename)
    ## remove first 4 lines
    counts.exprs <- counts.exprs[5:nrow(counts.exprs), , drop=FALSE]
}else if (count_tool == "FEATURECOUNTS"){
    ## Load FeatureCounts data
    exprs.in <- list.files(path=input_path, pattern="counts.csv$", full.names=TRUE, recursive=TRUE)
    counts.exprs <- lapply(exprs.in, function(f){
        z <- read.csv(f, sep="\t", row.names=1, comment.char="#", check.names=FALSE)
        data.frame(z[,6], row.names=rownames(z))
    })
    counts.exprs <- data.frame(counts.exprs)
    colnames(counts.exprs)<- gsub("_norRNA", "", gsub("_counts.csv","",sapply(exprs.in, basename)))
    
}else if (count_tool == "HTSEQ"){
    ## Load HTSeq data
  exprs.in <- list.files(path=input_path, pattern="_counts.csv$", full.names=TRUE, recursive=TRUE)
  counts.exprs <- lapply(exprs.in, read.csv, sep="\t", header=FALSE, row.names=1, check.names=FALSE)
  counts.exprs <- as.data.frame(counts.exprs)
  colnames(counts.exprs) <- gsub("_norRNA", "", gsub("_counts.csv","",sapply(exprs.in, basename)))
  counts.exprs <- counts.exprs[1:(nrow(counts.exprs)-5), , drop=FALSE]  
}

## export count table(s)
write.csv(ensembl2symbol(counts.exprs, gtf), file=file.path(odir, "tablecounts_raw.csv"))

## TPM
exonic.gene.sizes <- getExonicGeneSize(gtf)
dtpm <- getTPM(counts.exprs, exonic.gene.size=exonic.gene.sizes)
write.csv(ensembl2symbol(dtpm, gtf), file=file.path(odir, "tablecounts_tpm.csv"))
