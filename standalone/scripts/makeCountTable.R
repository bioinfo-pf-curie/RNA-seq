#rm(list=ls())

#args <- commandArgs(TRUE)
#la <- length(args)
#if (la > 0){
#  for (i in 1:la)
#    eval(parse(text=args[[i]]))
#}

source("lib_expan_rnaseq.R")
require(rtracklayer)

## Note cannot replace the gene_id by the gene_name because some genes are duplicated (with the same gene_name)
ensembl2symbol <- function(x, gtf.in){
  dgtf <- rtracklayer::import(gtf.in)
  my_genes <- dgtf[dgtf$type == "gene"]
  mcols(my_genes) <- mcols(my_genes)[c("gene_id", "gene_type","gene_name")]
  m <- match(rownames(x), my_genes$gene_id)
  if(length(!is.na(m)) != dim(x)[1]){
    warning("Unable to convert ENSEMBL to SYMBOL gene names ! ")
    return(x)
  }else{
    rownames(x) <- paste(my_genes$gene_id[m], my_genes$gene_name[m], sep="|")
    return(x)
  }
}



if (count_tool == "STAR"){
  ## Load STAR data
    message("loading STAR gene counts ...")
    exprs.in <- list.files(path=input_path, pattern="ReadsPerGene.out.tab", full.names=TRUE, recursive=TRUE)
    prefix <- gsub("_norRNAReadsPerGene.out.tab$","",basename(exprs.in))
    counts.exprs <- lapply(exprs.in, read.csv, sep="\t", header=FALSE, row.names=1)
    if (stranded == "reverse"){
        counts.exprs <- data.frame(lapply(counts.exprs, "[", 3))
    }else if (stranded == "yes"){
        counts.exprs <- data.frame(lapply(counts.exprs, "[", 2))
    }else{
        counts.exprs <- data.frame(lapply(counts.exprs, "[", 1))
    }    
    colnames(counts.exprs) <- gsub("_norRNA", "", gsub("ReadsPerGene.out.tab","",sapply(exprs.in, basename)))
    ## remove first 4 lines
    counts.exprs <- counts.exprs[5:nrow(counts.exprs), , drop=FALSE]
}else if (count_tool == "FEATURECOUNTS"){
    ## Load FeatureCounts data
    exprs.in <- list.files(path=input_path, pattern="featurecounts.csv$", full.names=TRUE, recursive=TRUE)
    counts.exprs <- lapply(exprs.in, function(f){
        z <- read.csv(f, sep="\t", row.names=1, comment.char="#")
        data.frame(z[,6], row.names=rownames(z))
    })
    counts.exprs <- data.frame(counts.exprs)
    colnames(counts.exprs)<- gsub("_norRNA", "", gsub("featurecounts.csv","",sapply(exprs.in, basename)))
}else if (count_tool == "HTSEQCOUNT"){
    ## Load HTSeq data
    exprs.in <- list.files(path=input_path, pattern="htseq.csv", full.names=TRUE, recursive=TRUE)
    counts.exprs <- lapply(exprs.in, read.csv, sep="\t", header=FALSE, row.names=1)
    counts.exprs <- as.data.frame(counts.exprs)
    colnames(counts.exprs)<- gsub("_norRNA", "", gsub("htseq.csv","",sapply(exprs.in, basename)))
}

## export count table(s)
write.csv(ensembl2symbol(counts.exprs, gtf), file=file.path(odir, "tablecounts_raw.csv"))

## TPM
exonic.gene.sizes <- getExonicGeneSize(gtf)
dtpm <- getTPM(counts.exprs, exonic.gene.size=exonic.gene.sizes)
write.csv(ensembl2symbol(dtpm, gtf), file=file.path(odir, "tablecounts_tpm.csv"))
