#!/usr/bin/env Rscript
# Command line argument processing
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("Usage: makeCountTable.r <inputList> <gtf> <count_tool> <stranded>", call.=FALSE)
}
inputFiles <- args[1]
gtf <- args[2]
count_tool <- toupper(args[3])
stranded <- args[4]

if (!require(rtracklayer)){
  source("http://bioconductor.org/biocLite.R")
  biocLite("rtracklayer", suppressUpdates=TRUE)
  library("rtracklayer") 
}

if (!require(GenomicFeatures)){
  source("http://bioconductor.org/biocLite.R")
  biocLite("GenomicFeatures", suppressUpdates=TRUE)
  library("GenomicFeatures")
}

##ensembl2symbol
## x : a count matrix with gene name as rownames
## gtf.in : GTF file used for annotation with ENSEMBL as gene_id and SYMBOL as gene_name (i.e GENCODE GTF)

ensembl2symbol <- function(x, gtf.in, lab.in="gene_id", lab.out=c("gene_name", "gene_symbol")){
  dgtf <- rtracklayer::import(gtf.in)
  myGenes <- dgtf[dgtf$type == "gene"]
  
  ## Check input label (gene_id by default)
  if (is.element(lab.in, colnames(elementMetadata(myGenes)))){
    colsOfInterest <- lab.in
  }else{
    warning("Unable to convert ID to SYMBOL gene names ! Input label not found !")
    return(x)
  }
  
  ## Try all output labels
  for (lab in lab.out){
    if (is.element(lab, colnames(elementMetadata(myGenes)))){
      colsOfInterest <- c(colsOfInterest, lab)
      break
    }
  }
  if(length(colsOfInterest) != 2){
    warning("Unable to convert ID to SYMBOL gene names ! Output label not found !")
    return(x)
  }
  
  mcols(myGenes) <- mcols(myGenes)[colsOfInterest]
  m <- match(rownames(x),  elementMetadata(myGenes)[,lab.in])
  if(length(!is.na(m)) != dim(x)[1]){
    warning("Unable to convert ENSEMBL to SYMBOL gene names ! ")
    return(x)
  }else{
    rownames(x) <- paste(elementMetadata(myGenes)[,lab.in][m], elementMetadata(myGenes)[,colsOfInterest[2]][m], sep="|")
    return(x)
  }
}


## getTPM
## Calculate TPM values from a gene expression matrix and a gtf file
## x : matrix of counts
## exonic.gene.size : vector of exonic sizes per gene. If not provided, gtf.in is used to build it
## gtf.in : path to gtf file
## Details : see http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/ for differences between RPKM and TPM

getTPM <- function(x, exonic.gene.sizes=NULL, gtf.in=NULL){

  stopifnot(require(GenomicFeatures))
  
  ## First, import the GTF-file that you have also used as input for htseq-count
  if (is.null(exonic.gene.sizes)){
    if (is.null(gtf.in))
      stop("Unable to calculate gene size")
    exonic.gene.sizes <- getExonicGeneSize(gtf.in)
  }
  
  if (length(setdiff(rownames(x), names(exonic.gene.sizes)))>0){
    warning(length(setdiff(rownames(x), names(exonic.gene.sizes))), " genes from table were not found in the gtf file ...")
  }
  
  exonic.gene.sizes <- exonic.gene.sizes[rownames(x)]
  
  ## Calculate read per kilobase
  rpk <- x * 10^3/matrix(exonic.gene.sizes, nrow=nrow(x), ncol=ncol(x), byrow=FALSE)
  ## Then normalize by lib size
  tpm <- rpk *  matrix(10^6 / colSums(rpk, na.rm=TRUE), nrow=nrow(rpk), ncol=ncol(rpk), byrow=TRUE)
  
  return(round(tpm,2))
}##getTPM


## getExonicGeneSize
## Calcule the exons size per gene for RPMKM/TPM normalization
## gtf.in : path to gtf file
getExonicGeneSize <- function(gtf.in){
  stopifnot(require(GenomicFeatures))
  txdb <- makeTxDbFromGFF(gtf.in,format="gtf")
  ## then collect the exons per gene id
  exons.list.per.gene <- exonsBy(txdb,by="gene")
  ## then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
  exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))
  return(exonic.gene.sizes)
}## getExonicGeneSize


##################################################################

exprs.in <- as.vector(read.table(inputFiles, header=FALSE)[,1])

if (count_tool == "STAR"){
  ## Load STAR data
  message("loading STAR gene counts ...")
  counts.exprs <- lapply(exprs.in, read.csv, sep="\t", header=FALSE, row.names=1, check.names=FALSE)
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
  counts.exprs <- lapply(exprs.in, function(f){
    z <- read.csv(f, sep="\t", row.names=1, comment.char="#", check.names=FALSE)
    data.frame(z[,6], row.names=rownames(z))
  })
  counts.exprs <- data.frame(counts.exprs)
  colnames(counts.exprs)<- gsub("_norRNA", "", gsub("_featurecounts.csv","",sapply(exprs.in, basename)))
}else if (count_tool == "HTSEQCOUNT"){
  ## Load HTSeq data
  counts.exprs <- lapply(exprs.in, read.csv, sep="\t", header=FALSE, row.names=1, check.names=FALSE)
  counts.exprs <- as.data.frame(counts.exprs)
  colnames(counts.exprs) <- gsub("_norRNA", "", gsub("_htseq.csv","",sapply(exprs.in, basename)))
}

## export count table(s)
write.csv(ensembl2symbol(counts.exprs, gtf), file="tablecounts_raw.csv")

## TPM
exonic.gene.sizes <- getExonicGeneSize(gtf)
dtpm <- getTPM(counts.exprs, exonic.gene.size=exonic.gene.sizes)
write.csv(ensembl2symbol(dtpm, gtf), file="tablecounts_tpm.csv")
