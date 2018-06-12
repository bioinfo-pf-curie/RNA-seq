##
## Library for RNA-seq pipeline
## To use it ; source("exploratory_analysis.R")
##

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


## getExonicGeneSize
## Calcule the exons size per gene for RPMKM/TPM normalization
## gtf.in : path to gtf file

getGeneCoordinatesFromGTF <- function(gtf.in, ids=c("gene_id", "gene_name", "gene_type")){
    ##stopifnot(require(GenomicFeatures))
    ##txdb <- makeTxDbFromGFF(gtf.in,format="gtf")
    ##genes(txdb)
    stopifnot(require(rtracklayer))
    geneannot<-import(gtf.in)
    geneannot <- geneannot[which(geneannot$type=="gene")]
    if (length(ids)>0)
        geneannot <- geneannot[,ids]
    sort(geneannot)
}##getGeneCoordinates

## getTPM
## Calculate TPM values from a gene expression matrix and a gtf file
## x : matrix of counts
## exonic.gene.size : vector of exonic sizes per gene. If not provided, gtf.in is used to build it
## gtf.in : path to gtf file
## Details : see http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/ for differences between RPKM and TPM

getTPM <- function(x, exonic.gene.sizes=NULL, gtf.in=NULL){

  ## First, import the GTF-file that you have also used as input for htseq-count
  stopifnot(require(GenomicFeatures))

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
  tpm <- rpk *  matrix(10^6 / colSums(rpk), nrow=nrow(rpk), ncol=ncol(rpk), byrow=TRUE)

  return(round(tpm,2))
}##getTPM


## getRPKM
## Calculate RPKM values from a gene expression matrix and a gtf file
## x : matrix of counts
## exonic.gene.size : vector of exonic sizes per gene. If not provided, gtf.in is used to build it
## gtf.in : path to gtf file
getRPKM <- function(x, exonic.gene.size=NULL, gtf.in=NULL){

  ## First, import the GTF-file that you have also used as input for htseq-count
  stopifnot(require(GenomicFeatures))

  if (is.null(exonic.gene.size)){
      if (is.null(gtf.in))
          stop("Unable to calculate gene size")
      exonic.gene.sizes <- getExonicGeneSize(gtf.in)
  }
  
  if (length(setdiff(rownames(x), names(exonic.gene.sizes)))>0){
      warning(ength(setdiff(rownames(x), names(exonic.gene.sizes))), " genes from table were not found in the gtf file ...")
  }
  
  exonic.gene.sizes <- exonic.gene.sizes[rownames(x)]

  ## Normalize by read size and by gene length
  rpkm <- x /
    matrix(colSums(x), nrow=nrow(x), ncol=ncol(x), byrow=TRUE)  /
      matrix(exonic.gene.sizes, nrow=nrow(x), ncol=ncol(x), byrow=FALSE) * 10^9

  ## or
  ## require(edgeR)
  ## rpkm <- rpkm(x, exonic.gene.sizes[rownames(x)])

  return(round(rpkm,2))
}##getRPKM


## estimate_saturation
## Estimate saturation of genes based on rarefaction of reads
## counts : a matrix of counts
## max_reads : maximum number of reads to downsample
## ndepths : resampling levels
## nreps : number of times the subsampling is performed to calculate a variance
## mincounts : minimum counts level to consider a gene as expressed. If NA, the threshold is fixed to 1 CPM
## extend.lines : If TRUE, the max number of detected genes is returned when the maximum number of reads is reached.
estimate_saturation <- function(counts, max_reads=Inf, ndepths=6, nreps=1, mincounts=NA, extend.lines=FALSE){
  stopifnot(require(S4Vectors))
  
  counts <- as.matrix(counts)
  readsums <- colSums(counts)
  max_reads <- min(max(readsums), max_reads)
  depths <- round(seq(from=0, to=max_reads, length.out=ndepths+1))
  
  saturation <- mclapply(1:ncol(counts), function(k){
    message("Processing sample ", colnames(counts)[k], "...")
    x <- counts[,k]
    nreads <- sum(x)
    
    ## minimum expression levels
    if (is.na(mincounts)){
      mincounts <- nreads / 1e6
    }
    ## max number of detected genes
    ngenes_detected <- length(which(x>=mincounts))
    
    probs <- x / nreads ## calculate gene probabilities for the library
    probs <- probs[probs > 0] ## zero counts add nothing but computational time!
    ngenes <- length(probs)
    res <- lapply(depths, function(dp, nreps=1, ...){
      rsim <- c(NA, NA)
      if (extend.lines)
        rsim <- c(ngenes_detected, NA)
      if (dp <= nreads){
        estim <- lapply(1:nreps, function(i, ngenes, dp, probs, mincounts){
          csim <- sample(x=ngenes, size=dp, replace=TRUE, prob=probs)
          length(which(runLength(Rle(sort(csim)))>=mincounts))
        }, ngenes=ngenes, dp=dp, probs=probs, mincounts=mincounts)
        rsim <- c(mean(unlist(estim)), var(unlist(estim)))
      }
      return(rsim)
    }, nreps=nreps, nreads=nreads,  ngenes=ngenes, probs=probs, mincounts=mincounts)
    data.frame(depths=depths,
               sat.estimates=sapply(res, "[[", 1),
               sat.var.estimates=sapply(res, "[[", 2))
  })
  names(saturation) <- colnames(counts)
  return(saturation)
}##estimate_saturation



##ensembl2symbol
## x : a count matrix with gene name as rownames
## gtf.in : GTF file used for annotation with ENSEMBL as gene_id and SYMBOL as gene_name (i.e GENCODE GTF)

ensembl2symbol <- function(x, gtf.in, lab.in="gene_id", lab.out=c("gene_name", "gene_symbol")){
  dgtf <- rtracklayer::import(gtf.in)
  my_genes <- dgtf[dgtf$type == "gene"]

  ## Check input label (gene_id by default)
  if (is.element(lab.in, colnames(elementMetadata(my_genes)))){
    colsOfInterest <- lab.in
  }else{
    warning("Unable to convert ID to SYMBOL gene names ! Input label not found !")
    return(x)
  }

  ## Try all output labels
  for (lab in lab.out){
    if (is.element(lab, colnames(elementMetadata(my_genes)))){
      colsOfInterest <- c(colsOfInterest, lab)
      break
    }
  }
  if(length(colsOfInterest) != 2){
    warning("Unable to convert ID to SYMBOL gene names ! Output label not found !")
    return(x)
  }

  mcols(my_genes) <- mcols(my_genes)[colsOfInterest]
  m <- match(rownames(x),  elementMetadata(my_genes)[,lab.in])
  if(length(!is.na(m)) != dim(x)[1]){
    warning("Unable to convert ENSEMBL to SYMBOL gene names ! ")
    return(x)
  }else{
    rownames(x) <- paste(elementMetadata(my_genes)[,lab.in][m], elementMetadata(my_genes)[,colsOfInterest[2]][m], sep="|")
    return(x)
  }
}

