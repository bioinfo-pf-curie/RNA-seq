##
## Library for exploratory analysis of RNA-seq data
## To use it ; source("exploratory_analysis.R")
##



## plotHeatCor
## Plot a correlation heatmap of all samples
## x : matrix with gene in rows and sample in columns
## cor.met : correlation metric, either spearman or pearson
## add.text : if true, add the correlation values
## col.range : color range for the heatmap

plotHeatCor <- function(x, cor.met="spearman", add.text=TRUE, col.range=c(-1, 1)){
  require(reshape2)
  require(ggplot2)
  
  cormat <- round(cor(x, method=cor.met),3)
  cormat[lower.tri(cormat)]<- NA

  ## Melt the correlation matrix
  meltedcor <- melt(cormat, na.rm = TRUE)
  
  ## Create a ggheatmap
  ggheatmap <- ggplot(meltedcor, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white") +
      scale_fill_gradient2(low = "blue", high = "red",
                           mid = "white", midpoint = mean(col.range), limit = col.range, space = "Lab",
                           name="Pearson\nCorrelation") +
                             theme_minimal() + # minimal theme
                               theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                       size = 12, hjust = 1)) + coord_fixed()
  if (add.text){
    ggheatmap <- ggheatmap + geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            legend.justification = c(1, 0),
            legend.position = c(0.6, 0.7),
            legend.direction = "horizontal")+
              guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                       title.position = "top", title.hjust = 0.5))
  }
  print(ggheatmap)
}##plotHeatCor



## plotPCA
## Perform a PCA analysis on the gene expression matrix
## x : matrix with gene in rows and sample in columns
## splan : data.frame describing the sample plan. Must have a condition column

plotPca <- function(x, splan){

  require(DESeq2)
  
  dds <- DESeqDataSetFromMatrix(countData=x, DataFrame(condition=splan$condition), ~ condition)
  dds <- estimateSizeFactors(dds)
  rld <- rlog(dds, blind=TRUE)
  
  res.pca <- plotPCA(rld, intgroup = c("condition"), returnData=TRUE, ntop=1000)
  percentVar <- round(100 * attr(res.pca, "percentVar"))
  ggplot(res.pca, aes(PC1, PC2, color=condition)) + geom_point(size=3) +
    geom_text(aes(label=name)) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance"))
}##plotPCA



## plotClutering
## Perform a Hierachical clutering (Pearson/Ward) on the gene expression matrix
## x : matrix with gene in rows and sample in columns
## splan : data.frame describing the sample plan. Must have a condition column

plotClustering <- function(x, splan){
  
  require("pheatmap")
  require("RColorBrewer")

  dds <- DESeqDataSetFromMatrix(countData=x, DataFrame(condition=splan$condition), ~ condition)
  dds <- estimateSizeFactors(dds)
  rld <- rlog(dds, blind=TRUE)
    
  sampleDists <- as.dist(1 - abs(cor(assay(rld), method = "spearman")))
  cl <- hclust(sampleDists, "ward.D2")
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(colnames(rld), rld$condition, sep="-" )
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           clustering_method="ward.D2",
           col=colors)
}##plotClustering


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



