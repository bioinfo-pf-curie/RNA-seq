#!/usr/bin/env Rscript

## Load libraries
library(proxy)
library(dendextend)
library(optparse)

option_list <- list(make_option(c("-i", "--input"), type="character", default=NULL, help="Input file with SNPs allelic frequencies", metavar="path"),
                    make_option(c("-s", "--splan"), type="character", default=NULL, help="Sample plan", metavar="path"),
                    make_option(c("-o", "--odir"), type="character", default="./", help="Output director", metavar="path"),
                    make_option(c("-d", "--dist"), type="character", default='ejaccard', help="Clustering distance", metavar="character"),
                    make_option(c("-m", "--minsnp"), type="integer", default='22', help="Minimum SNPs number to consider a sample", metavar="int"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#args<-commandArgs(trailingOnly = TRUE)
#if (length(args) < 2) {
#    stop("Usage: compute_clust.r <inputTable> <outputDir> <distance> [<minSNPsNumber>]", call.=FALSE)
#}

# Load arguments
inputTable <- opt$input
samplePlan <- opt$splan
outputDir <- opt$odir
distance <- opt$dist
minSNPsNumber <- opt$minsnp

# Handle path & Output Names
outputPlot <- paste0(outputDir,"/clustering_identito_mqc.png")
outputMatrix <- paste0(outputDir,"/clustering_identito.csv")

#Load data
d <- read.table(inputTable, header=TRUE, stringsAsFactors=FALSE, row.names=1)
d <- as.matrix(d)
d[which(d == "NEC")] <- NA

clust_mat <- matrix(as.numeric(d), ncol=ncol(d),
                    dimnames=list(rownames(d),colnames(d)))

# Update sample name with sample plan
if (!is.null(samplePlan)){
    splan <- read.csv(samplePlan, header=FALSE, row.names=1)
    if (length(setdiff(rownames(clust_mat), rownames(splan))) == 0){
        splan <- splan[rownames(clust_mat),]
        rownames(clust_mat) <- splan[,1]
    }
}

# Handle sample full of NA
toRemove <- c()
for (i in 1:nrow(as.matrix(clust_mat))){
    if (sum(!is.na(clust_mat[i,])) < minSNPsNumber){
        warning('One or more samples contain only NAs, those samples will be removed')
        toRemove <- c(toRemove, i)
    }
}

if (length(toRemove) > 0){
    clust_mat <- clust_mat[-toRemove,]
}

if (is.null(dim(clust_mat)) || nrow(clust_mat) == 0){
    warning("Not enough samples to run the clustering")
    #file.create(outputMatrix)
    #df <- data.frame()
    #write.csv(df, outputMatrix)
    quit(save="no", status = 0)
}

if (distance == "ejaccard"){
    dt <- as.matrix(dist(clust_mat, method="eJaccard"))
}else if (distance == "euclidean"){
    dt <- as.matrix(dist(clust_mat, method="Euclidean"))
}else{
    stop("Undefined distance metric")
}


## remove NA values in Jaccard distances
naToRemove <- which(apply(dt, 1, function(x){length(which(is.na(x)))}) >= 1)

if (length(naToRemove) > 0){
    dmat <- dt[-naToRemove, -naToRemove]
}else {
    dmat <- dt
}

if (!is.null(dim(dmat)) && dim(dmat)[1]>1){
    clust <- hclust(as.dist(dmat), "ward.D2")
    clust_order <- clust$order
    
    png(outputPlot, width=1400, height=600)
    par(mar=c(8,5,2,2))
    clust %>%
        as.dendrogram() %>%
        set("branches_col", "darkgrey") %>% set("branches_lwd", 2) %>%
        set("labels_cex", 1) %>%
        set("leaves_pch", 16)  %>% set("leaves_col", "orange") %>%
        plot()
    dev.off()
    
    write.csv(round(dmat[clust_order,clust_order],3), outputMatrix, quote=TRUE, row.names=TRUE)
}else{
    warning("Not enough samples to run the clustering")
    write.csv(round(jacdist,3), outputMatrix, quote=TRUE, row.names=TRUE)
}

print("script ended successfully")
