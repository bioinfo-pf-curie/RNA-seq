#!/usr/bin/env Rscript

args<-commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: compute_clust.r <inputTable> <outputDir> <figureName>", call.=FALSE)
}

# Load arguments
inputTable <- args[1]
outputDir <- args[2]
minSNPsNumber <- ifelse(is.na(args[3]),10,as.numeric(as.character(args[3])))

# Handle path & Output Names
outputPlot <- paste0(outputDir,"/clustering_identito.png")
outputMatrix <- paste0(outputDir,"/clustering_identito.csv")

# Load libraries
library(pheatmap)
library(proxy)

# Saving function
save_pheatmap_png <- function(x, filename, width=1000, height=1000) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    png(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
}

#Load data
d <- read.table(inputTable, header=TRUE, stringsAsFactors=FALSE, row.names=1)
d <- as.matrix(d)
d[which(d == "NEC")] <- NA

clust_mat <- matrix(as.numeric(d), ncol=ncol(d),
                    dimnames=list(rownames(d),colnames(d)))

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
    file.create(outputMatrix)
    df <- data.frame()
    write.csv(df, outputMatrix)
    quit(save="no", status = 0)
}

jacdist <- as.matrix(dist(clust_mat, method="eJaccard"))

## remove NA values in Jaccard distances
naToRemove <- which(apply(jacdist, 1, function(x){length(which(is.na(x)))}) > 1)
jacmat <- jacdist[-naToRemove, -naToRemove]

if (!is.null(dim(jacmat)) && dim(jacmat)[1]>1){
    clust <- pheatmap(jacmat,
                      method = "ward.D2",
                      clustering_distance_rows=as.dist(jacdist),
                      clustering_distance_cols=as.dist(jacdist),
                      show_rownames = TRUE,
                      show_colnames = TRUE,
                      display_numbers = TRUE,
                      legend = TRUE ,
                      fontsize = 20,
                      silent = TRUE)
    
    clust_order <- clust$tree_col$order
    save_pheatmap_png(clust, outputPlot)
    write.csv(round(jacmat[clust_order,clust_order],3), outputMatrix, quote=TRUE, row.names=TRUE)
}else{
    warning("Not enough samples to run the clustering")
    write.csv(round(jacdist,3), outputMatrix, quote=TRUE, row.names=TRUE)
}

# End message + Path:
print("script ended successfully")
