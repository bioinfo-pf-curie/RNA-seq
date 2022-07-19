#!/usr/bin/env Rscript

stopifnot(require(rtracklayer))
stopifnot(require(edgeR))

# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Usage: gene_type_expression.r <counts_table> <gtf> <output_file> <R-package-location (optional)>", call.=FALSE)
}
counts <- args[1]
gtf <- args[2]
ofile <- args[3]
if (length(args) > 3) { .libPaths( c( args[4], .libPaths() ) ) }

message("Count table (Arg 1):", counts)
message("GTF file (Arg 2):", gtf)
message("Output file (Arg 3):", ofile)

stopifnot(require(rtracklayer))
if (grepl(".tsv$", counts)){
    d <- read.table(counts, row.names=1, check.names=FALSE)
}else if (grepl(".csv$", counts)){
    d <- read.csv(counts, row.names=1, check.names=FALSE)
}else{
    stop("Unexpected counts file format. '.csv' or '.tsv' are supported")
}

## count number of expressed genes
y <- DGEList(counts=d)
keep <- filterByExpr(y)
df <- d[keep,,drop=FALSE]

## gene annotation
d.gtf <- rtracklayer::import(gtf)
my_genes <- d.gtf[d.gtf$type == "gene"]

## remove "." in ENSEMBL Ids
my_genes$gene_id <- gsub("\\.[0-9]+$","",my_genes$gene_id)
rownames(df) <- gsub("\\.[0-9]+$","",rownames(df))

## Droso samples
if (!is.element("gene_type", colnames(elementMetadata(my_genes))) && is.element("gene_biotype", colnames(elementMetadata(my_genes)))){
    md <- elementMetadata(my_genes)
    cn <- colnames(md)
    cn[which(cn=="gene_biotype")] <- "gene_type"
    colnames(md) <- cn
    elementMetadata(my_genes) <- md
}

if (length(my_genes) > 0 && is.element("gene_type", colnames(elementMetadata(my_genes)))){
   mcols(my_genes) <- mcols(my_genes)[c("gene_id", "gene_type","gene_name")]
   n_items <- 5
   d2p <- as.matrix(data.frame(lapply(as.list(df), function(x){
          n_ex <- length(which(x>1))
   	  ids <- sapply(strsplit(rownames(df)[which(x>1)], "\\|"), "[[", 1)
	  dt <- table(factor(my_genes$gene_type[match(ids, my_genes$gene_id)], levels=unique(sort(my_genes$gene_type))))
	  c(total=n_ex, dt)
	}), check.names=FALSE))
}else{
    n_items <- 1
    d2p <-  as.matrix(data.frame(lapply(as.list(df), function(x){
            n_ex <- length(which(x>1))
            c(total=n_ex)
            }), check.names=FALSE))
}

## remove zero values and sort by counts
line2remove <- c(which(rownames(d2p)=="total"), which(rowSums(d2p)==0))
if (length(line2remove) > 0 ){
   d2p <- d2p[-line2remove,,drop=FALSE]
}
d2p <- d2p[order(rowSums(d2p), decreasing=TRUE),,drop=FALSE]

## reduce
if (nrow(d2p) > (n_items + 1)){
   d2p <- rbind(d2p[1:n_items,,drop=FALSE], others=colSums(d2p[(n_items+1):nrow(d2p),,drop=FALSE]))
}

## export
write.csv(t(d2p), file=ofile)
