#!/usr/bin/env Rscript

args<-commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
    stop("Usage: apComputePolym.r <inputList> <output> <sample name><bed polym> <minDP (optional, default=10)>", call.=FALSE)
}

inputFiles <- args[1]
outputFile <- args[2]
sampleName <- args[3]
polym <- args[4]
minDP <- ifelse(is.na(args[5]),10,as.numeric(as.character(args[5])))

## load BED with polym
bed <- read.table(polym, stringsAsFactors=FALSE)[,c(1,3:4),]
bed <- cbind(bed,data.frame(Reduce(function(x,y) rbind(x,y),
                                   strsplit(gsub(">","_",as.character(bed[,"V4"])),"_"), accumulate=FALSE)))
colnames(bed) <- c("chr","pos","name","gene","rs","ref_polym","alt_polym")

## load tsv file from SnpSift with AD & DP:
var_tsv <- read.table(inputFiles, header=TRUE, stringsAsFactors=FALSE)
colnames(var_tsv) <- c("chr","pos","ref","alt","DP","AD")

## Merge polym & var_tsv
polym_table <- merge(bed, var_tsv, by=c("chr","pos"), all=TRUE)

VAF <- NULL
for (i in 1:nrow(polym_table)){
    polym_indice <- match(polym_table$alt_polym[i],
                          c(polym_table[i,"ref"],unlist(strsplit(as.character(polym_table$alt[i]),","))))
    polym_AD <- as.numeric(unlist(strsplit(as.character(polym_table[i,11]),","))[polym_indice])
    if (is.na(polym_AD)) {
        compute_vaf <- "NEC"
        VAF <- c(VAF,compute_vaf)
    } else {
        compute_vaf <- round(polym_AD/polym_table[i,"DP"],3)
        VAF <- c(VAF,compute_vaf)
    }
}

polym_table$VAF <- VAF

# Create matrix for clustering:
clust_mat <- as.data.frame(t(polym_table[,c(5,12)]))
rownames(clust_mat) <- c("rs",sampleName)
write.table(clust_mat, outputFile, quote=FALSE, col.names=FALSE, row.names=TRUE, sep="\t")
