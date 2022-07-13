#!/usr/bin/env Rscript

args<-commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: apComputePolym.r <inputList> <output> <sample name><bed polym> <minDP (optional, default=10)>", call.=FALSE)
}

inputFiles <- args[1]
outputFile<-args[2]
sampleName<-args[3]
polym<-args[4]
minDP<-ifelse(is.na(args[5]),10,as.numeric(as.character(args[5])))

## load BED with polym
bed <- read.table(polym,stringsAsFactors = FALSE)[,c(1,3:4),]
bed <- cbind(bed, data.frame(Reduce(function(x,y) rbind(x,y),
                                    strsplit(gsub(">","_",as.character(bed[,"V4"])),"_"), accumulate=FALSE)))
colnames(bed) <- c("chr","pos","name","gene","rs","ref_polym","alt_polym")
rownames(bed) <- bed$rs

## load tsv file from SnpSift with AD & DP:
var_tsv <- read.table(inputFiles, header=TRUE, stringsAsFactors = FALSE)
colnames(var_tsv) <- c("chr","pos","ref","alt","DP","AD")

## Merge polym & var_tsv
polym_table <- merge(bed, var_tsv, by=c("chr","pos"), all=TRUE)
rownames(polym_table) <- polym_table$rs

### reorder polym table
polym_table <- polym_table[rownames(bed),]

VAF <- rep(NA, nrow(polym_table))
for (i in 1:nrow(polym_table)){
    if (is.na(polym_table[i,"DP"]) || polym_table[i,"DP"] < minDP){
        VAF[i] <- "NEC"
    }else{
        polym_indice <- match(polym_table[i, "alt_polym"], c(polym_table[i,"ref"], unlist(strsplit(as.character(polym_table[i, "alt"]),","))))
        ## No alternative allele is detected
        if (is.na(polym_indice)){
            VAF[i] <- "0"
        }else{
            polym_AD <- as.numeric(unlist(strsplit(as.character(polym_table[i,"AD"]),","))[polym_indice])
            compute_vaf <- round(polym_AD/polym_table[i,"DP"]*100,3)
            VAF[i] <- compute_vaf
        }
    }
}
polym_table$VAF <- VAF

# Create matrix for clustering:
clust_mat <- as.data.frame(t(polym_table[,c(5,12)]))
rownames(clust_mat) <- c(" ", sampleName)
write.table(clust_mat, outputFile, quote=FALSE, col.names=FALSE, row.names=TRUE, sep="\t")
