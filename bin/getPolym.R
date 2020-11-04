#!/usr/bin/env Rscript

args<-commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: getPolym.r <inputList> <output> <bed polym> <minDP (optional, default=10)>", call.=FALSE)
}

inputFiles <- args[1]
outputFile<-args[2]
polym<-args[3]
minDP<-ifelse(is.na(args[4]),10,as.numeric(as.character(args[4])))

## load BED with polym
bed <- read.table(polym)[,c(1,3:4)]
bed <- cbind(bed,data.frame(Reduce(function(x,y) rbind(x,y), strsplit(gsub(">","_",as.character(bed[,4])),"_"), accumulate=FALSE)))
colnames(bed) <- c("chr","pos","name","gene","rs","ref","alt")

## look for files
lfile <- as.vector(read.table(inputFiles, header=FALSE)[,1])
res <- lapply(1:length(lfile),function(i){
    file <- lfile[i]
    name <- gsub("polym.g.vcf$","", basename(file))
    
    ## Load data
    tab <- read.table(file)[,c(1,2,4,5,9,10)]
    colnames(tab)=c("chr","pos","REF","ALT","FORMAT","GENO")
    tab <- merge(bed, tab, by=c("chr","pos"), all=TRUE)
    
    tab[,"idxAD"] <- apply(tab,1,function(x){
        xs <- unlist(strsplit(x["FORMAT"],":"))
        ifelse(length(which(xs=="AD"))>0, which(xs=="AD"), 0)})
    tab[,"idxDP"] <- apply(tab,1,function(x){
        xs <- unlist(strsplit(x["FORMAT"],":"))
        ifelse(length(which(xs=="DP"))>0, which(xs=="DP"), 0)})
    tab[,"idxALT"] <- apply(tab,1,function(x){ 
        lalt <- unlist(strsplit(as.character(unlist(x["ALT"])),","))
        v <- 0
        if(length(lalt)>1){
            v=which(lalt==x["alt"])
        }
        return(v)
    })
    
    #if homozygous ref : no AD then DP in 2nd position
    temp <- apply(tab,1,function(x){
     	if(x["idxALT"]!=0){
            #ad <- as.numeric(as.character(sapply(strsplit(sapply(strsplit(as.character(x["GENO"]),":"),"[",2),","),"[",as.numeric(as.character(x["idxALT"]))+1)))
            adall <- unlist(strsplit(x["GENO"],":"))[as.numeric(x["idxAD"])]
            ad <- as.numeric(as.character(unlist(strsplit(adall, ","))[as.numeric(x["idxALT"])+1]))
            dp <- as.numeric(as.character(sapply(strsplit(as.character(x["GENO"]),":"),"[",as.numeric(x["idxDP"]))))
     	}else{
            ad <- NA 
            #dp <- as.numeric(as.character(sapply(strsplit(as.character(x["GENO"]),":"),"[",2)))
            dp <- as.numeric(as.character(sapply(strsplit(as.character(x["GENO"]),":"),"[",as.numeric(x["idxDP"]))))
     	}
     	return(list(ad,dp))
    })

    tab[,"AD"] <- as.numeric(as.character(sapply(temp,"[",1)))
    tab[,"DP"] <- as.numeric(as.character(sapply(temp,"[",2)))
    tab[,name] <- round(tab[,"AD"]*100/tab[,"DP"],2)
    
    if(length(which(is.na(tab[,name]) & tab[,"DP"]>=minDP))>0){
        tab[which(is.na(tab[,name]) & tab[,"DP"]>=dp_min),name] <- 0
    }
    
    if(length(which(is.na(tab[,name]) & tab[,"DP"]<minDP))>0){
        tab[which(is.na(tab[,name]) & tab[,"DP"]<dp_min),name] <- paste0("NEC:",tab[which(is.na(tab[,name]) & tab[,"DP"]<dp_min),"DP"])
    }
    tab[,"Sample_ID"] <- paste(tab[,"gene"],tab[,"rs"],sep="_")
    #return(tab[,c("Sample_ID",name)])
    return(tab)
})

dres <- Reduce(function(x,y) merge(x,y,by="Sample_ID",all=TRUE), res, accumulate=FALSE)
write.table(data.frame(t(dres)), outputFile, quote=FALSE, col.names=FALSE, row.names=TRUE, sep="\t")
