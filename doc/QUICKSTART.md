---
title: "MANUAL"
author: "Nicolas Servant"
date: "7 septembre 2016"
output: html_document
---

# RNA-seq pipeline
# v0.1.0

## Quick start guide

This pipeline is a bioinformatics worflow to process RNA sequencing reads regardless the protocol used (mRNA or total RNA).  
The current version is gene-based oriented and do not include any transcript-based analysis.  
It is designed to take the raw sequence files (fastq) from a set of experiments (a sequencing run) or a single sample, and to process the data up to the final read counts per gene.  
A report with RNA-seq quality metrics is also generated.

### Main features

The RNA-seq pipeline is based on the following workflow :

![RNA-seq worflow][RNApip-wkf.png]

### Dependancies

The current version is designed to run on the Centos cluster.  
By default, the current dependencies (and paths) are used :

Tools | Path on cluster 
--- | --- 
FASTQC_PATH | /bioinfo/local/build/Centos/FastQC/FastQC_v0.11.5
JAVA_PATH | /bioinfo/local/build/Centos/java/jre1.8.0_101/bin
BOWTIE_PATH | /bioinfo/local/build/Centos/bowtie/bowtie-1.2/bin
TOPHAT2_PATH | /bioinfo/local/build/Centos/tophat/tophat_2.1.1/bin
BOWTIE2_PATH | /bioinfo/local/build/Centos/bowtie2/bowtie2-2.2.9
STAR_PATH | /bioinfo/local/build/Centos/STAR/STAR-2.5.3a/bin/Linux_x86_64
SAMTOOLS_PATH | /bioinfo/local/build/Centos/samtools/samtools-1.3/bin
HTSEQ_PATH | /bioinfo/local/build/Centos/python/python-2.7.13/bin
FEATURECOUNTS_PATH | /bioinfo/local/build/Centos/subread/subread-1.5.1-Linux-x86_64/bin
BEDTOOLS_PATH | /bioinfo/local/build/Centos/bedtools/bedtools-2.25.0/bin
R_PATH | /bioinfo/local/build/Centos/R/R-3.4.0/bin
PYTHON_PATH | /bioinfo/local/build/Centos/python/python-2.7.11/bin
RSEQC_PATH | /bioinfo/local/build/Centos/RSeQC/RSeQC-2.6.4/scripts
PICARD_PATH | /bioinfo/local/build/Centos/picard/2.6.0/picard.jar
PRESEQ_PATH | /bioinfo/local/build/Centos/preseq/preseq_v2.0/


### Setting the configuration files

All options used during the processing can be defined in a CONFIGURATION file.
In this case, simply copy the CONFIG file from the installation folder to your local folder, and edit it.

PARAMETERS | DESCRIPTION 
ORG | Organism
BUILD | Version of reference genome
RUN_FASTQC | Run Fastq for quality controls; 0=no, 1=yes (default: no)
MAPPING_TOOL | Mapping tool. Must be STAR or TOPHAT2 (default: STAR)
COUNT_TOOL | Tools to quantify expression per gene. Must be STAR, HTSEQ or FEATURECOUNTS (default: STAR)
STRANDED | The strandness of the protocol. Must be yes, no or reverse following HTSeq-Count standard. If not specified, the pipeline will detect it automatically. (default: ) 
BOWTIE_RRNA_IDX | Indexes for rRNA mapping. If not specified, this step is skipped
BOWTIE_OPTS | Option for bowtie mapping on rRNA reference (default: -v 2 -a -m 1 --best --strata --nomaqround -y)
STAR_IDX_PATH | Path to STAR indexes
STAR_OPTS | Options for STAR mapping (default: --runThreadN 8 --outSAMtype BAM SortedByCoordinate --runMode alignReads --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMprimaryFlag OneBestScore --outMultimapperOrder Random --outSAMattributes All)
TOPHAT2_IDX_PATH | Path to Tophat2 indexes
TOPHAT2_OPTS | Options for TopHat2 mapping (default: --b2-sensitive -p 6 -g 1 -N 2 --no-coverage-search)
TRANSCRIPTS_GTF | GTF file for annotation
HTSEQ_OPTS | Options for HTSeq read counts (default: -f bam -t exon -r pos) 
FEATURECOUNTS_OPTS | Options for FeatureCounts read counts (default: -t exon -C -T 8 -p)
BWT2_IDX_PATH | Path to bowtie2 indexes for RSeQC usage
GENE_BED | Path to BED file for RSeQC usage


### Run the pipeline on a single sample :

The RNA-seq pipeline is currently available at **/bioinfo/local/curie/ngs-data-analysis/centos/RNAseq_pip/** (hereafter refer as PIP_PATH)
 
```
PIP_PATH/bin/RNApip -f SAMPLE_1.fastq.gz [-r SAMPLE_2.fastq.gz] -o PIPELINE_OUTPUT_DIR -c CONFIG [-s RUN_ID] 
```

### Run the pipeline on a list of sample, using a cluster :

In cases of multiple samples, the user has to write a SAMPLE_PLAN file with the following information
SAMPLE_ID | SAMPLE_NAME | PATH_TO_R1.fastq | PATH_TO_R2.fastq   
 
Using the *RNApip_cluster* script, all samples will be processed in parallel:  

```
PIP_PATH/bin/RNApip_cluster -i SAMPLE_PLAN -o PIPELINE_OUTPUT_DIR -c CONFIG 
```

### Run the report manually

In cluster mode, the HTML report is generated automatically.  
If for any reason, the user wants to manually generate the report, the following command can be used :

```
PIP_PATH/bin/makeRNAreport -i PIPELINE_OUTPUT_DIR -c CONFIG
```
