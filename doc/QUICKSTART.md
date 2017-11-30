---
title: "MANUAL"
author: "Nicolas Servant"
date: "7 septembre 2016"
output: html_document
---

# RNA-seq pipeline

## Quick start guide

This pipeline is a bioinformatics worflow to process RNA sequencing reads regardless the protocol used (mRNA or total RNA).
The current version is gene-based oriented and do not include any transcript based analysis.
It is designed to take the raw sequence ouput from a set of experiments (a sequencing run) or a single sample, and to processed the data up to the final read counts per gene.
It will also produce a set of metrics and a report that can be used to assess the quality of the data.

### Main features

* Quality control
* RNA-seq read mapping
* Gene-based quantification
* Annotation

### Dependancies

The following dependencies are required :
* Unix-based operating system
* Tophat2
* STAR
* HTSeqCount
* samtools
* python2.XX
* FeatureCount
* R


### Run the pipeline on a single sample :

Run the analysis worflow  
 ./bin/RNApip -f /data/tmp/testdataset_pipelines/RNA-seq/Gendrel2012/SRR1106775_1M_1.fastq.gz \
	-r /data/tmp/testdataset_pipelines/RNA-seq/Gendrel2012/SRR1106775_1M_2.fastq.gz -o /data/tmp/testrnapip$$ -c CONFIG_star -s TEST_SAMPLE

Generate the Markdown report  
./bin/makeRNAreport -i /data/tmp/test_star -c CONFIG

### Run the pipeline on a list of sample, using a cluster :

./bin/RNApip_cluster -i /data/tmp/testdataset_pipelines/RNA-seq/SAMPLE_PLAN -o /data/tmp/test_star_p -c CONFIG_Mouse 
