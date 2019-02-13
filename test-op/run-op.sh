#!/bin/bash

## Small dataset - Gendrel et al. Mouse
../bin/RNApip_cluster -i /bioinfo/local/curie/ngs-data-analysis/testdataset_pipelines/RNA-seq/SAMPLE_PLAN -o /data/tmp/nservant/RNA_TESTDATA \
    -c /bioinfo/local/curie/ngs-data-analysis/testdataset_pipelines/RNA-seq/CONFIG -k

## Test dataset - polyA RNA
../bin/RNApip_cluster -i SAMPLE_PLAN_B249 -o /data/tmp/RNAPIP_OP_B249/ -c config-files/CONFIG_Human.hg19 -n testop_B249 -k

## Test dataset - total RNA
../bin/RNApip_cluster -i SAMPLE_PLAN_L81 -o /data/tmp/RNAPIP_OP_L81/ -c config-files/CONFIG_Human.hg38 -n testop_L81 -k

## Test dataset - Human Tophat - featureCounts
../bin/RNApip_cluster -i SAMPLE_PLAN_G448 -o /data/tmp/RNAPIP_OP_G448/ -c CONFIG_G448 -n testop_G448 -k

## Test dataset - Human hg38
../bin/RNApip_cluster -i SAMPLE_PLAN_D095 -o /data/tmp/RNAPIP_OP_D095/ -c CONFIG_D095 -n testop_D095 -k
