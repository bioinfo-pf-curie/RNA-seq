#!/usr/bin/python
# -*- coding: utf-8 -*-
# The script is used to generate SAMPLE.statFile.txt file for RNA analysis.
# J. Brayet
# 10/2016

import os
import sys
import argparse

# creation du parse des arguments
parser = argparse.ArgumentParser(description="Create statFile.txt file by sample. Ex : ./statFile.py -i A468T16.stats.txt -o A468T16.statFile.txt -s A468T16 -b A468T16-15days")
 
# declaration et configuration des arguments
parser.add_argument('-i', '--inputFile', type=str, action="store", default="", help="flagstat file")
parser.add_argument('-o', '--outputFile', type=str, action="store", default="", help="output file")
parser.add_argument('-s', '--sampleIndentifier', type=str, action="store", default="", help="sample ID")
parser.add_argument('-b', '--biologicalIndentifier', type=str, action="store", default="", help="biological ID")

# dictionnaire des arguments
dargs = vars(parser.parse_args())

sampleID = dargs["sampleIndentifier"]
biosampleID = dargs["biologicalIndentifier"]

flagstatFile = open(dargs["inputFile"],"r")

for line in flagstatFile:
    if "total" in line:
        totalReadsQCpassed = line.split(" ")[0]
        totalReadsQCfailed = line.split(" ")[2]
    if "duplicates" in line:
        duplicates = line.split(" ")[0]
    if ("mapped" in line) and ("%" in line):
        line = line.replace(" : N/A", ":0.00%")
        mappedReadsQCpassed = line.split(" ")[0]
        mappedReadsQCfailed = line.split(" ")[2]
        mappedPercents = line.split(" ")[4]

flagstatFile.close()

outputFile = open(dargs["outputFile"],"w")

sep="\t"
colnames=("Sample identifier", "Biological identifier", "Number of total reads (QC-passed)", "Number of total reads (QC-failed)", "Number of mapped reads (QC-passed)", "Percent of mapped reads (QC-passed)", "Number of mapped reads (QC-failed)", "Percent of mapped reads (QC-failed)", "Number of duplicates", "Percent of duplicates")

outputFile.write(sep.join(colnames)+"\n")

mappedPercentQCpassed = mappedPercents.split(":")[0].replace("(", "")
mappedPercentQCfailed = mappedPercents.split(":")[1].replace(")", "").replace("\n","")
duplicatePercent = "{0:.2f}".format(float(duplicates)*100/(float(totalReadsQCpassed)+float(totalReadsQCfailed)))

values=(sampleID, biosampleID, totalReadsQCpassed, totalReadsQCfailed, mappedReadsQCpassed, mappedPercentQCpassed, mappedReadsQCfailed, mappedPercentQCfailed, duplicates, str(duplicatePercent)+"%")

outputFile.write(sep.join(values))
outputFile.close()

