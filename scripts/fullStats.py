#!/usr/bin/env python
# The script is used to generate SAMPLE.statFile.txt file for RNA analysis.
# J. Brayet - 10/2016
# Modified N. Servant - 29-05-17

import os
import sys
import argparse

from os import listdir
from os.path import isdir, isfile, join

# creation du parse des arguments
parser = argparse.ArgumentParser(description="Combine statFile.txt file for several samples")
 
# declaration et configuration des arguments
parser.add_argument('-i', '--inputs', type=str, nargs="+", help="Input stats files")
parser.add_argument('-o', '--odir', type=str, help="Output directory")
parser.add_argument('-r', '--runID', type=str, help="runID")

# dictionnaire des arguments
dargs = vars(parser.parse_args())

############## Function - Create xls file #############

def createXLSfile(txtFile, xlsFile):

    import xlwt

    book = xlwt.Workbook()
    ws = book.add_sheet('Stat file')

    inputFile = open(txtFile, "r+")

    data = inputFile.readlines() # read all lines at once
    for i in range(len(data)):
        row = data[i].split(",")  # This will return a line of string data, you may need to convert to other formats depending on your use case
        for j in range(len(row)):
            ws.write(i, j, row[j].replace("\n",""))

    book.save(xlsFile)
    inputFile.close() 

#######################################################

outputPath = dargs["odir"]
samples = dargs["inputs"]
runID = dargs["runID"]

outputFile = open(outputPath + "/fullStatFile.txt", "w")

write_header=0
for sample in samples:
    print sample
    statFileBySample = open(sample, "r")
    for i, line in enumerate(statFileBySample):
        if i == 0 and write_header == 0:
            outputFile.write(line) 
            write_header=1
        if i == 1:
            outputFile.write(line) 
    statFileBySample.close()
outputFile.close()

createXLSfile(outputPath+"/fullStatFile.txt", outputPath+"/fullStatFile.xls")




