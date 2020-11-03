#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re
from os import path

regexes = {
    'rnaseq': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'STAR': ['v_star.txt', r"(\S+)"],
    'HISAT2': ['v_hisat2.txt', r"version (\S+)"],
    'bowtie': ['v_bowtie.txt', r"version (\S+)"],
    'bowtie2': ['v_bowtie2.txt', r"version (\S+)"],
    'Picard': ['v_picard.txt', r"Version:(\S+)"],
    'Samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'featureCounts': ['v_featurecounts.txt', r"featureCounts v(\S+)"],
    'HTseqCounts': ['v_htseq.txt', r"version (\S+)"],
    'Preseq': ['v_preseq.txt', r"Version: (\S+)"],
    'RSeQC': ['v_rseqc.txt', r"infer_experiment.py ([\d\.]+)"],
    'R': ['v_R.txt', r"R version (\S+)"]
}
results = OrderedDict()
results['rnaseq'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['STAR'] = '<span style="color:#999999;\">N/A</span>' 
results['HISAT2'] = '<span style="color:#999999;\">N/A</span>' 
results['Picard'] = '<span style="color:#999999;\">N/A</span>'
results['Samtools'] = '<span style="color:#999999;\">N/A</span>'
results['featureCounts'] = '<span style="color:#999999;\">N/A</span>'
results['HTseqCounts'] = '<span style="color:#999999;\">N/A</span>'
results['Preseq'] = '<span style="color:#999999;\">N/A</span>'
results['RSeQC'] = '<span style="color:#999999;\">N/A</span>'
results['R'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    if path.isfile(v[0]):
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))

# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'Software Versions'
section_href: 'https://gitlab.curie.fr/rnaseq'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k,v))
print ("    </dl>")
