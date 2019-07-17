# Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [rRNA Mapping](#rrna-mapping) - Alignment on ribosomal RNAs
* [Genome Mapping](#genome-mapping) - Alignment on the reference genome
* [Strandness](#strandness) - Sequencing orientation
* [Read Distribution](#read-distribution) - Distribution of sequencing reads
* [Complexity Curves](#complexity-curves) - Complexity of the librairies
* [Gene-based Saturation](#gene-based-saturation) - Libraries complexity based on number of detected genes
* [Counts](#counts) - Counts reads per gene
* [Expressed genes](#expressed-genes) - Number of expressed genes per sample
* [DupRadar](#dupradar) - Reads duplication level 
* [PCA](#pca) - Principal Component Analysis
* [Correlation](#correlation) - Sample Pearson's correlation
* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline


## FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows the input reads. In theory, they should be already trimmed for adapter sequence and potentially regions with low quality. 
For details about reads trimming, see the `raw_qc` pipeline.

**Output directory: `results/fastqc`**

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

## rRNA Mapping
When the annotation are available, the pipeline first aligns the sequencing reads on a database of ribosomal RNA sequence with the bowtie1 software.
This metric is mainly useful for protocols based on rRNA depletion. Samples with a high level of rRNA contamination should be carefully considered.

For detail about the bowtie1 mapper, see the [bowtie help](http://bowtie-bio.sourceforge.net/index.shtml)

> **NB:** Note that the fastq files after rRNA cleaning are exported only if the `--saveAlignedIntermediates` parameter is turn on.

**Output directory: `results/rRNA_mapping`**

* `logs/sample.log`
  * Log file of bowtie1 mapping with the number of aligned reads on rRNA sequences
* `sample_norRNA_R{1,2}.fastq.gz`
  * The fastq file after rRNA cleaning

## Genome Mapping
Raw (or rRNA-cleaned) sequencing data are then aligned on the reference genome.
The current version includes the following RNA reads mappers: TopHat2, STAR, HiSat2.

For details about these mapper, see their help page:
- [`STAR`](https://github.com/alexdobin/STAR) 
- [`tophat2`](http://ccb.jhu.edu/software/tophat/index.shtml) 
- [`hisat2`](http://ccb.jhu.edu/software/hisat2/index.shtml)

> **NB:** Note that the TopHat2 software was included for historical reason, but is now deprecated and replaced by the HiSat2 mapper.

Number of aligned reads (unique and multiple) are then extracted from the aligned (.bam) files.
In practice, we usually observe at least 70% of reads aligned on the reference genome. Sample with lower mapper rate should be carrefully considered and check for contamination, adapter content, sequencing issues, etc.

> **NB:** The final results of the alignment step is a sorted bam file per sample, with its index file (.bai). By default, all other files such as unsorted bam files are not saved.
Use the `--saveAlignedIntermediates` options to save all files.

**Output directory: `results/mapping`**

* `logs/`
  * Log files with mapping statistics
* `sample.bam`
  * The final sorted bam file
* `sample.bam.bai`
  * The index file of the sorted bam file

## Strandness

## Read Distribution

## Complexity Curves

## Gene-based Saturation

## Counts

## Expressed genes

## DupRadar

## PCA

## Correlation


## MultiQC
[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see http://multiqc.info
