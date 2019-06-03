# RNA-seq 

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

### Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. 
It comes with docker / singularity containers making installation trivial and results highly reproducible.

The current workflow is based on the nf-core best practice. See the nf-core project from details on [guidelines](https://nf-co.re/).

### Pipline summary

1. Run quality control of raw sequencing reads (fastqc)
2. Align reads on ribosomal RNAs sequences when available (bowtie1)
3. Align reads on reference genome (STAR/Tophat2/hisat)
4. Infer reads orientation (rseqc)
5. Dedicated quality controls
  - Saturation curves (preseq)
  - Duplicates (dupRadar)
  - Reads annotation (rseqc)
6. Generate counts table (STAR/featureCounts/HTSeqCounts)
7. Exploratory analysis (R)

### Quick help

```bash
nextflow run main.nf --help
N E X T F L O W  ~  version 19.04.0
Launching `main.nf` [stupefied_darwin] - revision: aa905ab621
rnaseq v2.0.0dev
=======================================================

Usage:
nextflow run rnaseq --reads '*_R{1,2}.fastq.gz' --genome 'hg19' 

Mandatory arguments:
  --reads                       Path to input data (must be surrounded with quotes)
  -profile                      Configuration profile to use. test / curie / conda / docker / singularity

Options:
  --genome                      Name of genomes reference
  --singleEnd                   Specifies that the input is single end reads

Strandedness:
  --stranded                    Library strandness ['auto', 'yes', 'reverse', 'no']. Default: 'auto'

Mapping:
  --aligner                     Tool for read alignments ['star', 'hisat2', 'tophat2']. Default: 'star'

Counts:
  --counts                      Tool to use to estimate the raw counts per gene ['star', 'featureCounts', 'HTseqCounts']. Default: 'star'

References:                     If not specified in the configuration file or you wish to overwrite any of the references.
  --star_index                  Path to STAR index
  --hisat2_index                Path to HiSAT2 index
  --tophat2_index		    Path to TopHat2 index
  --gtf                         Path to GTF file
  --bed12                       Path to gene bed12 file
  --saveAlignedIntermediates    Save the BAM files from the Aligment step  - not done by default

Other options:
  --metadata                    Add metadata file for multiQC report
  --outdir                      The output directory where the results will be saved
  -w/--work-dir                 The temporary directory where intermediate data will be saved
  --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
  -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

QC options:
  --skip_qc                     Skip all QC steps apart from MultiQC
  --skip_rrna                   Skip rRNA mapping
  --skip_fastqc                 Skip FastQC
  --skip_preseq                 Skip Preseq
  --skip_dupradar               Skip dupRadar (and Picard MarkDups)
  --skip_read_dist              Skip read distribution step
  --skip_expan                  Skip exploratory analysis
  --skip_multiqc                Skip MultiQC

```

### Quick run

The pipeline can be run on any infrastructure as follow

#### Run the pipeline locally

```
nextflow run rnaseq --reads '*_R{1,2}.fastq.gz' --genome 'hg19' --outdir MY_OUTPUT_DIR

```

#### Run the pipeline on the Institut Curie cluster

```
echo "nextflow run main.nf --reads '*.R{1,2}.fastq.gz' --genome 'hg19' --outdir MY_OUTPUT_DIR --queue 'batch'" | qsub -q mpi -N rnaseq-2.0

```

### Full Documentation

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Reference genomes](docs/configuration/reference_genomes.md)  
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

### Credits
This pipeline has been written by the bioinformatics platform of the Institut Curie (P. La Rosa, N. Servant)
