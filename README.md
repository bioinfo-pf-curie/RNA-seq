# RNA-seq 

**Institut Curie - Nextflow rna-seq analysis pipeline**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![MultiQC](https://img.shields.io/badge/MultiQC-1.11-blue.svg)](https://multiqc.info/)
[![Install with conda](https://img.shields.io/badge/install%20with-conda-brightgreen.svg)](https://conda.anaconda.org/anaconda)
[![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)](https://singularity.lbl.gov/)
[![Docker Container available](https://img.shields.io/badge/docker-available-003399.svg)](https://www.docker.com/)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7443721.svg)](https://doi.org/10.5281/zenodo.7443721)


### Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. 
It comes with conda / singularity containers making installation easier and results highly reproducible.

The first version of this pipeline was modified from the [nf-core/rnaseq](https://github.com/nf-core/rnaseq) pipeline. 
See the [nf-core](https://nf-co.re/) project for more details.

### Pipline summary

1. Trim adapters from sequencing reads ([`TrimGalore!`](https://github.com/FelixKrueger/TrimGalore)
2. Separate host/graft reads for PDX model ([`xengsort`](https://gitlab.com/genomeinformatics/xengsort))
3. Run quality control of sequencing reads ([`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. Align reads on ribosomal RNAs sequences when available ([`bowtie1`](http://bowtie-bio.sourceforge.net/index.shtml))
5. Align reads on reference genome ([`STAR`](https://github.com/alexdobin/STAR) / [`hisat2`](http://ccb.jhu.edu/software/hisat2/index.shtml))
6. Infer reads orientation ([`rseqc`](http://rseqc.sourceforge.net/))
7. Dedicated quality controls
    - Saturation curves ([`preseq`](http://smithlabresearch.org/software/preseq/) / [`R`](https://www.r-project.org/))
    - Duplicates ([`picard`](https://broadinstitute.github.io/picard/) / [`dupRadar`](https://bioconductor.org/packages/release/bioc/html/dupRadar.html))
    - Reads annotation ([`qualimap`](http://qualimap.conesalab.org/) / [`R`](https://www.r-project.org/))
    - Gene body coverage ([`qualimap`](http://qualimap.conesalab.org/))
8. Identito monitoring based on a list of known polymorphism ([`bcftools`](http://samtools.github.io/bcftools/bcftools.html) / [`R`](https://www.r-project.org/))
9. Generate counts table from aligned data or pseudo-alignment ([`STAR`](https://github.com/alexdobin/STAR) / [`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/) / [`HTSeqCounts`](https://htseq.readthedocs.io/en/release_0.11.1/count.html)/[`salmon`](https://salmon.readthedocs.io/en/latest/salmon.html))
10. Exploratory analysis ([`R`](https://www.r-project.org/))
11. Reference-guided de novo transcripts assembly ([`stringtie`](https://ccb.jhu.edu/software/stringtie/), [`scallop`](https://github.com/Kingsford-Group/scallop))
12. Present all QC results in a final report ([`MultiQC`](http://multiqc.info/))

### Quick help

```bash
N E X T F L O W  ~  version 20.10.0
Launching `main.nf` [awesome_archimedes] - revision: 7f7a25de60
------------------------------------------------------------------------
   ______ _   _   ___ 
   | ___ \ \ | | / _ \
   | |_/ /  \| |/ /_\ \______ ___  ___  __ _ 
   |    /| . ` ||  _  |______/ __|/ _ \/ _` |
   | |\ \| |\  || | | |      \__ \  __/ (_| |
   \_| \_\_| \_/\_| |_/      |___/\___|\__, |
                                          | |
                                          |_|
                     v4.1.0
------------------------------------------------------------------------
------------------------------------------------------------------------
Usage:

The typical command for running the pipeline is as follows:

nextflow run main.nf --reads PATH --samplePlan PATH -profile STRING --genome STRING

MANDATORY ARGUMENTS:
    --genome        STRING                                                                            Name of the reference genome.
    --reads         PATH                                                                              Path to input data (must be surrounded with quotes)
    --samplePlan    PATH                                                                              Path to sample plan (csv format) with raw reads (if `--reads` is not specified)
    --aligner       STRING [star, hisat2]                                                             Tool for reads alignment
    --counts        STRING [star, featureCounts, HTseqCounts, salmon]                                 Tool to use to estimate the raw counts per gene
    --pseudoAligner STRING [salmon]                                                                   Tool for reads pseudo-alignment
		 
INPUTS:
    --singleEnd                                       For single-end input data
    --stranded  STRING [auto, forward, reverse, no]   Library strandness

PREPROCESSING:
    --pdx                Deconvolute host/graft reads for PDX samples
	--trimming           Trim adapters with TrimGalore	

MAPPING:
    --bowtieOpts               STRING    Options for rRNA mapping with bowtie
    --hisat2Opts               STRING    Options for genome mapping with Hisat2
    --saveAlignedIntermediates           Save intermediates alignment files
    --starOpts                 STRING    Options for STAR mapping
    --starTwoPass                        Run STAR in two pass mode

COUNTS:
    --featurecountsOpts STRING   Options for featureCounts quantification
    --htseqOpts         STRING   Options for HTSeq quantification
    --salmonQuantOpts   STRING   Options for Salmon quantification

DE NOVO ASSEMBLY:
    --denovo STRING [stringtie, scallop]  Tool for reference-guided assembly of RNA transcripts
    --scallopOpts  STRING                 Options for Scallop analysis
    --stringtieOps STRING                 Options for Stringtie analysis
		
REFERENCES:
    --bed12                PATH   Path to gene file (BED12)
    --fasta                PATH   Path to genome fasta file
    --fastaFai             PATH   Path to genome index fasta file
    --genomeAnnotationPath PATH   Path to genome annotations folder
    --gtf                  PATH   Path to GTF annotation file
    --hisat2Index          PATH   Path to Hisat2 indexes
    --bowtie2Index         PATH   Path to Bowtie2 indexes
    --polym                PATH   Path to BED file with polymorphisms for identito monitoring
    --rrna                 PATH   Path to Bowtie indexes for rRNA mapping
    --salmonIndex          PATH   Path to Salmon indexes
    --starIndex            PATH   Path to STAR indexes
    --transcriptsFasta     PATH   Path to transcriptome fasta file

SKIP OPTIONS:
    --skipBigWig                       Disable bigwig generation with Deeptools
    --skipDupradar                     Disable duplicates analysis with DupRadar
    --skipFastqc                       Disable Fastqc
    --skipGeneCountsAnalysis           Disable exporatory analysis of genes count
    --skipIdentito                     Disable Identito
    --skipMultiqc                      Disable MultiQC
    --skipQC                           Disable quality controls on raw and aligned reads [fastqc, qualimap, preseq]
    --skipQualimap                     Disable RNA-seq quality controls with Qualimap
    --skipRrna                         Disable rRNA mapping
    --skipSaturation                   Disable saturation analysis with Preseq

OTHER OPTIONS:
    --metadata      PATH     Specify a custom metadata file for MultiQC
    --multiqcConfig PATH     Specify a custom config file for MultiQC
    --name          STRING   Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    --outDir        PATH     The output directory where the results will be saved
	
=======================================================
Available Profiles
  -profile test                        Run the test dataset
  -profile conda                       Build a new conda environment before running the pipeline. Use `--condaCacheDir` to define the conda cache path
  -profile multiconda                  Build a new conda environment per process before running the pipeline. Use `--condaCacheDir` to define the conda cache path
  -profile path                        Use the installation path defined for all tools. Use `--globalPath` to define the insallation path
  -profile multipath                   Use the installation paths defined for each tool. Use `--globalPath` to define the insallation path
  -profile docker                      Use the Docker images for each process
  -profile singularity                 Use the Singularity images for each process. Use `--singularityImagePath` to define the insallation path
  -profile cluster                     Run the workflow on the cluster, instead of locally
						  
```

### Quick run

The pipeline can be run on any infrastructure from a list of input files or from a sample plan as follow

#### Run the pipeline on a test dataset
See the conf/test.conf to set your test dataset.

```
nextflow run main.nf -profile test,conda
```

#### Run the pipeline from a sample plan

```
nextflow run main.nf --samplePlan MY_SAMPLE_PLAN --aligner 'star' --counts 'star' --genome 'hg38' --outDir MY_OUTPUT_DIR -profile conda
```

#### Run the pipeline on a computational cluster

```
echo "nextflow run main.nf --reads '*.R{1,2}.fastq.gz' --aligner 'star' --counts 'star' --genome 'hg19' --outDir MY_OUTPUT_DIR -profile singularity,cluster" | qsub -N rnaseq
```

### Defining the '-profile'

By default (whithout any profile), Nextflow will excute the pipeline locally, expecting that all tools are available from your `PATH` variable.

In addition, we set up a few profiles that should allow you i/ to use containers instead of local installation, ii/ to run the pipeline on a cluster instead of on a local architecture.
The description of each profile is available on the help message (see above).

Here are a few examples of how to set the profile option. See the [full documentation](docs/profiles.md) for details.

```
## Run the pipeline locally, using the paths defined in the configuration for each tool (see conf/path.config)
-profile path --globalPath INSTALLATION_PATH 

## Run the pipeline on the cluster, using the Singularity containers
-profile cluster,singularity --singularityImagePath SINGULARITY_IMAGE_PATH 

## Run the pipeline on the cluster, building a new conda environment
-profile cluster,conda --condaCacheDir CONDA_CACHE 
```

### Sample Plan

A sample plan is a csv file (comma separated) that list all samples with their biological IDs, **with no header**.


SAMPLE_ID,SAMPLE_NAME,PATH_TO_R1_FASTQ,[PATH_TO_R2_FASTQ]


### Full Documentation

1. [Installation](docs/installation.md)
2. [Reference genomes](docs/referenceGenomes.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

#### Credits

This pipeline has been written by the bioinformatics platform of the Institut Curie (P. La Rosa, N. Servant)

#### Citation

If you use this pipeline for your project, please cite it using the following doi: [10.5281/zenodo.7443721](https://doi.org/10.5281/zenodo.7443721)  
Do not hesitate to use the Zenodo doi corresponding to the version you used !

#### Contacts

For any question, bug or suggestion, please use the issues system or contact the bioinformatics core facility.
