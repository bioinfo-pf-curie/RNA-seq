# Usage

## Table of contents

* [Introduction](#general-nextflow-info)
* [Running the pipeline](#running-the-pipeline)
    * [For gene-based quantification](#for-gene-based-quantification)
	* [For isoform-based quantification](#for-isoform-based-quentification)
	* [For reference-guided de-novo assembly](#for-reference-guided-de-novo-assembly)
* [Main arguments](#main-arguments)
    * [`--reads`](#--reads)
	* [`--samplePlan`](#--samplePlan)
    * [`--aligner`](#--aligner)
	* [`--pseudoAligner`](#--pseudoAligner)
    * [`--counts`](#--counts)
    * [`--singleEnd`](#--singleend)
	* [`--stranded`](#--stranded)
* [Reference genomes](#reference-genomes)
    * [`--genome`](#--genome)
	* [`--genomeAnnotationPath`](#--genomeAnnotationPath)
* [Preprocessing](#preprocessing)
    * [`--trimming`](#--trimming)
	* [`--pdx`](#--pdx)
* [Alignment](#aligment)
	* [`--bowtieOpts`](#--bowtieOpts)
    * [`--hisat2Opts`](#--hisat2Opts)
    * [`--starOpts`](#--starOpts)
	* [`--starTwoPass`](#--starTwoPass)
* [Counts](#counts)
	* [`--htseqOpts`](#--htseqOpts)
	* [`--featurecountsOpts`](#--featurecountsOpts)
	* [`--salmonQuantOpts`](#--salmonQuantOpts)
* [Reference-guided Transcriptome Assembly](#reference-guided-transcriptome-assembly)
    * [`--denovo`](#--denovo)
	* [`--scallopOpts`](#--scallopOpts)
	* [`--stringtieOpts`](#--stringtieOpts)
* [Profiles](#profiles)
* [Job resources](#job-resources)
* [Automatic resubmission](#automatic-resubmission)
* [Custom resource requests](#custom-resource-requests)
* [Other command line parameters](#other-command-line-parameters)
    * [`--skip*`](#--skip*)
	* [`--metadata`](#--metadta)
	* [`--outDir`](#--outdir)
    * [`-name`](#-name-single-dash)
    * [`-resume`](#-resume-single-dash)
    * [`-c`](#-c-single-dash)
    * [`--maxMemory`](#--max_memory)
    * [`--maxTime`](#--max_time)
    * [`--maxCpus`](#--max_cpus)
    * [`--multiqcConfig`](#--multiqc_config)

## General Nextflow info

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
nextflow run main.nf --reads '*_R{1,2}.fastq.gz' --aligner 'star' --counts 'star' -profile 'singularity'
```

This will launch the pipeline with the `singularity` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

You can change the output director using the `--outDir/-w` options.

### For gene-based quantification

Gene-based quantification is a very common question in RNA-seq analysis. Broadly speaking, two different strategies can be used:

- Sequencing reads are first aligned on a reference genome, and gene counts are estimated using
tools such as `STAR`, `HTSeqCounts` or `featureCounts`.

```bash
nextflow run main.nf --reads '*_R{1,2}.fastq.gz' --aligner 'star' --counts 'star' -profile 'singularity,cluster'
```

- Or gene abundance can be directly infered from raw sequencing data using pseudo-alignment (or selective-alignment) methods such as `salmon`

```bash
nextflow run main.nf --reads '*_R{1,2}.fastq.gz' --pseudoAligner 'salmon' -profile 'singularity,cluster'
```

In this case, there is no aligned (BAM) file, and the pipeline just take raw fastq files and directly extract counts table.

Currently, [several studies](https://f1000research.com/articles/4-1521/v1) have shown that transcript quantification approaches (such as `salmon`) are more sensitive than read 
counting methods (such as `HTSeqCounts` or `featureCounts`), although most of these demonstrations are made from simulated data. 
So far, the current best practice would thus be to favor the usage of `salmon` over other tools.

Then regarding the differences / benefits of running `salmon` in alignment mode (from a BAM) versus in selective-alignment mode (from raw reads), this point is discussed in the [mapping and alignment methodology paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8). The most obvious difference is that modes which are not based on a BAM file are much faster. On the other hand, it is usually useful to have a BAM file to perform additional analysis and visualization of the data. Regarding accuracy, it doesn't seem to have huge difference and as long as the `salmon` indexes are properly built (considering the genome as "decoy" sequence), both methods usually lead to accurate quantification estimates.


### For isoform-based quantification

One of the main interest of `salmon` is its ability to estimate the abundance at both genes and (known) transcripts levels.  
If you are interested in isoform analysis, it could also be useful to run `STAR` in a two-pass mode. Here, the idea is to run a first alignment with usual parameters, 
then collect the junctions detected and use them as "annotated" junctions for the second mapping pass.  
In addition, all default tools' parameters can be updated on the command line. As an example here, we would like to update the `--numBootstraps` parameters which is
required to run tools such as `sleuth` for differential analysis as the isoform levels.

The typical command line to estimate both known isoform and gene abundances with `salmon` would be ;

```bash
nextflow run main.nf --reads '*_R{1,2}.fastq.gz' --aligner 'star' --starTwoPass --counts 'salmon' --salmonQuantOpts '--numBootstraps 200' -profile 'singularity,cluster'
```

### For reference-guided de-novo assembly

Since v4.0.0, the pipeline now includes tools for reference-guided de-novo assembly.  
The goal of such analysis is to detect new isoform/genes from short reads data. The typical output is a new `gtf` file with known and new transcripts.

In the current version, `scallop` and `stringtie` are available and can be specified using the `--denovo` option. 
Note that both methods require a BAM file as input. Multiple tools can be specificed (comma separated).  

The results are then assessed using the `gffCompare` utility which will compare a know `gtf` with the one(s) generated by the pipeline.
In most of the case, a high fraction is detection transcripts should correspond to known ones. `GffCompare` thus proposes sensitivity/specificy metrics for [accuracy estimation](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml).

Here again, it is recommanded to run `STAR` in two-pass mode to [improve novel splice junction detection](https://academic.oup.com/bioinformatics/article/32/1/43/1744001).

```bash
nextflow run main.nf --reads '*_R{1,2}.fastq.gz' --aligner 'star' --starTwoPass --denovo 'stringtie,scallop' -profile 'singularity,cluster'
```

## Main arguments

### `--reads`

Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### `--samplePlan`

Use this to specify a sample plan file instead of a regular expression to find fastq files. For example :

```bash
--samplePlan 'path/to/data/sample_plan.csv
```

The sample plan is a csv file with the following information :

Sample ID | Sample Name | Path to R1 fastq file | Path to R2 fastq file

### `--aligner`

The current version of the pipeline supports two different aligners;
- [`STAR`](https://github.com/alexdobin/STAR)
- [`hisat2`](http://ccb.jhu.edu/software/hisat2/index.shtml)

Since version `4.0.0`, no default value is defined. You can specify the tool to use as follows:

```bash
--aligner 'STAR'
```

### `--pseudoAligner`

Recent advances in the field of RNA-seq data analysis promote the usage of pseudo-aligner which are able to 
estimate the gene/transcripts abundance directly from the raw fastq files.

The following tools are currently available;
- [`salmon`](https://salmon.readthedocs.io/en/latest/salmon.html)

```bash
--pseudoAligner 'salmon'
```

### `--counts`

The raw count table for all samples can be generated from alignment file (bam) using one of the following tool:
- [`STAR`](https://github.com/alexdobin/STAR). Require `--aligner 'STAR'`. Default value.
- [`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/)
- [`HTSeqCounts`](https://htseq.readthedocs.io/en/release_0.11.1/count.html)
- [`salmon`](https://salmon.readthedocs.io/en/latest/salmon.html)

Since version `4.0.0`, no default value is defined. You can specify one of these tools using:
```bash
--counts 'salmon`
```

### `--singleEnd`

By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--singleEnd --reads '*.fastq.gz'
```

It is not possible to run a mixture of single-end and paired-end files in one run.

### `--stranded`

Several parts of the RNA-seq data analysis rely on the strandness information.
If you already have the information, you should specifiy the strandness using either `forward` (ie. stranded), `no`, (ie. unstranded), `reverse` (ie. reverse stranded), as follows:

```bash
--stranded 'reverse'
```

If you do not have the information, you can the automatic detection mode (default mode) as follows:

```bash
--stranded 'auto'
```

In the case, the pipeline will the run the [`rseqc`](http://rseqc.sourceforge.net/) tool to automatically detect the strandness parameter.

## Reference genomes

The pipeline config files come bundled with paths to the genomes reference files. 

### `--genome`

There are different species supported in the genomes references file. To run the pipeline, you must specify which to use with the `--genome` flag.

You can find the keys to specify the genomes in the [genomes config file](../conf/genomes.config). Common genomes that are supported are:

* Human
  * `--genome hg38`
* Mouse
  * `--genome mm10`

> There are numerous others - check the config file for more.

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the genomes resource. 
See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'hg19' {
      fasta            = '<path to genome fasta file for identito monitoring>'
      fastaFai         = '<path to genome index file for identito monitoring>'
      bowtie2          = '<path to the bowtie2 index files>' 
      star             = '<path to the STAR index files>'
      hisat2           = '<path to the HiSat2 index files>'
      salmon           = '<path to the Salmon index files>'
      rrna             = '<path to bowtie1 mapping on rRNA reference>'
      bed12            = '<path to Bed12 annotation file>'
      gtf              = '<path to GTF annotation file>'
      transcriptsFasta = '<path to fasta transcriptome file for pseudo-alignment>'
      gencode          = '<boolean - is the annotation file based on Gencode ?'
      polym            = '<path to bed file for identito monitoring>'
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

Note that these paths can be updated on command line using the following parameters:
- `--bowtie2Index` - Path to bowtie2 index
- `--starIndex` - Path to STAR index
- `--hisat2Index` - Path to HiSAT2 index
- `--salmonIndex` - Path to Salmon index
- `--gtf` - Path to GTF file
- `--gencode` - Specify that the GTF file is from Gencode
- `--bed12` - Path to gene bed12 file
- `--transcriptsFasta` - Path to transcriptome in fasta format
- `--saveAlignedIntermediates` - Save the BAM files from the Aligment step  - not done by default

## Preprocessing

### `--trimming`

The raw data trimming can be performed with [TrimGalore](https://github.com/FelixKrueger/TrimGalore) if the `--trimming` option is specified.
By default, `TrimGalore` should be able to automatically detect the Illumna 3' adapter. If more advance trimming parameter are required, plese use the option `--trimmingOpts`.
The trimmed fastq files are then used for downstream analysis but are not exported by default. Use the `--saveIntermediates` parameters to export them.

### `--pdx`

In the context of Mouse xenograft samples, it is strongly recommanded to distinguish Mouse from Human reads in order to avoid data misalignment.
To do so, we implemented the [`xengsort`](https://gitlab.com/genomeinformatics/xengsort) tool (`--pdx`) which generates in output distinct fastq files for both genomes.
These new fastq files will then be used for downstream alignment and analysis.

According to the `--genome` option specified, the analysis will be either performed from the Human or from the Mouse fastq files.

## Alignment

### `--bowtieOpts`

Change default bowtie1 mapping options for rRNA cleaning. See the `nextflow.config` file for details.

### `--hisat2Opts`

Change default Hisat2 mapping option. See the `nextflow.config` file for details.

### `--starOpts`
 
Change default STAR mapping options for mapping.
By default STAR is run with the following options

```
--outSAMmultNmax 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --outSAMprimaryFlag OneBestScore \
--outMultimapperOrder Random --outSAMattributes All
```

In other words, it means that only one alignment will be reported in the output, randomly chosen from the top scoring alignments (in case of multiple alignments).
The allowed number of mismatches is indexed on the read length (0.04 * read length). 
And all common SAM attributes will be added.

Note that the default STAR options can vary from an organism to another. 
For instance, for Human data, the pipeline adds the ENCODE recommanded options to the default ones 

```
--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000
```

See the `nextflow.config` file for details.

### `--starTwoPass`

Run the STAR aligner in two pass mode with the `--twopassMode Basic` option


## Counts

### `--htseqOpts`

Change default HTSeq options. See the `nextflow.config` file for details.

### `--featurecountsOpts`

Change default featureCounts options. See the `nextflow.config` file for details.

### `--salmonQuantOpts`

Change default options for Salmon quantification from either BAM file or FASTQ files. See the `nextflow.config` file for details.

## Reference-guided Transcriptome Assembly

### `--denovo`

Specify which tools to use for reference-guided transcriptome assembly. Several tools can be specified (comma separated)

```
--denovo 'scallop,stringtie'
```

### `--scallopOpts`

Change default scallop options. See the `nextflow.config` file for details.

### `--stringtieOps`

Change default stringtie options. See the `nextflow.config` file for details.			


## Profiles

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. 
Note that multiple profiles can be loaded, for example: `-profile docker` - the order of arguments is important!
Look at the [Profiles documentation](profiles.md) for details.

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time (see the `conf/base.conf` file). 
For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

## Other command line parameters

### `--skip*`

The pipeline is made with a few *skip* options that allow to skip optional steps in the workflow.
The following options can be used:
- `--skipQC` - Skip all QC steps apart from MultiQC
- `--skipRrna` - Skip rRNA mapping
- `--skipFastqc` - Skip FastQC
- `--skipQualimap` - Skip genebody coverage step
- `--skipSaturation` - Skip Saturation qc
- `--skipDupradar` - Skip dupRadar (and Picard MarkDups)
- `--skipExpan` - Skip exploratory analysis
- `--skipBigwig` - Do not generate bigwig files
- `--skipIdentito` - Skip identito monitoring
- `--skipMultiqc` - Skip MultiQC
- `--skipSoftVersions` - Skip software versions reporting
			
			
### `--metadata`

Specify a two-columns (tab-delimited) metadata file to diplay in the final Multiqc report.

### `--outDir`

The output directory where the results will be saved.

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.
You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--maxMemory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--maxTime`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--maxCpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--multiqcConfig`

Specify a path to a custom MultiQC configuration file.
