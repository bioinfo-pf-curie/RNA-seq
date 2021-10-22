// RNA-seq Analysis Pipeline : functions 
//  helpMessage 

def helpMessage() {
  if ("${workflow.manifest.version}" =~ /dev/ ){
    devMess = file("$baseDir/assets/devMessage.txt")
    log.info devMess.text
  }

  log.info """
  rnaseq v${workflow.manifest.version}
  ======================================================================

  Usage:
  nextflow run rnaseq --reads '*_R{1,2}.fastq.gz' --genome hg19 -profile conda
  nextflow run rnaseq --samplePlan sample_plan --genome hg19 -profile conda

  Mandatory arguments:
    --reads [file]                       Path to input data (must be surrounded with quotes)
    --samplePlan [file]                  Path to sample plan input file (cannot be used with --reads)
    --genome [str]                       Name of genome reference
    -profile [str]                       Configuration profile to use. test / conda / toolsPath / singularity / cluster (see below)

  Inputs:
    --singleEnd [bool]                   Specifies that the input is single end reads

  Strandness:
    --stranded [bool]                    Library strandness ['auto', 'forward', 'reverse', 'no']. Default: 'auto'

  Mapping:
    --aligner [str]                      Tool for read alignments ['star', 'hisat2']. Default: 'star'

  Counts:
    --counts [str]                       Tool to use to estimate the raw counts per gene ['star', 'featureCounts', 'HTseqCounts']. Default: 'star'

  References: If not specified in the configuration file or you wish to overwrite any of the references.
    --genomeAnnotationPath [file]        Path  to genome annotation folder
    --fasta [file]                       Path the genome fasta file
    --starIndex [dir]                    Path to STAR index
    --hisat2Index [file]                 Path to HiSAT2 index
    --gtf [file]                         Path to GTF file
    --bed12 [file]                       Path to gene bed12 file
    --polym [file]                       Path to BED file with polym to check
    --saveAlignedIntermediates [bool]    Save the intermediate files from the Aligment step. Default: false

  Other options:
    --metadata [file]                    Add metadata file for multiQC report
    --outDir [dir]                       The output directory where the results will be saved
    -w/--work-dir [dir]                  The temporary directory where intermediate data will be saved
    -name [str]                          Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

  Skip options:
    --skipQC [bool]                      Skip all QC steps apart from MultiQC
    --skipRrna [bool]                    Skip rRNA mapping
    --skipFastqc [bool]                  Skip FastQC
    --skipSaturation [bool]              Skip Saturation qc
    --skipDupradar [bool]                Skip dupRadar (and Picard MarkDups)
    --skipQualimap [bool]                Skip Qualimap analysis
    --skipExpan [bool]                   Skip exploratory analysis
    --skipBigwig [bool]                  Skip bigwig files 
    --skipIdentito [bool]                Skip identito checks
    --skipMultiQC [bool]                 Skip MultiQC
    --skipSoftVersions [bool]            Skip getSoftwareVersion

  =======================================================
  Available Profiles
    -profile test                        Run the test dataset
    -profile conda                       Build a new conda environment before running the pipeline. Use `--condaCacheDir` to define the conda cache path
    -profile multiconda                  Build a new conda environment per process before running the pipeline. Use `--condaCacheDir` to define the conda cache path
    -profile path                        Use the installation path defined for all tools. Use `--globalPath` to define the insallation path
    -profile multipath                   Use the installation paths defined for each tool. Use `--globalPath` to define the insallation path
    -profile docker                      Use the Docker images for each process
    -profile singularity                 Use the Singularity images for each process. Use `--singularityPath` to define the insallation path
    -profile cluster                     Run the workflow on the cluster, instead of locally

  """.stripIndent()
}

