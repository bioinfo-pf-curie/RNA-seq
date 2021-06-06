#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2019-2020
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data.
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

This script is based on the nf-core guidelines. See https://nf-co.re/ for more information
*/


/*
========================================================================================
                         RNA-seq
========================================================================================
 RNA-seq Analysis Pipeline.
 #### Homepage / Documentation
 https://gitlab.curie.fr/data-analysis/rnaseq
----------------------------------------------------------------------------------------
*/


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

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
  helpMessage()
  exit 0
}

// Configurable reference genomes
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
  exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Reference index path configuration
// Define these here - after the profiles are loaded with the genomes paths
params.starIndex = params.genome ? params.genomes[ params.genome ].star ?: false : false
params.bowtie2Index = params.genome ? params.genomes[ params.genome ].bowtie2 ?: false : false
params.hisat2Index = params.genome ? params.genomes[ params.genome ].hisat2 ?: false : false
params.rrna = params.genome ? params.genomes[ params.genome ].rrna ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.bed12 = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.polym = params.genome ? params.genomes[ params.genome ].polym ?: false : false

// Tools option configuration
// Add here the list of options that can change from a reference genome to another
if (params.genome){
  params.starOptions = params.genomes[ params.genome ].starOpts ?: params.starOpts
}
// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
customRunName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  customRunName = workflow.runName
}

// Stage config files
chMultiqcConfig = Channel.fromPath(params.multiqcConfig)
chOutputDocs = Channel.fromPath("$baseDir/docs/output.md")
chOutputDocsImages = file("$baseDir/docs/images/", checkIfExists: true)
chPcaHeader = Channel.fromPath("$baseDir/assets/pcaHeader.txt")
chHeatmapHeader = Channel.fromPath("$baseDir/assets/heatmapHeader.txt")

/*
 * CHANNELS
 */

// Validate inputs
if (params.aligner != 'star' && params.aligner != 'hisat2'){
  exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star', 'hisat2'"
}
if (params.counts != 'star' && params.counts != 'featureCounts' && params.counts != 'HTseqCounts'){
  exit 1, "Invalid counts option: ${params.counts}. Valid options: 'star', 'featureCounts', 'HTseqCounts'"
}
if (params.counts == 'star' && params.aligner != 'star'){
  exit 1, "Cannot run STAR counts without STAR aligner. Please check the '--aligner' and '--counts' parameters."
}
if (params.stranded != 'auto' && params.stranded != 'reverse' && params.stranded != 'forward' && params.stranded != 'no'){
  exit 1, "Invalid stranded option: ${params.stranded}. Valid options: 'auto', 'reverse', 'forward', 'no'"
}
if (params.stranded == 'auto' && !params.bed12){  
  exit 1, "Strandness detection is not possible without gene bed file. Please specify a stranded option: 'reverse', 'forward', 'no'"
}

if ((params.reads && params.samplePlan) || (params.readPaths && params.samplePlan)){
  exit 1, "Input reads must be defined using either '--reads' or '--samplePlan' parameter. Please choose one way"
}

if( params.starIndex && params.aligner == 'star' ){
  Channel
    .fromPath(params.starIndex)
    .ifEmpty { exit 1, "STAR index not found: ${params.starIndex}" }
    .set {chStarIndex}
}
else if ( params.hisat2Index && params.aligner == 'hisat2' ){
  Channel
    .fromPath("${params.hisat2Index}*")
    .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2Index}" }
    .set{chHisat2Index}
}
else {
    exit 1, "No reference genome specified!"
}

if( params.gtf ){
  Channel
    .fromPath(params.gtf)
    .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
    .into { chGtfStar; chGtfDupradar; chGtfFeatureCounts; chGtfGenetype; chGtfHTseqCounts; chGtfTable; chGtfMakeHisatSplicesites; chGtfQualimap }
}else {
  log.warn "No GTF annotation specified - dupRadar, table counts - will be skipped !" 
  Channel
    .empty()
    .into { chGtfStar; chGtfDupradar; chGtfFeatureCounts; chGtfGenetype; chGtfHTseqCounts; chGtfTable; chGtfMakeHisatSplicesites; chGtfQualimap } 
}

if( params.bed12 ){
  Channel
    .fromPath(params.bed12)
    .ifEmpty { exit 1, "BED12 annotation file not found: ${params.bed12}" }
    .set { chBedRseqc }
}else{
  log.warn "No BED gene annotation specified - strandness detection, gene body coverage, read distribution - will be skipped !"
  Channel
    .empty()
    .set { chBedRseq }
}

if( params.rrna ){
  Channel
    .fromPath(params.rrna)
    .ifEmpty { exit 1, "rRNA annotation file not found: ${params.rrna}" }
    .set { chRrnaAnnot }
}else{
  log.warn "No rRNA fasta file available - rRNA mapping - will be skipped !"
  chRrnaAnnot = Channel.empty()
}

if ( params.metadata ){
  Channel
    .fromPath( params.metadata )
    .ifEmpty { exit 1, "Metadata file not found: ${params.metadata}" }
    .set { chMetadata }
}

if ( params.fasta ){
  Channel
    .fromPath( params.fasta )
    .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
    .set { chFasta }
}else {
  chFasta = Channel.empty()
}

if ( params.polym ){
  Channel
    .fromPath( params.polym )
    .ifEmpty { exit 1, "Polym BED file not found: ${params.polym}" }
    .set { chPolymBed }
}else{
  chPolymBed = Channel.empty()
}

/*
 * Create a channel for input read files
 */

if(params.samplePlan){
  if(params.singleEnd){
    Channel
      .from(file("${params.samplePlan}"))
      .splitCsv(header: false)
      .map{ row -> [ row[0], [file(row[2])]] }
      .into { chRawReadsFastqc; chRawReadsStar; chRawReadsHisat2; chRawReadsRnaMapping; chRawReadsPrepRseqc; chRawReadsStrandness; chSaveStrandness }
  }else{
     Channel
       .from(file("${params.samplePlan}"))
       .splitCsv(header: false)
       .map{ row -> [ row[0], [file(row[2]), file(row[3])]] }
       .into { chRawReadsFastqc; chRawReadsStar; chRawReadsHisat2; chRawReadsRnaMapping; chRawReadsPrepRseqc; chRawReadsStrandness; chSaveStrandness }
   }
   params.reads=false
}
else if(params.readPaths){
  if(params.singleEnd){
    Channel
      .from(params.readPaths)
      .map { row -> [ row[0], [file(row[1][0])]] }
      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
      .into { chRawReadsFastqc; chRawReadsStar; chRawReadsHisat2; chRawReadsRnaMapping; chRawReadsPrepRseqc; chRawReadsStrandness; chSaveStrandness }
  } else {
    Channel
      .from(params.readPaths)
      .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
      .into { chRawReadsFastqc; chRawReadsStar; chRawReadsHisat2; chRawReadsRnaMapping; chRawReadsPrepRseqc; chRawReadsStrandness; chSaveStrandness }
  }
} else {
  Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { chRawReadsFastqc; chRawReadsStar; chRawReadsHisat2; chRawReadsRnaMapping; chRawReadsPrepRseqc; chRawReadsStrandness; chSaveStrandness }
}

/*
 * Make sample plan if not available
 */

if (params.samplePlan){
  chSplan = Channel.fromPath(params.samplePlan)
}else if(params.readPaths){
  if (params.singleEnd){
    Channel
      .from(params.readPaths)
      .collectFile() {
        item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
       }
      .set{ chSplan }
  }else{
     Channel
       .from(params.readPaths)
       .collectFile() {
         item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
        }
       .set{ chSplan }
  }
}else{
  if (params.singleEnd){
    Channel
      .fromFilePairs( params.reads, size: 1 )
      .collectFile() {
         item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
      }     
      .set { chSplan }
  }else{
    Channel
      .fromFilePairs( params.reads, size: 2 )
      .collectFile() {
         item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
      }     
      .set { chSplan }
  }
}


// Header log info
if ("${workflow.manifest.version}" =~ /dev/ ){
   devMess = file("$baseDir/assets/devMessage.txt")
   log.info devMess.text
}

log.info """=======================================================

 rnaseq : RNA-Seq workflow v${workflow.manifest.version}
======================================================="""
def summary = [:]
summary['Run Name']     = customRunName ?: workflow.runName
summary['Command Line'] = workflow.commandLine
summary['Metadata']	= params.metadata
if (params.samplePlan) {
   summary['SamplePlan']   = params.samplePlan
}else{
   summary['Reads']        = params.reads
}
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Genome']       = params.genome
summary['Strandedness'] = params.stranded
if(params.aligner == 'star'){
  summary['Aligner'] = "star"
  if(params.starIndex) summary['STAR Index'] = params.starIndex
} else if(params.aligner == 'hisat2') {
  summary['Aligner'] = "HISAT2"
  if(params.hisat2Index) summary['HISAT2 Index'] = params.hisat2Index
}
summary['Counts'] = params.counts
if(params.gtf)  summary['GTF Annotation']  = params.gtf
if(params.bed12) summary['BED Annotation']  = params.bed12
if(params.polym) summary['BED Polym']  = params.polym
summary['Save Intermeds'] = params.saveAlignedIntermediates ? 'Yes' : 'No'
summary['Max Memory']     = params.maxMemory
summary['Max CPUs']       = params.maxCpus
summary['Max Time']       = params.maxTime
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outDir
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


/*
 * FastQC
 */
process fastqc {
  tag "${prefix}"
  label 'fastqc'
  label 'medCpu'
  label 'lowMem'
  publishDir "${params.outDir}/fastqc", mode: 'copy',
    saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

  when:
  !params.skipQC && !params.skipFastqc

  input:
  set val(prefix), file(reads) from chRawReadsFastqc

  output:
  file "*_fastqc.{zip,html}" into chFastqcResults
  file("v_fastqc.txt") into chFastqcVersion

  script:
  pbase = reads[0].toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
  """
  fastqc --version &> v_fastqc.txt
  fastqc -q $reads
  mv ${pbase}_fastqc.html ${prefix}_fastqc.html
  mv ${pbase}_fastqc.zip ${prefix}_fastqc.zip
  """
}


/*
 * rRNA mapping 
 */
process rRNAMapping {
  tag "${prefix}"
  label 'bowtie'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outDir}/rRNAmapping", mode: 'copy',
    saveAs: {filename ->
      if (filename.indexOf("fastq.gz") > 0 &&  params.saveAlignedIntermediates) filename
      else if (filename.indexOf(".log") > 0) "logs/$filename"
      else null
    }

  when:
  !params.skipRrna && params.rrna

  input:
  set val(prefix), file(reads) from chRawReadsRnaMapping
  file annot from chRrnaAnnot.collect()

  output:
  set val(prefix), file("*fastq.gz") into chRrnaMappingRes
  set val(prefix), file("*.sam") into chRrnaSam
  file "*.log" into chRrnaLogs
  file("v_bowtie.txt") into chBowtieVersion

  script:
  inputOpts = params.singleEnd ? "${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
  """
  bowtie --version &> v_bowtie.txt
  bowtie ${params.bowtieOpts} \\
         -p ${task.cpus} \\
         --un ${prefix}_norRNA.fastq \\
         --sam ${params.rrna} \\
         ${inputOpts} \\
         ${prefix}.sam  2> ${prefix}.log && \
  gzip -f ${prefix}_norRNA*.fastq 
  """
}


/*
 * Strandness
 */

process saveStrandness {
  label 'unix'
  label 'lowCpu'
  label 'lowMem'
  publishDir "${params.outDir}/strandness" , mode: 'copy',
    saveAs: {filename ->
      if (filename.indexOf(".txt") > 0) "$filename"
      else null
    }
  
  when:
  params.stranded == 'reverse' || params.stranded == 'no' || params.stranded == 'yes' || (params.stranded == 'auto' && !params.bed12)

  input:
  set val(prefix), file(reads) from chSaveStrandness

  output:
  file "*.txt" into chSavedStrandness

  script:
  """
  echo ${params.stranded} > ${prefix}_strandness.txt
  """
}


// User defined
if (params.stranded == 'reverse' || params.stranded == 'forward' || params.stranded == 'no'){
  chRawReadsStrandness
    .map { file ->
           def key = params.stranded
           return tuple(key)
    }
    .into { chStrandedResultsBigwig ; chStrandedResultsFeatureCounts; chStrandedResultsGenetype; chStrandedResultsHTseqCounts;
            chStrandedResultsDupradar; chStrandedResultsHisat; chStrandedResultsTable; chStrandedResultsQualimap }
   chBowtie2Version = Channel.empty()
   chRseqcVersionInferExperiment = Channel.empty()  
}else if (params.stranded == 'auto' && params.bed12){
 
  // auto
  process prepRseqc {
    tag "${prefix}"
    label 'bowtie2'
    label 'medCpu'
    label 'medMem'

    input:
    set val(prefix), file(reads) from chRawReadsPrepRseqc

    output:
    set val("${prefix}"), file("${prefix}_subsample.bam") into chBamRseqc
    file("v_bowtie2.txt") into chBowtie2Version

    script:
    inputOpts = params.singleEnd ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    bowtie2 --version &> v_bowtie2.txt
    bowtie2 --fast --end-to-end --reorder \\
            -p ${task.cpus} \\
            -u ${params.nCheck} \\
            -x ${params.bowtie2Index} \\
            ${inputOpts} > ${prefix}_subsample.bam 
     """
   }

  process rseqc {
    tag "${prefix - '_subsample'}"
    label 'rseqc'
    label 'medCpu'
    label 'lowMem'
    publishDir "${params.outDir}/strandness" , mode: 'copy',
      saveAs: {filename ->
        if (filename.indexOf(".txt") > 0) "$filename"
        else null
      }

    input:
    set val(prefix), file(bamRseqc) from chBamRseqc
    file bed12 from chBedRseqc.collect()

    output:
    file "${prefix}*.{txt,pdf,r,xls}" into chRseqcResults
    stdout into ( chStrandedResultsBigwig, chStrandedResultsFeatureCounts, chStrandedResultsGenetype, chStrandedResultsHTseqCounts,
                  chStrandedResultsDupradar, chStrandedResultsHisat, chStrandedResultsTable, chStrandedResultsQualimap)
    file("v_rseqc.txt") into chRseqcVersionInferExperiment

    script:
    """
    infer_experiment.py --version &> v_rseqc.txt    
    infer_experiment.py -i $bamRseqc -r $bed12 > ${prefix}.txt
    parse_rseq_output.sh ${prefix}.txt > ${prefix}_strandness.txt
    cat ${prefix}_strandness.txt
    """  
  }
}

chStrandnessResults = Channel.empty()
if (params.stranded == 'auto' && params.bed12){
  chStrandnessResults = chRseqcResults
}else{
  chStrandnessResults = chSavedStrandness
}


/*
 * Reads mapping
 */

// From nf-core
// Function that checks the alignment rate of the STAR output
// and returns true if the alignment passed and otherwise false
skippedPoorAlignment = []
def checkStarLog(logs) {
  def percentAligned = 0;
  logs.eachLine { line ->
    if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
      percentAligned = matcher[0][1]
    }else if ((matcher = line =~ /Uniquely mapped reads number\s*\|\s*([\d\.]+)/)) {
      numAligned = matcher[0][1]
    }
  }
  logname = logs.getBaseName() - 'Log.final'
  if(percentAligned.toFloat() <= '2'.toFloat() || numAligned.toInteger() <= 1000.toInteger() ){
      log.info "#################### VERY POOR ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($logname)    >> ${percentAligned}% <<"
      skippedPoorAlignment << logname
      return false
  } else {
      log.info "          Passed alignment > star ($logname)   >> ${percentAligned}% <<"
      return true
  }
}

// Update input channel
chStarRawReads = Channel.empty()
if( params.rrna && !params.skipRrna){
  chStarRawReads = chRrnaMappingRes
} else {  
  chStarRawReads = chRawReadsStar
}


// STAR
if(params.aligner == 'star'){
  chHisat2Version = Channel.empty()

  process star {
    tag "$prefix"
    label 'star'
    label 'highCpu'
    label 'highMem'
    publishDir "${params.outDir}/mapping", mode: 'copy',
      saveAs: {filename ->
        if (filename.indexOf(".bam") == -1) "logs/$filename"
        else if (params.saveAlignedIntermediates) filename
        else null
      }
    publishDir "${params.outDir}/counts", mode: 'copy',
      saveAs: {filename ->
        if (filename.indexOf("ReadsPerGene.out.tab") > 0) "$filename"
        else null
      }

    input:
    set val(prefix), file(reads) from chStarRawReads
    file index from chStarIndex.collect()
    file gtf from chGtfStar.collect().ifEmpty([])

    output:
    set val(prefix), file ("*Log.final.out"), file ('*.bam') into chStarSam
    file "*.out" into chAlignmentLogs
    file "*.out.tab" into chStarLogCounts
    file "*Log.out" into chStarLog
    file "*ReadsPerGene.out.tab" optional true into chStarCountsToMerge, chStarCountsToR
    file("v_star.txt") into chStarVersion

    script:
    def starCountOpt = params.counts == 'star' && params.gtf ? params.starOptsCounts : ''
    def starGtfOpt = params.gtf ? "--sjdbGTFfile $gtf" : ''
    """
    STAR --version &> v_star.txt
    STAR --genomeDir $index \\
         ${starGtfOpt} \\
         --readFilesIn $reads  \\
         --runThreadN ${task.cpus} \\
         --runMode alignReads \\
         --outSAMtype BAM Unsorted  \\
         --readFilesCommand zcat \\
         --runDirPerm All_RWX \\
         --outTmpDir "${params.starTmpDir}"\\
         --outFileNamePrefix $prefix  \\
         --outSAMattrRGline ID:$prefix SM:$prefix LB:Illumina PL:Illumina  \\
         ${params.starOptions} \\
	 --limitOutSJcollapsed 5000000 \\
	 ${starCountOpt}
    """
  }

  process starSort {
    tag "$prefix"
    label 'samtools'
    label 'medCpu'
    label 'medMem'
    publishDir "${params.outDir}/mapping", mode: 'copy'
 
    input:
    set val(prefix), file(LogFinalOut), file (starBam) from chStarSam

    output:
    set file("${prefix}Log.final.out"), file ("*.{bam,bam.bai}") into chStarAligned
    file "${prefix}_sorted.bam.bai"
    file("v_samtools.txt") into chSamtoolsVersionSort

    script:
    """
    samtools --version &> v_samtools.txt
    samtools sort  \\
        -@  ${task.cpus}  \\
        -m ${params.sortMaxMemory} \\
        -o ${prefix}_sorted.bam  \\
        ${starBam}
    samtools index ${prefix}_sorted.bam
    """
    }

    // Filter removes all 'aligned' channels that fail the check
    chStarAligned
      .filter { logs, bams -> checkStarLog(logs) }
      .map { logs, bams -> bams }
      .dump (tag:'starbams')
      .into { chBamBigwig; chBamCount; chBamPreseq; chBamMarkduplicates; chBamFeaturecounts; chBamQualimap;
              chBamGenetype; chBamHTseqCounts }
}


// HiSat2
chHisat2RawReads = Channel.empty()
if( params.rrna && !params.skipRrna ){
    chHisat2RawReads = chRrnaMappingRes
}else {
    chHisat2RawReads = chRawReadsHisat2 
}

if(params.aligner == 'hisat2'){
  chStarLog = Channel.empty()
  chStarVersion = Channel.empty()  

  process makeHisatSplicesites {
     label 'hisat2'
     label 'lowCpu'
     label 'lowMem'
     publishDir "${params.outDir}/mapping", mode: 'copy',
       saveAs: { filename ->
         if (params.saveAlignedIntermediates) filename
         else null
       }

     input:
     file gtf from chGtfMakeHisatSplicesites

     output:
     file "${gtf.baseName}.hisat2SpliceSites.txt" into chAlignmentSplicesites

     script:
     """
     hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.hisat2SpliceSites.txt
     """
  }

  process hisat2Align {
    tag "$prefix"
    label 'hisat2'
    label 'highCpu'
    label 'highMem'
    publishDir "${params.outDir}/mapping", mode: 'copy',
      saveAs: {filename ->
        if (filename.indexOf(".hisat2_summary.txt") > 0) "logs/$filename"
        else if (params.saveAlignedIntermediates) filename
        else null
      }

    input:
    set val(prefix), file(reads) from chHisat2RawReads 
    file hs2Index from chHisat2Index.collect()
    file alignmentSplicesites from chAlignmentSplicesites.collect()
    val parseRes from chStrandedResultsHisat

    output:
    file "${prefix}.bam" into chHisat2Bam
    file "${prefix}.hisat2_summary.txt" into chAlignmentLogs
    file("v_hisat2.txt") into chHisat2Version

    script:
    indexBase = hs2Index[0].toString() - ~/.\d.ht2/
    def rnastrandness = ''
    if (parseRes=='forward'){
        rnastrandness = params.singleEnd ? '--rna-strandness F' : '--rna-strandness FR'
    } else if (parseRes=='reverse'){
        rnastrandness = params.singleEnd ? '--rna-strandness R' : '--rna-strandness RF'
    }
    inputOpts = params.singleEnd ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    hisat2 --version &> v_hisat2.txt
    hisat2 -x $indexBase \\
           ${inputOpts} \\
           $rnastrandness \\
           --known-splicesite-infile $alignmentSplicesites \\
           -p ${task.cpus} \\
           --met-stderr \\
           --new-summary \\
	   --rg-id ${prefix} \\
           --summary-file ${prefix}.hisat2_summary.txt \\
           | samtools view -bS -F 4 -F 256 - > ${prefix}.bam
    """
  }

  process hisat2Sort {
    tag "${hisat2Bam.baseName}"
    label 'samtools'
    label 'medCpu'
    label 'medMem'  
    publishDir "${params.outDir}/mapping", mode: 'copy'

    input:
    file hisat2Bam from chHisat2Bam

    output:
    file ('*sorted.{bam,bam.bai}') into chBamBigwig, chBamCount, chBamPreseq, chBamMarkduplicates, 
                                        chBamFeaturecounts, chBamGenetype, chBamHTseqCounts, 
                                        chBamQualimap
    file "${hisat2Bam.baseName}_sorted.bam.bai"
    file("v_samtools.txt") into chSamtoolsVersionSort 

    script:
    def availMem = task.memory ? "-m ${task.memory.toBytes() / task.cpus}" : ''
    """
    samtools --version &> v_samtools.txt
    samtools sort \\
             ${hisat2Bam} \\
             -@ ${task.cpus} $availMem \\
             -m ${params.sortMaxMemory} \\
             -o ${hisat2Bam.baseName}_sorted.bam
    samtools index ${hisat2Bam.baseName}.sorted.bam
    """
  }
}

/*
 * Generate bigwig file
 */

process bigWig {
  tag "${prefix}"
  label 'deeptools'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outDir}/bigWig", mode: "copy",
    saveAs: {filename ->
    	     if ( filename.endsWith(".bigwig") ) "$filename"
             else null}

  when:
  !params.skipBigwig

  input:
  file(bam) from chBamBigwig
  val parseRes from chStrandedResultsBigwig

  output:
  file('*.bigwig') into chBigWig
  file("v_deeptools.txt") into chDeeptoolsVersion

  script:
  prefix = bam[0].toString() - ~/(_sorted)?(.bam)?$/
  strandOpt = parseRes == 'forward' ? '--filterRNAstrand forward' : parseRes == 'reverse' ? '--filterRNAstrand reverse' : ''
  """
  bamCoverage --version &> v_deeptools.txt
  bamCoverage -b ${bam[0]} \\
              -o ${prefix}_cpm.bigwig \\
              -p ${task.cpus} \\
              ${strandOpt} \\
	      --normalizeUsing CPM \\
	      --skipNonCoveredRegions
  """
}

/* 
 * Qualimap
 */

process qualimap {
  tag "${bam[0].baseName}"
  label 'qualimap'
  label 'minCpu'
  label 'medMem'
  publishDir "${params.outDir}/qualimap/" , mode: 'copy'

  when:
  !params.skipQualimap

  input:
  file bam from chBamQualimap
  file gtf from chGtfQualimap.collect()
  val stranded from chStrandedResultsQualimap

  output:
  file ("${bam[0].baseName}") into chQualimapResults
  file ("v_qualimap.txt") into chQualimapVersion

  script:
  peOpts = params.singleEnd ? '' : '-pe'
  memory     = task.memory.toGiga() + "G"
  strandnessOpts = 'non-strand-specific'
  if (stranded == 'forward') {
    strandnessOpts = 'strand-specific-forward'
  } else if (stranded == 'reverse') {
    strandnessOpts = 'strand-specific-reverse'
  }
  """
  mkdir tmp
  export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
  qualimap \\
    --java-mem-size=$memory \\
    rnaseq \\
    -bam $bam[0] \\
    -gtf $gtf \\
    -p $strandnessOpts \\
    $peOpts \\
    -outdir ${bam[0].baseName}
  echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//' > v_qualimap.txt
  """
}


/*
 * Saturation Curves
 */

process preseq {
  tag "${bam[0]}"
  label 'preseq'
  label 'lowCpu'
  label 'medMem'
  publishDir "${params.outDir}/preseq", mode: 'copy'

  when:
  !params.skipQC && !params.skipSaturation

  input:
  file bam from chBamPreseq

  output:
  file "*ccurve.txt" into chPreseqResults
  file("v_preseq.txt") into chPreseqVersion

  script:
  peOpts = params.singleEnd ? '' : '-pe'
  """
  preseq &> v_preseq.txt
  preseq lc_extrap -seed 1 -v -B ${bam[0]} ${peOpts} -o ${bam[0].baseName}_extrap_ccurve.txt -e 200e+06 -seg_len 100000000
  """
}

/*
 * Duplicates
 */

process markDuplicates {
  tag "${bam[0]}"
  label 'picard'
  label 'lowCpu'
  label 'medMem'
  publishDir "${params.outDir}/markDuplicates", mode: 'copy',
    saveAs: {filename -> 
      if (filename.indexOf("_metrics.txt") > 0) "metrics/$filename" 
      else if (params.saveAlignedIntermediates) filename
    }

  when:
  !params.skipQC && !params.skipDupradar

  input:
  file bam from chBamMarkduplicates

  output:
  file('*markDups.{bam,bam.bai}') into ( chBamMd, chBamMdPolym )
  file('*markDups_metrics.txt') into chPicardResults
  file('v_picard.txt') into chPicardVersion

  script:
  markdupMemOption = "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
  """
  echo \$(picard MarkDuplicates --version 2>&1) &> v_picard.txt
  picard ${markdupMemOption} MarkDuplicates \\
      MAX_RECORDS_IN_RAM=50000 \\
      INPUT=${bam[0]} \\
      OUTPUT=${bam[0].baseName}.markDups.bam \\
      METRICS_FILE=${bam[0].baseName}.markDups_metrics.txt \\
      REMOVE_DUPLICATES=false \\
      ASSUME_SORTED=true \\
      PROGRAM_RECORD_ID='null' \\
      VALIDATION_STRINGENCY=LENIENT
  samtools index ${bam[0].baseName}.markDups.bam
  """
}

process dupradar {
  tag "${bamMd[0]}"
  label 'dupradar'
  label 'lowCpu'
  label 'lowMem'
  publishDir "${params.outDir}/dupradar", mode: 'copy',
    saveAs: {filename ->
      if (filename.indexOf("_duprateExpDens.pdf") > 0) "scatter_plots/$filename"
      else if (filename.indexOf("_duprateExpBoxplot.pdf") > 0) "box_plots/$filename"
      else if (filename.indexOf("_expressionHist.pdf") > 0) "histograms/$filename"
      else if (filename.indexOf("_dupMatrix.txt") > 0) "gene_data/$filename"
      else if (filename.indexOf("_duprateExpDensCurve.txt") > 0) "scatter_curve_data/$filename"
      else if (filename.indexOf("_intercept_slope.txt") > 0) "intercepts_slopes/$filename"
      else "$filename"
    }

  when:
  !params.skipQC && !params.skipDupradar

  input:
  file bamMd from chBamMd
  file gtf from chGtfDupradar.collect()
  val parseRes from chStrandedResultsDupradar

  output:
  file "*.{pdf,txt}" into chDupradarResults

  script: 
  def dupradarDirection = 0
  if (parseRes == 'forward'){
      dupradarDirection = 1
  } else if ((parseRes == 'reverse')){
      dupradarDirection = 2
  }
  def paired = params.singleEnd ? 'single' :  'paired'
  """
  dupRadar.r ${bamMd[0]} ${gtf} ${dupradarDirection} ${paired} ${task.cpus}
  """
}

/* 
 * Identito - polym
 */

process getPolym {
  label 'lowCpu'
  label 'medMem'
  label 'identito'

  publishDir "${params.outDir}/identito", mode: 'copy'

  when:
  !params.skipQC && !params.skipIdentito

  input:
  file(fasta) from chFasta.collect()
  file(polyms) from chPolymBed.collect()
  file(bam) from chBamMdPolym

  output:
  file("v_bcftools.txt") into chBcftoolsVersion
  file("*matrix.tsv") into clustPolymCh

  script:
  """
  bcftools --version &> v_bcftools.txt 2>&1 || true
  bcftools mpileup -R ${polyms} -f ${fasta} -x -A -B -q 20 -I -Q 0 -d 1000 --annotate FORMAT/DP,FORMAT/AD ${bam[0]} > ${bam[0].baseName}_bcftools.vcf
  SnpSift extractFields -e "."  -s ";" ${bam[0].baseName}_bcftools.vcf CHROM POS REF ALT GEN[*].DP GEN[*].AD > ${bam[0].baseName}_bcftools.tsv
  computePolym.R ${bam[0].baseName}_bcftools.tsv ${bam[0].baseName}_matrix.tsv ${bam[0].baseName} ${polyms}
  """
}

process combinePolym {
  label 'lowCpu'
  label 'lowMem'
  label 'identito'

  publishDir "${params.outDir}/identito", mode: 'copy'

  when:
  !params.skipQC && !params.skipIdentito

  input:
  file(matrix) from clustPolymCh.collect()

  output:
  file("*.csv") into clustPolymResultsCh
  file("v_R.txt") into chCombinePolymVersion

  script:
  """
  R --version &> v_R.txt
  (head -1 "${matrix[0]}"; tail -n +2 -q *matrix.tsv) > clust_mat.tsv
  computeClust.R clust_mat.tsv ./ 10
  """
}


/*
 * Counts
 */

process featureCounts {
  tag "${bamFeaturecounts.baseName - 'Aligned.sortedByCoord.out'}"
  label 'featurecounts'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outDir}/counts", mode: 'copy',
    saveAs: {filename ->
      if (filename.indexOf("_counts.csv.summary") > 0) "gene_count_summaries/$filename"
      else if (filename.indexOf("_counts.csv") > 0) "gene_counts/$filename"
      else "$filename"
   }

  when:
  params.counts == 'featureCounts'

  input:
  file bam from chBamFeaturecounts
  file gtf from chGtfFeatureCounts.collect()
  val parseRes from chStrandedResultsFeatureCounts

  output:
  file "${bam.baseName}_counts.csv" into chFeatureCountsCountsToMerge, chFeatureCountsCountsToR
  file "${bam.baseName}_counts.csv.summary" into chFeatureCountsLogs
  file("v_featurecounts.txt") into chFeaturecountsVersion

  script:
  def featureCountsDirection = 0
  if (parseRes == 'forward'){
      featureCountsDirection = 1
  } else if ((parseRes == 'reverse')){
      featureCountsDirection = 2
  }
  """
  featureCounts -v &> v_featurecounts.txt
  featureCounts ${params.featurecountsOpts} -T ${task.cpus} -a ${gtf} -o ${bam[0].baseName}_counts.csv -p -s ${featureCountsDirection} ${bamFeaturecounts}
  """
}

process HTseqCounts {
  tag "${bamHTseqCounts}"
  label 'htseq'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outDir}/counts", mode: 'copy',
    saveAs: {filename ->
      if (filename.indexOf("_gene.HTseqCounts.txt.summary") > 0) "gene_count_summaries/$filename"
      else if (filename.indexOf("_gene.HTseqCounts.txt") > 0) "gene_counts/$filename"
      else "$filename"
    }
  
  when:
  params.counts == 'HTseqCounts'

  input:
  file bam from chBamHTseqCounts
  file gtf from chGtfHTseqCounts.collect()
  val parseRes from  chStrandedResultsHTseqCounts

  output: 
  file "*_counts.csv" into chHtseqCountsToMerge, chHtseqCountsToR, chHtseqCountsLogs
  file("v_htseq.txt") into chHtseqVersion 

  script:
  def strandedOpt = '-s no' 
  if (parseRes == 'forward'){
      strandedOpt= '-s yes'
  } else if ((parseRes == 'reverse')){
      strandedOpt= '-s reverse'
  }
  """
  htseq-count -h | grep version  &> v_htseq.txt
  htseq-count ${params.htseqOpts} ${strandedOpt} ${bam[0]} $gtf > ${bam[0].baseName}_counts.csv
  """
}

chCountsToMerge = Channel.empty()
chCountsToR = Channel.empty()
if( params.counts == 'featureCounts' ){
  chCountsToMerge = chFeatureCountsCountsToMerge
  chCountsToR = chFeatureCountsCountsToR
} else if (params.counts == 'HTseqCounts'){
  chCountsToMerge = chHtseqCountsToMerge
  chCountsToR = chHtseqCountsToR	
}else if (params.counts == 'star'){
  chCountsToMerge = chStarCountsToMerge
  chCountsToR = chStarCountsToR
}

process mergeCounts {
  publishDir "${params.outDir}/counts", mode: 'copy'
  label 'r'
  label 'lowCpu'
  label 'medMem'

  input:
  file inputCounts from chCountsToMerge.collect()
  file gtf from chGtfTable.collect()
  val parseRes from chStrandedResultsTable.collect()

  output:
  file 'tablecounts_raw.csv' into chRawCounts, chCountsSaturation
  file 'tablecounts_tpm.csv' into chTpmCounts, chTpmGenetype
  file 'tableannot.csv' into chGenesAnnot
  file("v_R.txt") into chMergeCountsVersion

  script:
  """
  R --version &> v_R.txt
  echo -e ${inputCounts} | tr " " "\n" > listofcounts.tsv
  echo -n "${parseRes}" | sed -e "s/\\[//" -e "s/\\]//" -e "s/,//g" | tr " " "\n" > listofstrandness.tsv
  makeCountTable.r listofcounts.tsv ${gtf} ${params.counts} listofstrandness.tsv
  """
}

chCountsLogs = Channel.empty()
if( params.counts == 'featureCounts' ){
  chCountsLogs = chFeatureCountsLogs
} else if (params.counts == 'HTseqCounts'){
  chCountsLogs = chHtseqCountsLogs
}else if (params.counts == 'star'){
  chCountsLogs = chStarLogCounts
}


/*
 * Gene-based saturation
 */

process geneSaturation {
  label 'r'
  label 'lowCpu'
  label 'medMem'
  publishDir "${params.outDir}/geneSaturation" , mode: 'copy'

  when:
  !params.skipQC && !params.skipSaturation

  input:
  file inputCounts from chCountsSaturation.collect()

  output:
  file "*gcurve.txt" into chGenesatResults
  file("v_R.txt") into chGeneSaturationVersion

  script:
  """
  R --version &> v_R.txt
  gene_saturation.r $inputCounts counts.gcurve.txt
  """
}


/*
 * Reads distribution
 */

process getCountsPerGeneType {
  label 'r'
  label 'lowCpu'
  label 'lowMem'
  publishDir "${params.outDir}/readDistribution", mode: 'copy'

  when:
  !params.skipReaddist

  input:
  file tpmGenetype from chTpmGenetype
  file gtf from chGtfGenetype.collect()
 
  output:
  file "*genetype.txt" into chCountsPerGenetype
  file("v_R.txt") into chGeneTypeVersion

  script:
  """
  R --version &> v_R.txt
  gene_type_expression.r ${tpmGenetype} ${gtf} counts_genetype.txt 
  """
}


/*
 * Exploratory analysis
 */

process exploratoryAnalysis {
  label 'r'
  label 'lowCpu'
  label 'lowMem'
  publishDir "${params.outDir}/exploratoryAnalysis", mode: 'copy'

  when:
  !params.skipExpan && numSample > 1

  input:
  file tableRaw from chRawCounts.collect()
  file tableTpm from chTpmCounts.collect()
  val numSample from chCountsToR.count()
  file pcaHeader from chPcaHeader
  file heatmapHeader from chHeatmapHeader

  output:
  file "*.{txt,pdf,csv}" into chExploratoryAnalysisResults
  file("v_R.txt") into chAnaExpVersion

  script:
  """
  R --version &> v_R.txt
  exploratory_analysis.r ${tableRaw}
  cat $pcaHeader deseq2_pca_coords_mqc.csv >> tmp_file
  mv tmp_file deseq2_pca_coords_mqc.csv 
  cat $heatmapHeader vst_sample_cor_mqc.csv >> tmp_file
  mv tmp_file vst_sample_cor_mqc.csv
  """
}

/*
 * MultiQC
 */

process getSoftwareVersions{
  label 'python'
  label 'lowCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/softwareVersions", mode: "copy"

  when:
  !params.skipSoftVersions

  input:
  file 'v_fastqc.txt' from chFastqcVersion.first().ifEmpty([])
  file 'v_star.txt' from chStarVersion.first().ifEmpty([])
  file 'v_hisat2.txt' from chHisat2Version.first().ifEmpty([])
  file 'v_bowtie.txt' from chBowtieVersion.first().ifEmpty([])
  file 'v_bowtie2.txt' from chBowtie2Version.first().ifEmpty([])
  file 'v_samtools.txt' from chSamtoolsVersionSort.first().ifEmpty([])
  file 'v_picard.txt' from chPicardVersion.first().ifEmpty([])
  file 'v_preseq.txt' from chPreseqVersion.first().ifEmpty([])
  file 'v_R.txt' from chMergeCountsVersion.concat(chCombinePolymVersion,chGeneSaturationVersion,chGeneTypeVersion,chAnaExpVersion).first().ifEmpty([])
  file 'v_rseqc.txt' from chRseqcVersionInferExperiment.first().ifEmpty([])
  file 'v_featurecounts.txt' from chFeaturecountsVersion.first().ifEmpty([])
  file 'v_deeptools.txt' from chDeeptoolsVersion.first().ifEmpty([])
  file 'v_bcftools.txt' from chBcftoolsVersion.first().ifEmpty([])
  file 'v_htseq.txt' from chHtseqVersion.first().ifEmpty([])
  file 'v_qualimap.txt' from chQualimapVersion.first().ifEmpty([])

  output:
  file 'software_versions_mqc.yaml' into softwareVersionsYaml

  script:
  """
  echo $workflow.manifest.version &> v_pipeline.txt
  echo $workflow.nextflow.version &> v_nextflow.txt
  scrape_software_versions.py &> software_versions_mqc.yaml
  """
}

process workflowSummaryMqc {
  label 'unix'
  label 'minCpu'
  label 'minMem'

  when:
  !params.skipMultiQC

  output:
  file 'workflow_summary_mqc.yaml' into workflowSummaryYaml

  exec:
  def yaml_file = task.workDir.resolve('workflow_summary_mqc.yaml')
  yaml_file.text  = """
  id: 'summary'
  description: " - this information is collected when the pipeline is started."
  section_name: 'Workflow Summary'
  section_href: 'https://gitlab.curie.fr/data-analysis/RNA-seq'
  plot_type: 'html'
  data: |
      <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
      </dl>
  """.stripIndent()
}


if (skippedPoorAlignment.size() > 0){
  Channel.fromList(skippedPoorAlignment)
         .flatMap{ it -> it + ": Poor alignment rate. Sample discarded"}
         .collectFile(name: 'warnings.txt', newLine: true)
         .set{chWarn}
}else{
  chWarn = Channel.empty()
}

process multiqc {
  label 'multiqc'
  label 'lowCpu'
  label 'lowMem'
  publishDir "${params.outDir}/MultiQC", mode: 'copy'

  when:
  !params.skipMultiQC

  input:
  file splan from chSplan.collect()
  file metadata from chMetadata.ifEmpty([])
  file multiqcConfig from chMultiqcConfig.ifEmpty([])
  file (fastqc:'fastqc/*') from chFastqcResults.collect().ifEmpty([])
  file ('rrna/*') from chRrnaLogs.collect().ifEmpty([])
  file ('alignment/*') from chAlignmentLogs.collect().ifEmpty([])
  file ('strandness/*') from chStrandnessResults.collect().ifEmpty([])
  file ('qualimap/*') from chQualimapResults.collect().ifEmpty([])
  file ('preseq/*') from chPreseqResults.collect().ifEmpty([])
  file ('genesat/*') from chGenesatResults.collect().ifEmpty([])
  file ('dupradar/*') from chDupradarResults.collect().ifEmpty([])
  file ('picard/*') from chPicardResults.collect().ifEmpty([])	
  file ('counts/*') from chCountsLogs.collect().ifEmpty([])
  file ('genetype/*') from chCountsPerGenetype.collect().ifEmpty([])
  file ('identito/*') from clustPolymResultsCh.collect().ifEmpty([])
  file ('exploratoryAnalysis_results/*') from chExploratoryAnalysisResults.collect().ifEmpty([]) 
  file ('softwareVersions/*') from softwareVersionsYaml.collect().ifEmpty([])
  file ('workflowSummary/*') from workflowSummaryYaml.collect().ifEmpty([])
  file ('workflowSummary/*') from chWarn.collect().ifEmpty([]) 

  output:
  file splan
  file "*report.html" into multiqcReport
  file "*_data"

  script:
  rtitle = customRunName ? "--title \"$customRunName\"" : ''
  rfilename = customRunName ? "--filename " + customRunName + "_rnaseq_report" : "--filename rnaseq_report"
  metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
  splanOpts = params.samplePlan ? "--splan ${params.samplePlan}" : ""
  isPE = params.singleEnd ? 0 : 1
    
  modulesList = "-m custom_content -m preseq -m rseqc -m bowtie1 -m hisat2 -m star -m cutadapt -m fastqc -m qualimap"
  modulesList = params.counts == 'featureCounts' ? "${modulesList} -m featureCounts" : "${modulesList}"  
  modulesList = params.counts == 'HTseqCounts' ? "${modulesList} -m htseq" : "${modulesList}"  
 
  warn=skippedPoorAlignment.size() > 0 ? "--warn workflowSummary/warnings.txt" : ""
  """
  stats2multiqc.sh ${splan} ${params.aligner} ${isPE}
  medianReadNb="\$(sort -t, -k3,3n mq.stats | awk -F, '{a[i++]=\$3;} END{x=int((i+1)/2); if (x<(i+1)/2) printf "%.0f", (a[x-1]+a[x])/2; else printf "%.0f",a[x-1];}')"
  mqc_header.py --name "RNA-seq" --version ${workflow.manifest.version} ${metadataOpts} ${splanOpts} --nbreads \${medianReadNb} ${warn} > multiqc-config-header.yaml
  multiqc . -f $rtitle $rfilename -c $multiqcConfig -c multiqc-config-header.yaml $modulesList
  """    
}


/*
 * Sub-routine
 */


process outputDocumentation {
  label 'python'
  label 'lowCpu'
  label 'lowMem'
  publishDir "${params.summaryDir}/", mode: 'copy'

  input:
  file outputDocs from chOutputDocs
  file images from chOutputDocsImages

  output:
  file "results_description.html"

  script:
  """
  markdown_to_html.py $outputDocs -o results_description.html
  """
}                                                                                                                                                                                                           


workflow.onComplete {

  /*pipeline_report.html*/

  def report_fields = [:]
  report_fields['version'] = workflow.manifest.version
  report_fields['runName'] = customRunName ?: workflow.runName
  report_fields['success'] = workflow.success
  report_fields['dateComplete'] = workflow.complete
  report_fields['duration'] = workflow.duration
  report_fields['exitStatus'] = workflow.exitStatus
  report_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
  report_fields['errorReport'] = (workflow.errorReport ?: 'None')
  report_fields['commandLine'] = workflow.commandLine
  report_fields['projectDir'] = workflow.projectDir
  report_fields['summary'] = summary
  report_fields['summary']['Date Started'] = workflow.start
  report_fields['summary']['Date Completed'] = workflow.complete
  report_fields['summary']['Pipeline script file path'] = workflow.scriptFile
  report_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
  if(workflow.repository) report_fields['summary']['Pipeline repository Git URL'] = workflow.repository
  if(workflow.commitId) report_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
  if(workflow.revision) report_fields['summary']['Pipeline Git branch/tag'] = workflow.revision

  report_fields['skippedPoorAlignment'] = skippedPoorAlignment

  // Render the TXT template
  def engine = new groovy.text.GStringTemplateEngine()
  def tf = new File("$baseDir/assets/onCompleteTemplate.txt")
  def txt_template = engine.createTemplate(tf).make(report_fields)
  def report_txt = txt_template.toString()
    
  // Render the HTML template
  def hf = new File("$baseDir/assets/onCompleteTemplate.html")
  def html_template = engine.createTemplate(hf).make(report_fields)
  def report_html = html_template.toString()
  // Write summary e-mail HTML to a file
  def output_d = new File( "${params.summaryDir}/" )
  if( !output_d.exists() ) {
    output_d.mkdirs()
  }
  def output_hf = new File( output_d, "pipelineReport.html" )
  output_hf.withWriter { w -> w << report_html }
  def output_tf = new File( output_d, "pipelineReport.txt" )
  output_tf.withWriter { w -> w << report_txt }

  /*oncomplete file*/
  File woc = new File("${params.outDir}/workflowOnComplete.txt")
  Map endSummary = [:]
  endSummary['Completed on'] = workflow.complete
  endSummary['Duration']     = workflow.duration
  endSummary['Success']      = workflow.success
  endSummary['exit status']  = workflow.exitStatus
  endSummary['Error report'] = workflow.errorReport ?: '-'
  String endWfSummary = endSummary.collect { k,v -> "${k.padRight(30, '.')}: $v" }.join("\n")
  println endWfSummary
  String execInfo = "Execution summary\n${endWfSummary}\n"
  woc.write(execInfo)

  /*final logs*/
  if(skippedPoorAlignment.size() > 0){
    log.info "[rnaseq] WARNING - ${skippedPoorAlignment.size()} samples skipped due to poor alignment scores!"
  }
  if(workflow.success){
    log.info "[rnaseq] Pipeline Complete"
  }else{
    log.info "[rnaseq] FAILED: $workflow.runName"
  } 
}
