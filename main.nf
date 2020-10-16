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
    --aligner [str]                      Tool for read alignments ['star', 'hisat2', 'tophat2']. Default: 'star'

  Counts:
    --counts [str]                       Tool to use to estimate the raw counts per gene ['star', 'featureCounts', 'HTseqCounts']. Default: 'star'

  References: If not specified in the configuration file or you wish to overwrite any of the references.
    --genomeAnnotationPath [file]        Path  to genome annotation folder
    --starIndex [dir]                    Path to STAR index
    --hisat2Index [file]                 Path to HiSAT2 index
    --tophat2Index [file]                Path to TopHat2 index
    --gtf [file]                         Path to GTF file
    --bed12 [file]                       Path to gene bed12 file
    --saveAlignedIntermediates [bool]    Save the intermediate files from the Aligment step  - not done by default

  Other options:
    --metadata [file]                    Add metadata file for multiQC report
    --outdir [dir]                       The output directory where the results will be saved
    -w/--work-dir [dir]                  The temporary directory where intermediate data will be saved
    -name [str]                          Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

  Skip options:
    --skipQC [bool]                      Skip all QC steps apart from MultiQC
    --skipRrna [bool]                    Skip rRNA mapping
    --skipFastqc [bool]                  Skip FastQC
    --skipGenebodyCoverage [bool]        Skip calculating genebody coverage
    --skipSaturation [bool]              Skip Saturation qc
    --skipDupradar [bool]                Skip dupRadar (and Picard MarkDups)
    --skipReaddist [bool]                Skip read distribution steps
    --skipExpan [bool]                   Skip exploratory analysis
    --skipMultiQC [bool]                 Skip MultiQC

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
if (params.aligner != 'star' && params.aligner != 'hisat2' && params.aligner != 'tophat2'){
  exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star', 'hisat2', 'tophat2'"
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
else if ( params.bowtie2Index && params.aligner == 'tophat2' ){
  Channel.fromPath("${params.bowtie2Index}*")
    .ifEmpty { exit 1, "TOPHAT2 index not found: ${params.bowtie2Index}" }
    .set {chTophat2Index}
}
else {
    exit 1, "No reference genome specified!"
}

if( params.gtf ){
  Channel
    .fromPath(params.gtf)
    .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
    .into { chGtfStar; chGtfDupradar; chGtfFeatureCounts; chGtfGenetype; chGtfHTseqCounts; chGtfTophat; chGtfTable; chGtfMakeHisatSplicesites }
}else {
  log.warn "No GTF annotation specified - dupRadar, table counts - will be skipped !" 
  Channel
    .empty()
    .into { chGtfStar; chGtfDupradar; chGtfFeatureCounts; chGtfGenetype; chGtfHTseqCounts; chGtfTophat; chGtfTable; chGtfMakeHisatSplicesites } 
}

if( params.bed12 ){
  Channel
    .fromPath(params.bed12)
    .ifEmpty { exit 1, "BED12 annotation file not found: ${params.bed12}" }
    .into { chBedRseqc; chBedReadDist; chBedGenebodyCoverage} 
}else{
  log.warn "No BED gene annotation specified - strandness detection, gene body coverage, read distribution - will be skipped !"
  Channel
    .empty()
    .into { chBedRseqc; chBedReadDist; chBedGenebodyCoverage}
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

/*
 * Create a channel for input read files
 */

if(params.samplePlan){
  if(params.singleEnd){
    Channel
      .from(file("${params.samplePlan}"))
      .splitCsv(header: false)
      .map{ row -> [ row[0], [file(row[2])]] }
      .into { chRawReadsFastqc; chRawReadsStar; chRawReadsHisat2; chRawReadsTophat2; chRawReadsRnaMapping; chRawReadsPrepRseqc; chRawReadsStrandness; chSaveStrandness }
  }else{
     Channel
       .from(file("${params.samplePlan}"))
       .splitCsv(header: false)
       .map{ row -> [ row[0], [file(row[2]), file(row[3])]] }
       .into { chRawReadsFastqc; chRawReadsStar; chRawReadsHisat2; chRawReadsTophat2; chRawReadsRnaMapping; chRawReadsPrepRseqc; chRawReadsStrandness; chSaveStrandness }
   }
   params.reads=false
}
else if(params.readPaths){
  if(params.singleEnd){
    Channel
      .from(params.readPaths)
      .map { row -> [ row[0], [file(row[1][0])]] }
      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
      .into { chRawReadsFastqc; chRawReadsStar; chRawReadsHisat2; chRawReadsTophat2; chRawReadsRnaMapping; chRawReadsPrepRseqc; chRawReadsStrandness; chSaveStrandness }
  } else {
    Channel
      .from(params.readPaths)
      .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
      .into { chRawReadsFastqc; chRawReadsStar; chRawReadsHisat2; chRawReadsTophat2; chRawReadsRnaMapping; chRawReadsPrepRseqc; chRawReadsStrandness; chSaveStrandness }
  }
} else {
  Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { chRawReadsFastqc; chRawReadsStar; chRawReadsHisat2; chRawReadsTophat2; chRawReadsRnaMapping; chRawReadsPrepRseqc; chRawReadsStrandness; chSaveStrandness }
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
} else if(params.aligner == 'tophat2') {
  summary['Aligner'] = "Tophat2"
  if(params.bowtie2Index) summary['Tophat2 Index'] = params.bowtie2Index
} else if(params.aligner == 'hisat2') {
  summary['Aligner'] = "HISAT2"
  if(params.hisat2Index) summary['HISAT2 Index'] = params.hisat2Index
}
summary['Counts'] = params.counts
if(params.gtf)  summary['GTF Annotation']  = params.gtf
if(params.bed12) summary['BED Annotation']  = params.bed12
summary['Save Intermeds'] = params.saveAlignedIntermediates ? 'Yes' : 'No'
summary['Max Memory']     = params.maxMemory
summary['Max CPUs']       = params.maxCpus
summary['Max Time']       = params.maxTime
summary['Container Engine'] = workflow.containerEngine
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


/*
 * FastQC
 */
process fastqc {
  tag "${prefix}"
  publishDir "${params.outdir}/fastqc", mode: 'copy',
    saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

  when:
  !params.skipQC && !params.skipFastqc

  input:
  set val(prefix), file(reads) from chRawReadsFastqc

  output:
  file "*_fastqc.{zip,html}" into chFastqcResults

  script:
  pbase = reads[0].toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
  """
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
  publishDir "${params.outdir}/rRNA_mapping", mode: 'copy',
    saveAs: {filename ->
      if (filename.indexOf("fastq.gz") > 0 &&  params.saveAlignedIntermediates) filename
      else if (filename.indexOf(".log") > 0) "logs/$filename"
      else null
    }

  when:
  !params.skipRrna && paramsRrna

  input:
  set val(prefix), file(reads) from chRawReadsRnaMapping
  file annot from chRrnaAnnot.collect()

  output:
  set val(prefix), file("*fastq.gz") into chRrnaMappingRes
  set val(prefix), file("*.sam") into chRrnaSam
  file "*.log" into chRrnaLogs

  script:
  if (params.singleEnd) {
  """
  bowtie ${params.bowtieOpts} \\
         -p ${task.cpus} \\
         --un ${prefix}_norRNA.fastq \\
         --sam ${params.rrna} \\
         ${reads} \\
         ${prefix}.sam  2> ${prefix}.log && \
  gzip -f ${prefix}_norRNA*.fastq 
 """
  } else {
  """
  bowtie ${params.bowtieOpts} \\
         -p ${task.cpus} \\
         --un ${prefix}_norRNA.fastq \\
         --sam ${params.rrna} \\
         -1 ${reads[0]} \\
         -2 ${reads[1]} \\
         ${prefix}.sam  2> ${prefix}.log && \
  gzip -f ${prefix}_norRNA_*.fastq 
  """
  }  
}


/*
 * Strandness
 */

// User defined

if (params.stranded == 'reverse' || params.stranded == 'forward' || params.stranded == 'no'){
  chRawReadsStrandness
    .map { file ->
           def key = params.stranded
           return tuple(key)
    }
    .into { chStrandedResultsFeatureCounts; chStrandedResultsGenetype; chStrandedResultsHTseqCounts;
            chStrandedResultsDupradar; chStrandedResultsTophat; chStrandedResultsHisat; chStrandedResultsTable }
}else if (params.stranded == 'auto' && !params.bed12){
  log.warn "Strandness ('auto') cannot be run without GTF annotation - will be skipped !"
}

process saveStrandness {
  publishDir "${params.outdir}/strandness" , mode: 'copy',
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
  if (params.stranded == 'auto' && !params.bed12){
  """
  echo "NA" > ${prefix}_strandness.txt
  """
  }else{
  """
  echo ${params.stranded} > ${prefix}_strandness.txt
  """
  }
}

// auto

process prepRseqc {
  tag "${prefix}"

  input:
  set val(prefix), file(reads) from chRawReadsPrepRseqc

  when:
  params.stranded == 'auto' && params.bed12

  output:
  set val("${prefix}"), file("${prefix}_subsample.bam") into chBamRseqc

  script:
  if (params.singleEnd) {
  """
  bowtie2 --fast --end-to-end --reorder \\
          -p ${task.cpus} \\
          -u ${params.nCheck} \\
          -x ${params.bowtie2Index} \\
          -U ${reads} > ${prefix}_subsample.bam 
   """
  } else {
  """
  bowtie2 --fast --end-to-end --reorder \\
          -p ${task.cpus} \\
          -u ${params.nCheck} \\
          -x ${params.bowtie2Index} \\
          -1 ${reads[0]} \\
          -2 ${reads[1]} > ${prefix}_subsample.bam
  """
  }
}

process rseqc {
  tag "${prefix - '_subsample'}"
  publishDir "${params.outdir}/strandness" , mode: 'copy',
    saveAs: {filename ->
      if (filename.indexOf(".txt") > 0) "$filename"
      else null
    }

  input:
  set val(prefix), file(bamRseqc) from chBamRseqc
  file bed12 from chBedRseqc.collect()

  output:
  file "*.{txt,pdf,r,xls}" into chRseqcResults
  stdout into ( chStrandedResultsFeatureCounts, chStrandedResultsGenetype, chStrandedResultsHTseqCounts,
                chStrandedResultsDupradar, chStrandedResultsTophat, chStrandedResultsHisat, chStrandedResultsTable )

  when:
  params.stranded == 'auto' && params.bed12
  
  script:
  """
  infer_experiment.py -i $bamRseqc -r $bed12 > ${prefix}.txt
  parse_rseq_output.sh ${prefix}.txt > ${prefix}_strandness.txt
  cat ${prefix}_strandness.txt
  """  
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
def checkLog(logs) {
  def percentAligned = 0;
  logs.eachLine { line ->
    if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
      percentAligned = matcher[0][1]
    }
  }
  logname = logs.getBaseName() - 'Log.final'
  if(percentAligned.toFloat() <= '2'.toFloat() ){
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
  hisatStdout = Channel.from(false)
  process star {
    tag "$prefix"
    publishDir "${params.outdir}/mapping", mode: 'copy',
      saveAs: {filename ->
        if (filename.indexOf(".bam") == -1) "logs/$filename"
        else if (params.saveAlignedIntermediates) filename
        else null
      }
    publishDir "${params.outdir}/counts", mode: 'copy',
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

    script:
    def starCountOpt = params.counts == 'star' && params.gtf ? params.starOptsCounts : ''
    def starGtfOpt = params.gtf ? "--sjdbGTFfile $gtf" : ''
    """
    STAR --genomeDir $index \\
         ${starGtfOpt} \\
         --readFilesIn $reads  \\
         --runThreadN ${task.cpus} \\
         --runMode alignReads \\
         --outSAMtype BAM Unsorted  \\
         --readFilesCommand zcat \\
         --runDirPerm All_RWX \\
         --outTmpDir /local/scratch/rnaseq_\$(date +%d%s%S%N) \\
         --outFileNamePrefix $prefix  \\
         --outSAMattrRGline ID:$prefix SM:$prefix LB:Illumina PL:Illumina  \\
         ${params.starOptions} \\
	 --limitOutSJcollapsed 5000000 \\
	 ${starCountOpt}
    """
  }

  process starSort {
    tag "$prefix"
    publishDir "${params.outdir}/mapping", mode: 'copy'
 
    input:
    set val(prefix), file(LogFinalOut), file (starBam) from chStarSam

    output:
    set file("${prefix}Log.final.out"), file ('*.bam') into  chStarAligned  
    file "${prefix}Aligned.sortedByCoord.out.bam.bai" into chBamIndexStar

    script:
    """
    samtools sort  \\
        -@  ${task.cpus}  \\
        -m ${params.sortMaxMemory} \\
        -o ${prefix}Aligned.sortedByCoord.out.bam  \\
        ${starBam}   

    samtools index ${prefix}Aligned.sortedByCoord.out.bam
    """
    }

    // Filter removes all 'aligned' channels that fail the check
    chStarAligned
      .filter { logs, bams -> checkLog(logs) }
      .flatMap {  logs, bams -> bams }
      .into { chBamCount; chBamPreseq; chBamMarkduplicates; chBamFeaturecounts; chBamGenetype; chBamHTseqCounts; chBamReadDist; chBamForSubsamp; chBamSkipSubsamp }
}


// HiSat2

chHisat2RawReads = Channel.empty()
if( params.rrna && !params.skipRrna ){
    chHisat2RawReads = chRrnaMappingRes
}else {
    chHisat2RawReads = chRawReadsHisat2 
}

if(params.aligner == 'hisat2'){
  chStarLog = Channel.from(false)
  
  process makeHisatSplicesites {
     tag "$gtf"
     publishDir "${params.outdir}/mapping", mode: 'copy',
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
     hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.hisat2_splice_sites.txt
     """
  }

  process hisat2Align {
    tag "$prefix"
    publishDir "${params.outdir}/mapping", mode: 'copy',
      saveAs: {filename ->
        if (filename.indexOf(".hisat2_summary.txt") > 0) "logs/$filename"
        else if (params.saveAlignedIntermediates) filename
        else null
      }

    input:
    set val(prefix), file(reads) from chHisat2RawReads 
    file hs2Index from chHisat2Indices.collect()
    file alignmentSplicesites from chAlignmentSplicesites.collect()
    val parseRes from chStrandedResultsHisat

    output:
    file "${prefix}.bam" into chHisat2Bam
    file "${prefix}.hisat2_summary.txt" into chAlignmentLogs

    script:
    indexBase = hs2Index[0].toString() - ~/.\d.ht2/
    seqCenter = params.seqCenter ? "--rg-id ${prefix} --rg CN:${params.seqCenter.replaceAll('\\s','_')}" : ''
    def rnastrandness = ''
    if (parseRes=='forward'){
        rnastrandness = params.singleEnd ? '--rna-strandness F' : '--rna-strandness FR'
    } else if (parseRes=='reverse'){
        rnastrandness = params.singleEnd ? '--rna-strandness R' : '--rna-strandness RF'
    }
    if (params.singleEnd) {
    """
    hisat2 -x $indexBase \\
           -U $reads \\
           $rnastrandness \\
           --known-splicesite-infile $alignmentSplicesites \\
           -p ${task.cpus} \\
           --met-stderr \\
           --new-summary \\
           --summary-file ${prefix}.hisat2_summary.txt $seqCenter \\
           | samtools view -bS -F 4 -F 256 - > ${prefix}.bam
    """
    } else {
    """
    hisat2 -x $indexBase \\
           -1 ${reads[0]} \\
           -2 ${reads[1]} \\
           $rnastrandness \\
           --known-splicesite-infile $alignmentSplicesites \\
           --no-mixed \\
           --no-discordant \\
           -p ${task.cpus} \\
           --met-stderr \\
           --new-summary \\
           --summary-file ${prefix}.hisat2_summary.txt $seqCenter \\
           | samtools view -bS -F 4 -F 8 -F 256 - > ${prefix}.bam
     """
    }
  }

  process hisat2Sort {
      tag "${hisat2Bam.baseName}"
      publishDir "${params.outdir}/mapping", mode: 'copy'

      input:
      file hisatBam from chHisat2Bam

      output:
      file "${hisatBam.baseName}.sorted.bam" into chBamCount, chBamPreseq, chBamMarkduplicates, bchBmFeaturecounts, chBamGenetype, chBamHTseqCounts, chBamReadDist, chBamForSubsamp, chBamSkipSubsamp
      file "${hisatBam.baseName}.sorted.bam.bai" into chBamIndexHisat
 
      script:
      def availMem = task.memory ? "-m ${task.memory.toBytes() / task.cpus}" : ''
      """
      samtools sort \\
          $hisat2_bam \\
          -@ ${task.cpus} $availMem \\
          -m ${params.sortMaxMemory} \\
          -o ${hisat2Bam.baseName}.sorted.bam
      samtools index ${hisatBam.baseName}.sorted.bam
      """
  }
}

// Update channel for TopHat2
chTophat2RawReads = Channel.empty()
if( params.rrna && !params.skipRrna ){
  chTophat2RawReads = chRrnaMappingRes
}
else {
  chTophat2RawReads = chRawReadsTophat2
}

if(params.aligner == 'tophat2'){
  process tophat2 {
   tag "${prefix}"
   publishDir "${params.outdir}/mapping", mode: 'copy',
     saveAs: {filename ->
       if (filename.indexOf(".align_summary.txt") > 0) "logs/$filename"
       else filename
     }

   input:
   set val(prefix), file(reads) from chTophat2RawReads 
   file "tophat2" from chTophat2Index.collect()
   file gtf from chGtfTophat.collect()
   val parseRes from chStrandedResultsTophat

   output:
   file "${prefix}.bam" into chBamCount, chBamPreseq, chBamMarkduplicates, chBamFeaturecounts, chBamGenetype, chBamHTseqCounts, chBamReadDist, chBamForSubsamp, chBamSkipSubsamp
   file "${prefix}.align_summary.txt" into chAlignmentLogs
   file "${prefix}.bam.bai" into chBamIndexTophat

   script:
   def availMem = task.memory ? "-m ${task.memory.toBytes() / task.cpus}" : ''
   def strandedOpt = '--library-type fr-unstranded'
   if (parseRes == 'forward'){
     strandedOpt = '--library-type fr-secondstrand'
   }else if ((parseRes == 'reverse')){
     strandedOpt = '--library-type fr-firststrand'
   }
   def out = './mapping'
   def sample = "--rg-id ${prefix} --rg-sample ${prefix} --rg-library Illumina --rg-platform Illumina --rg-platform-unit ${prefix}"
   """
   mkdir -p ${out}
   tophat2 -p ${task.cpus} \\
           ${sample} \\
           ${params.tophat2Opts} \\
          --GTF $gtf \\
          ${strandedOpt} \\
          -o ${out} \\
          ${params.bowtie2Index} \\
          ${reads} 

   mv ${out}/accepted_hits.bam ${prefix}.bam
   mv ${out}/align_summary.txt ${prefix}.align_summary.txt
   samtools index ${prefix}.bam
   """
  }
}


/*
 * Subsample the BAM files if necessary
 */
chBamForSubsamp
  .filter { it.size() > params.subsampFilesizeThreshold }
  .map { [it, params.subsampFilesizeThreshold / it.size() ] }
  .set{ chBamForSubsampFiltered }
chBamSkipSubsamp
   .filter { it.size() <= params.subsampFilesizeThreshold }
   .set{ chBamSkipSubsampFiltered }

process bamSubsample {
  tag "${bam.baseName - '.sorted'}"

  input:
  set file(bam), val(fraction) from chBamForSubsampFiltered

  output:
  file "*_subsamp.bam" into chBamSubsampled

  script:
  """
  samtools view -s $fraction -b $bam | samtools sort -o ${bam.baseName}_subsamp.bam
  """
}

/*
 * Rseqc genebody_coverage
 */

process genebodyCoverage {
  tag "${bam.baseName - '.sorted'}"
  publishDir "${params.outdir}/genecov" , mode: 'copy',
     saveAs: {filename ->
       if (filename.indexOf("geneBodyCoverage.curves.pdf") > 0)       "geneBodyCoverage/$filename"
       else if (filename.indexOf("geneBodyCoverage.r") > 0)           "geneBodyCoverage/rscripts/$filename"
       else if (filename.indexOf("geneBodyCoverage.txt") > 0)         "geneBodyCoverage/data/$filename"
       else if (filename.indexOf("log.txt") > -1) false
       else filename
      }

  when:
  !params.skipQC && !params.skipGenebodyCoverage

  input:
  file bam from chBamSubsampled.concat(chBamSkipSubsampFiltered)
  file bed12 from chBedGenebodyCoverage.collect()

  output:
  file "*.{txt,pdf,r}" into chGenebodyCoverageResults

  script:
  """
  samtools index $bam
  geneBody_coverage.py \\
    -i $bam \\
    -o ${bam.baseName}.rseqc \\
    -r $bed12
  mv log.txt ${bam.baseName}.rseqc.log.txt
  """
}

/*
 * Saturation Curves
 */

process preseq {
  tag "${bamPreseq}"
  publishDir "${params.outdir}/preseq", mode: 'copy'

  when:
  !params.skipQC && !params.skipSaturation

  input:
  file bamPreseq from chBamPreseq

  output:
  file "*ccurve.txt" into chPreseqResults

  script:
  prefix = bamPreseq.toString() - ~/(.bam)?$/
  """
  preseq lc_extrap -v -B $bamPreseq -o ${prefix}.extrap_ccurve.txt -e 200e+06
  """
}

/*
 * Duplicates
 */

process markDuplicates {
  tag "${bam}"
  publishDir "${params.outdir}/markDuplicates", mode: 'copy',
    saveAs: {filename -> 
      if (filename.indexOf("_metrics.txt") > 0) "metrics/$filename" 
      else if (params.saveAlignedIntermediates) filename
    }

  when:
  !params.skipQC && !params.skipDupradar

  input:
  file bam from chBamMarkduplicates

  output:
  file "${bam.baseName}.markDups.bam" into chBamMd
  file "${bam.baseName}.markDups_metrics.txt" into chPicardResults
  file "${bam.baseName}.markDups.bam.bai"

  script:
  if( !task.memory ){
    log.info "[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this."
    availMem = 3
  } else {
    availMem = task.memory.toGiga()
  }

  markdupJavaOptions = availMem > 8 ? params.markdupJavaOptions : "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""  
  """
  picard ${markdupJavaOptions} -Djava.io.tmpdir=/local/scratch MarkDuplicates \\
      MAX_RECORDS_IN_RAM=50000 \\
      INPUT=$bam \\
      OUTPUT=${bam.baseName}.markDups.bam \\
      METRICS_FILE=${bam.baseName}.markDups_metrics.txt \\
      REMOVE_DUPLICATES=false \\
      ASSUME_SORTED=true \\
      PROGRAM_RECORD_ID='null' \\
      VALIDATION_STRINGENCY=LENIENT
  samtools index ${bam.baseName}.markDups.bam
  """
}

process dupradar {
  tag "${bamMd}"
  publishDir "${params.outdir}/dupradar", mode: 'copy',
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
  dupRadar.r $bamMd $gtf $dupradarDirection $paired ${task.cpus}
  """
}

/*
 * Counts
 */

process featureCounts {
  tag "${bamFeaturecounts.baseName - 'Aligned.sortedByCoord.out'}"
  publishDir "${params.outdir}/counts", mode: 'copy',
    saveAs: {filename ->
      if (filename.indexOf("_counts.csv.summary") > 0) "gene_count_summaries/$filename"
      else if (filename.indexOf("_counts.csv") > 0) "gene_counts/$filename"
      else "$filename"
   }

  when:
  params.counts == 'featureCounts'

  input:
  file bamFeaturecounts from chBamFeaturecounts
  file gtf from chGtfFeatureCounts.collect()
  val parseRes from chStrandedResultsFeatureCounts

  output:
  file "${bamFeaturecounts.baseName}_counts.csv" into chFeatureCountsCountsToMerge, chFeatureCountsCountsToR
  file "${bamFeaturecounts.baseName}_counts.csv.summary" into chFeatureCountsLogs

  script:
  def featureCountsDirection = 0
  if (parseRes == 'forward'){
      featureCountsDirection = 1
  } else if ((parseRes == 'reverse')){
      featureCountsDirection = 2
  }
  """
  featureCounts ${params.featurecountsOpts} -T ${task.cpus} -a $gtf -o ${bamFeaturecounts.baseName}_counts.csv -p -s $featureCountsDirection $bamFeaturecounts
  """
}

process HTseqCounts {
  tag "${bamHTseqCounts}"
  publishDir "${params.outdir}/counts", mode: 'copy',
    saveAs: {filename ->
      if (filename.indexOf("_gene.HTseqCounts.txt.summary") > 0) "gene_count_summaries/$filename"
      else if (filename.indexOf("_gene.HTseqCounts.txt") > 0) "gene_counts/$filename"
      else "$filename"
    }
  
  when:
  params.counts == 'HTseqCounts'

  input:
  file bamHTseqCounts from chBamHTseqCounts
  file gtf from chGtfHTseqCounts.collect()
  val parseRes from  chStrandedResultsHTseqCounts

  output: 
  file "*_counts.csv" into chHtseqCountsToMerge, chHtseqCountsToR, chHTSeqCountsLogs

  script:
  def strandedOpt = '-s no' 
  if (parseRes == 'forward'){
      strandedOpt= '-s yes'
  } else if ((parseRes == 'reverse')){
      strandedOpt= '-s reverse'
  }
  """
  htseq-count ${params.htseqOpts} $strandedOpt $bamHTseqCounts $gtf > ${bamHTseqCounts.baseName}_counts.csv
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
  publishDir "${params.outdir}/counts", mode: 'copy'

  input:
  file inputCounts from chCountsToMerge.collect()
  file gtf from chGtfTable.collect()
  val parseRes from chStrandedResultsTable.collect()

  output:
  file 'tablecounts_raw.csv' into chRawCounts, chCountsSaturation
  file 'tablecounts_tpm.csv' into chTpmCounts, chTpmGenetype
  file 'tableannot.csv' into chGenesAnnot

  script:
  """
  echo -e ${inputCounts} | tr " " "\n" > listofcounts.tsv
  echo -n "${parseRes}" | sed -e "s/\\[//" -e "s/\\]//" -e "s/,//g" | tr " " "\n" > listofstrandness.tsv
  makeCountTable.r listofcounts.tsv ${gtf} ${params.counts} listofstrandness.tsv
  """
}

chCountsLogs = Channel.empty()
if( params.counts == 'featureCounts' ){
  chCountsLogs = chFeatureCountsLogs
} else if (params.counts == 'HTseqCounts'){
  chCountsLogs = chHTSeqCountsLogs
}else if (params.counts == 'star'){
  chCountsLogs = chStarLogCounts
}


/*
 * Gene-based saturation
 */

process geneSaturation {
  publishDir "${params.outdir}/gene_saturation" , mode: 'copy'

  when:
  !params.skipQC && !params.skipSaturation

  input:
  file inputCounts from chCountsSaturation.collect()

  output:
  file "*gcurve.txt" into chGenesatResults

  script:
  """
  gene_saturation.r $inputCounts counts.gcurve.txt
  """
}


/*
 * Reads distribution
 */

process readDistribution {
  tag "${bamReadDist}"
  publishDir "${params.outdir}/read_distribution" , mode: 'copy'

  when:
  !params.skipReaddist

  input:
  file bamReadDist from chBamReadDist
  file bed12 from chBedReadDist.collect()

  output:
  file "*.txt" into chReadDistResults

  script:
  """
  read_distribution.py -i ${bamReadDist} -r ${bed12} > ${bamReadDist.baseName}.read_distribution.txt
  """
}


process getCountsPerGeneType {
  publishDir "${params.outdir}/read_distribution", mode: 'copy'

  when:
  !params.skipReaddist

  input:
  file tpmGenetype from chTpmGenetype
  file gtf from chGtfGenetype.collect()
 
  output:
  file "*genetype.txt" into chCountsPerGenetype

  script:
  """
  gene_type_expression.r ${tpmGenetype} ${gtf} counts_genetype.txt 
  """
}


/*
 * Exploratory analysis
 */

process exploratoryAnalysis {
  publishDir "${params.outdir}/exploratory_analysis", mode: 'copy'

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

  script:
  """
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

process getSoftwareVersions {
  output:
  file 'software_versions_mqc.yaml' into softwareVersionsYaml

  script:
  """
  echo $workflow.manifest.version &> v_rnaseq.txt
  echo $workflow.nextflow.version &> v_nextflow.txt
  fastqc --version &> v_fastqc.txt
  STAR --version &> v_star.txt
  tophat2 --version &> v_tophat2.txt
  hisat2 --version &> v_hisat2.txt
  preseq &> v_preseq.txt
  infer_experiment.py --version &> v_rseqc.txt
  read_duplication.py --version &> v_read_duplication.txt
  featureCounts -v &> v_featurecounts.txt
  htseq-count -h | grep version  &> v_htseq.txt
  picard MarkDuplicates --version &> v_markduplicates.txt || true
  samtools --version &> v_samtools.txt
  multiqc --version &> v_multiqc.txt
  scrape_software_versions.py &> software_versions_mqc.yaml
  """
}

process workflowSummaryMqc {
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

process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    when:
    !params.skipMultiQC

    input:
    file splan from chSplan.collect()
    file metadata from chMetadata.ifEmpty([])
    file multiqcConfig from chMultiqcConfig    
    file (fastqc:'fastqc/*') from chFastqcResults.collect().ifEmpty([])
    file ('rrna/*') from chRrnaLogs.collect().ifEmpty([])
    file ('alignment/*') from chAlignmentLogs.collect()
    file ('strandness/*') from chStrandnessResults.collect().ifEmpty([])
    file ('rseqc/*') from chReadDistResults.collect().ifEmpty([])
    file ('rseqc/*') from chGenebodyCoverageResults.collect().ifEmpty([])
    file ('preseq/*') from chPreseqResults.collect().ifEmpty([])
    file ('genesat/*') from chGenesatResults.collect().ifEmpty([])
    file ('dupradar/*') from chDupradarResults.collect().ifEmpty([])
    file ('picard/*') from chPicardResults.collect().ifEmpty([])	
    file ('counts/*') from chCountsLogs.collect()
    file ('genetype/*') from chCountsPerGenetype.collect().ifEmpty([])
    file ('exploratory_analysis_results/*') from chExploratoryAnalysisResults.collect().ifEmpty([]) 
    file ('software_versions/*') from softwareVersionsYaml.collect()
    file ('workflow_summary/*') from workflowSummaryYaml.collect()

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
    
    modulesList = "-m custom_content -m preseq -m rseqc -m bowtie1 -m hisat2 -m star -m tophat -m cutadapt -m fastqc"
    modulesList = params.counts == 'featureCounts' ? "${modulesList} -m featureCounts" : "${modulesList}"  
    modulesList = params.counts == 'HTseqCounts' ? "${modulesList} -m htseq" : "${modulesList}"  
 
    """
    stats2multiqc.sh ${splan} ${params.aligner} ${isPE}
    medianReadNb="\$(sort -t, -k3,3n mq.stats | awk -F, '{a[i++]=\$3;} END{x=int((i+1)/2); if (x<(i+1)/2) printf "%.0f", (a[x-1]+a[x])/2; else printf "%.0f",a[x-1];}')"
    mqc_header.py --name "RNA-seq" --version ${workflow.manifest.version} ${metadataOpts} ${splanOpts} --nbreads \${medianReadNb} > multiqc-config-header.yaml
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
  publishDir "${params.outdir}/pipeline_info", mode: 'copy'

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
    def tf = new File("$baseDir/assets/oncompleteTemplate.txt")
    def txt_template = engine.createTemplate(tf).make(report_fields)
    def report_txt = txt_template.toString()
    
    // Render the HTML template
    def hf = new File("$baseDir/assets/oncompleteTemplate.html")
    def html_template = engine.createTemplate(hf).make(report_fields)
    def report_html = html_template.toString()

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << report_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << report_txt }

    /*oncomplete file*/

    File woc = new File("${params.outdir}/workflow.oncomplete.txt")
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
