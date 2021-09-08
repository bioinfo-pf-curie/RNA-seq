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
                         RNA-seq DSL2
========================================================================================
 RNA-seq Analysis Pipeline.
 #### Homepage / Documentation
 https://gitlab.curie.fr/data-analysis/rnaseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

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
      .set { chRawReads }
  }else{
     Channel
       .from(file("${params.samplePlan}"))
       .splitCsv(header: false)
       .map{ row -> [ row[0], [file(row[2]), file(row[3])]] }

      .set { chRawReads } 
   }
   params.reads=false
}
else if(params.readPaths){
  if(params.singleEnd){
    Channel
      .from(params.readPaths)
      .map { row -> [ row[0], [file(row[1][0])]] }
      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
      .set { chRawReads }
    } else {
    Channel
      .from(params.readPaths)
      .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
      .set { chRawReads } 
  }
} else {
  Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
.set { chRawReads }
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
  label 'minCpu'
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
  label 'minCpu'
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
  label 'minCpu'
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
  label 'minCpu'
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

                                                                                                                                                                                                  
// Workflows
include { qcFlow } from './nf-modules/local/subworkflow/qc'
include { rseqFlow } from './nf-modules/local/subworkflow/rseq'
include { mappingFlow } from './nf-modules/local/subworkflow/mapping'
include { markdupFlow } from './nf-modules/local/subworkflow/markdup'
include { polymFlow } from './nf-modules/local/subworkflow/polym'

// Processes
include { getSoftwareVersions } from './nf-modules/local/process/getSoftwareVersions'
include { workflowSummaryMqc } from './nf-modules/local/process/workflowSummaryMqc'
include { multiqc } from './nf-modules/local/process/multiqc'
include { outputDocumentation } from './nf-modules/local/process/outputDocumentation'
include { bigWig } from './nf-modules/local/process/bigWig'
include { qualimap } from './nf-modules/local/process/qualimap'
include { preseq } from './nf-modules/local/process/preseq'

workflow {
    main:

      // subroutines
      outputDocumentation(
        chOutputDocs,
        chOutputDocsImages
      )

      // SUBWORKFLOW: QC : check factqc
      qcFlow(
        chRawReads
      )

      // SUBWORKFLOW: Strandness Rseq
      rseqFlow(
        chRawReads,
        chBedRseqc
      )


      // SUBWORKFLOW: mapping (rRNA and Read mapping)
      mappingFlow(
        chRawReads,
        chRrnaAnnot,
        chStarIndex,
        chGtf,
        chHisat2Index,
        rseqFlow.out.chStrandnessResults
      )

      // Generate bigwig file
      bigWig(
        mappingFlow.out.chBam,
        rseqFlow.out.chStrandedResults
      )

      // Qualimap
      qualimap(
        mappingFlow.out.chBam,
        chGtf.collect(),
        rseqFlow.out.chStrandedResults
      )

      // Saturation Curves
      preseq(
        mappingFlow.out.chBam
      )
      
      // SUBWORKFLOW: Duplicates
      markdupFlow(
        mappingFlow.out.chBam,
        chGtf,
        rseqFlow.out.chStrandedResults
      )

      // SUBWORKFLOW: Identito - polym
      polymFlow(
        chFasta,
        chPolymBed,
        markdupFlow.out.chBamMd
      )
      // SUBWORKFLOW: Counts

      // MultiQC

      getSoftwareVersions(
        qcFlow.out.chFastqcVersion.first().ifEmpty([]),
        mappingFlow.out.chStarVersion.first().ifEmpty([]),
        mappingFlow.out.chHisat2Version.first().ifEmpty([]),
        mappingFlow.out.chBowtieVersion.first().ifEmpty([]),
        rseqFlow.out.chBowtie2Version.first().ifEmpty([]),
        mappingFlow.out.chSamtoolsVersionSort.first().ifEmpty([]),
        chPicardVersion.first().ifEmpty([]),
        preseq.out.chPreseqVersion.first().ifEmpty([]),
        chMergeCountsVersion.concat(polymFlow.out.chCombinePolymVersion,chGeneSaturationVersion,chGeneTypeVersion,chAnaExpVersion).first().ifEmpty([]),
        rseqFlow.out.chRseqcVersionInferExperiment.first().ifEmpty([]),
        chFeaturecountsVersion.first().ifEmpty([]),
        chDeeptoolsVersion.first().ifEmpty([]),
        polymFlow.out.chBcftoolsVersion.first().ifEmpty([]),
        chHtseqVersion.first().ifEmpty([]),
        qualimap.out.chQualimapVersion.first().ifEmpty([]
      )

      workflowSummaryMqc(
        summary
      )

      if (skippedPoorAlignment.size() > 0){
        Channel.fromList(skippedPoorAlignment)
               .flatMap{ it -> it + ": Poor alignment rate. Sample discarded"}
               .collectFile(name: 'warnings.txt', newLine: true)
               .set{chWarn}
      }else{
         chWarn = Channel.empty()
      }

      multiqc(
        custom_runName,
        chSplan.collect(),
        chMetadata.ifEmpty([]),
        chMultiqcConfig.ifEmpty([]),
        qcFlow.out.chFastqcResults.collect().ifEmpty([]),
        chRrnaLogs.collect().ifEmpty([]),
        mappingFlow.out.chAlignmentLogs.collect().ifEmpty([]),
        rseqFlow.out.chStrandnessResults.collect().ifEmpty([]),
        qualimap.out.chQualimapResults.collect().ifEmpty([]),
        preseq.out.chPreseqResults.collect().ifEmpty([]),
        chGenesatResults.collect().ifEmpty([]),
        chDupradarResults.collect().ifEmpty([]),
        chPicardResults.collect().ifEmpty([]),
        chCountsLogs.collect().ifEmpty([]),
        chCountsPerGenetype.collect().ifEmpty([]),
        polymFlow.out.chClustPolymResults.collect().ifEmpty([]),
        chExploratoryAnalysisResults.collect().ifEmpty([]),
        chSoftwareVersionsYaml.collect().ifEmpty([]),
        chWorkflowSummaryYaml.collect().ifEmpty([]),
        chWarn.collect().ifEmpty([]) 
      }
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
