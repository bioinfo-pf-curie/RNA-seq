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

include { helpMessage } from './functions'
// Show help emssage
if (params.help){
  helpMessage()
  exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

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

skippedPoorAlignment = [] 

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

chHisat2Index = Channel.empty()
chStarIndex = Channel.empty()
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
    .set { chGtf }
}else {
  log.warn "No GTF annotation specified - dupRadar, table counts - will be skipped !" 
  Channel
    .empty()
    .set { chGtf } 
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

// Workflows
include { qcFlow } from './nf-modules/local/subworkflow/qc'
include { rseqFlow } from './nf-modules/local/subworkflow/rseq'
include { mappingFlow } from './nf-modules/local/subworkflow/mapping'
include { markdupFlow } from './nf-modules/local/subworkflow/markdup'
include { polymFlow } from './nf-modules/local/subworkflow/polym'
include { countsFlow } from './nf-modules/local/subworkflow/counts'

// Processes
include { getSoftwareVersions } from './nf-modules/local/process/getSoftwareVersions'
include { workflowSummaryMqc } from './nf-modules/local/process/workflowSummaryMqc'
include { multiqc } from './nf-modules/local/process/multiqc'
include { outputDocumentation } from './nf-modules/local/process/outputDocumentation'
include { bigWig } from './nf-modules/local/process/bigWig'
include { qualimap } from './nf-modules/local/process/qualimap'
include { preseq } from './nf-modules/local/process/preseq'

// WORKFLOW principal
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
        rseqFlow.out.chStrandnessResults,
        skippedPoorAlignment
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
      countsFlow(
        mappingFlow.out.chBam,
        chGtf,
        rseqFlow.out.chStrandedResults,
        mappingFlow.out.chStarCounts,
        mappingFlow.out.chStarLogCounts,
        chPcaHeader,
        chHeatmapHeader
      )

      // MultiQC

      getSoftwareVersions(
        qcFlow.out.chFastqcVersion.first().ifEmpty([]),
        mappingFlow.out.chStarVersion.first().ifEmpty([]),
        mappingFlow.out.chHisat2Version.first().ifEmpty([]),
        mappingFlow.out.chBowtieVersion.first().ifEmpty([]),
        rseqFlow.out.chBowtie2Version.first().ifEmpty([]),
        mappingFlow.out.chSamtoolsVersionSort.first().ifEmpty([]),
        markdupFlow.out.chPicardVersion.first().ifEmpty([]),
        preseq.out.chPreseqVersion.first().ifEmpty([]),
        countsFlow.out.chMergeCountsVersion.concat(
          polymFlow.out.chCombinePolymVersion,
          countsFlow.out.chGeneSaturationVersion,
          countsFlow.out.chGeneTypeVersion,
          countsFlow.out.chAnaExpVersion).first().ifEmpty([]),
        rseqFlow.out.chRseqcVersionInferExperiment.first().ifEmpty([]),
        countsFlow.out.chFeaturecountsVersion.first().ifEmpty([]),
        bigWig.out.chDeeptoolsVersion.first().ifEmpty([]),
        polymFlow.out.chBcftoolsVersion.first().ifEmpty([]),
        countsFlow.out.chHtseqVersion.first().ifEmpty([]),
        qualimap.out.chQualimapVersion.first().ifEmpty([])
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
        customRunName,
        chSplan.collect(),
        chMetadata.ifEmpty([]),
        chMultiqcConfig.ifEmpty([]),
        qcFlow.out.chFastqcResults.collect().ifEmpty([]),
        mappingFlow.out.chRrnaLogs.collect().ifEmpty([]),
        mappingFlow.out.chAlignmentLogs.collect().ifEmpty([]),
        rseqFlow.out.chStrandnessResults.collect().ifEmpty([]),
        qualimap.out.chQualimapResults.collect().ifEmpty([]),
        preseq.out.chPreseqResults.collect().ifEmpty([]),
        countsFlow.out.chGenesatResults.collect().ifEmpty([]),
        markdupFlow.out.chDupradarResults.collect().ifEmpty([]),
        markdupFlow.out.chPicardResults.collect().ifEmpty([]),
        countsFlow.out.chCountsLogs.collect().ifEmpty([]),
        countsFlow.out.chCountsPerGenetype.collect().ifEmpty([]),
        polymFlow.out.chClustPolymResults.collect().ifEmpty([]),
        countsFlow.out.chExploratoryAnalysisResults.collect().ifEmpty([]),
        getSoftwareVersions.out.chSoftwareVersionsYaml.collect().ifEmpty([]),
        workflowSummaryMqc.out.chWorkflowSummaryYaml.collect().ifEmpty([]),
        chWarn.collect().ifEmpty([]),
        skippedPoorAlignment
      )
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
