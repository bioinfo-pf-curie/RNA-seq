#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2019-2020
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data.
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.
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

include { helpMessage; checkAlignmentPercent } from './lib/functions'
// Show help message
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
params.fastaFai = params.genome ? params.genomes[ params.genome ].fastaFai ?: false : false
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

if( params.starIndex && params.aligner == 'star' ){
  Channel
    .fromPath(params.starIndex)
    .ifEmpty { exit 1, "STAR index not found: ${params.starIndex}" }
    .set {chStarIndex}
  chHisat2Index = Channel.empty()
}
else if ( params.hisat2Index && params.aligner == 'hisat2' ){
  Channel
    .fromPath("${params.hisat2Index}*")
    .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2Index}" }
    .set{chHisat2Index}
  chStarIndex = Channel.empty()
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

if ( params.fastaFai ){
  Channel
    .fromPath( params.fastaFai )
    .ifEmpty { exit 1, "Genome fastaFai file not found: ${params.fastaFai}" }
    .set { chFastaFai }
}else {
  chFastaFai = Channel.empty()
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
 * INCLUDE
 */ 

// Workflows
include { strandnessFlow } from './nf-modules/local/subworkflow/strandness'
include { mappingStarFlow } from './nf-modules/local/subworkflow/mappingStar'
include { mappingHisat2Flow } from './nf-modules/local/subworkflow/mappingHisat2'
include { markdupFlow } from './nf-modules/local/subworkflow/markdup'
include { identitoFlow } from './nf-modules/local/subworkflow/identito'
include { featureCountsFlow } from './nf-modules/local/subworkflow/featureCounts'
include { htseqCountsFlow } from './nf-modules/local/subworkflow/htseqCounts'
include { starCountsFlow } from './nf-modules/local/subworkflow/starCounts'
include { geneCountsAnalysisFlow } from './nf-modules/local/subworkflow/geneCountsAnalysis'

// Processes
include { getSoftwareVersions } from './nf-modules/local/process/getSoftwareVersions'
include { workflowSummaryMqc } from './nf-modules/local/process/workflowSummaryMqc'
include { multiqc } from './nf-modules/local/process/multiqc'
include { outputDocumentation } from './nf-modules/local/process/outputDocumentation'
include { bigWig } from './nf-modules/local/process/bigWig'
include { qualimap } from './nf-modules/local/process/qualimap'
include { preseq } from './nf-modules/local/process/preseq'
include { fastqc } from './nf-modules/local/process/fastqc'
include { rRNAMapping } from './nf-modules/local/process/rRNAMapping'

/*
 * WORKFLOW 
 */

workflow {
  chVersions = Channel.empty()

  main:

    // subroutines
    outputDocumentation(
      chOutputDocs,
      chOutputDocsImages
    )

    // PROCESS: fastqc
    if (!params.skipQC && !params.skipFastqc){
      fastqc(
        chRawReads
      )
      chVersions = chVersions.mix(fastqc.out.versions)
    }

    // SUBWORKFLOW: Strandness rseqc
    strandnessFlow(
      chRawReads,
      chBedRseqc
    )
    chVersions = chVersions.mix(strandnessFlow.out.versions)

    // PROCESS: rRNA mapping 
    if (!params.skipRrna && params.rrna){
      rRNAMapping(
        chRawReads,
        chRrnaAnnot.collect()
      )
      chFilteredReads = rRNAMapping.out.filteredReads
      chVersions = chVersions.mix(rRNAMapping.out.versions)
    }else{
      chFilteredReads = chRawReads
    }

    // SUBWORKFLOW: STAR mapping
    if (params.aligner == "star"){
      mappingStarFlow(
        chFilteredReads,
        chStarIndex,
        chGtf
      )
      chAlignedBam = mappingStarFlow.out.bam
      chAlignedBai = mappingStarFlow.out.bai
      chAlignedLogs = mappingStarFlow.out.logs
      chAlignedFlagstat = mappingStarFlow.out.flagstat
      chVersions = chVersions.mix(mappingStarFlow.out.versions)
    }

    // SUBWORKFLOW: HISAT2 mapping      
    if (params.aligner == "hisat2"){
      mappingHisat2Flow(
        chFilteredReads,
        chHisat2Index,
        chGtf,
        strandnessFlow.out.strandnessResults
      )
      chAlignedBam = mappingHisat2Flow.out.bam
      chAlignedBai = mappingHisat2Flow.out.bai
      chAlignedLogs = mappingHisat2Flow.out.logs
      chAlignedFlagstat = mappingHisat2Flow.out.flagstat
      chVersions = chVersions.mix(mappingHisat2Flow.out.versions)
    }

    // Filter removes all 'aligned' channels that fail the check
    chAlignedFlagstat.join(chAlignedBam).join(chAlignedBai)
      .filter { prefix, logs, bam, bai -> checkAlignmentPercent(prefix, logs) }
      .map { prefix, logs, bam, bai -> [ prefix, bam, bai ] }
      .set { chBamPassed }
  
    // PROCESS : bigwig file
    if (!params.skipBigwig){
      bigWig(
        chBamPassed,
        strandnessFlow.out.strandnessResults
      )
      chVersions = chVersions.mix(bigWig.out.versions)
    }

    // PROCESS : Qualimap
    if (!params.skipQC && !params.skipQualimap){
      qualimap(
        chBamPassed,
        chGtf.collect(),
        strandnessFlow.out.strandnessResults
      )
      chVersions = chVersions.mix(qualimap.out.versions)
    }

    // PROCESS : Saturation curves
    if (!params.skipQC && !params.skipSaturation){ 
      preseq(
        chBamPassed
      )
      chVersions = chVersions.mix(preseq.out.versions)
    }
      
    // SUBWORKFLOW: Duplicates
    markdupFlow(
        chBamPassed,
        chGtf,
        strandnessFlow.out.strandnessResults
    )
    chVersions = chVersions.mix(markdupFlow.out.versions)

    // SUBWORKFLOW: Identito - polym and Monitoring
    if (!params.skipIdentito){
      identitoFlow(
          markdupFlow.out.bam,
          chFasta,
          chFastaFai,
          chPolymBed
      )
      chVersions = chVersions.mix(identitoFlow.out.versions)
    }

    // SUBWORKFLOW: Counts
    if(params.counts == 'featureCounts'){
      featureCountsFlow(
        chBamPassed,
        chGtf.collect(),
        strandnessFlow.out.strandnessResults
      )
      chCounts = featureCountsFlow.out.counts
      chCountsTpm = featureCountsFlow.out.tpm
      chCountsLogs = featureCountsFlow.out.logs
      chVersions = chVersions.mix(featureCountsFlow.out.versions)
    } else if (params.counts == 'HTseqCounts'){
      htseqCountsFlow (
        chBamPassed,
        chGtf.collect(),
        strandnessFlow.out.strandnessResults
      )
      chCounts = htseqCountsFlow.out.counts
      chCountsTpm = htseqCountsFlow.out.tpm
      chCountsLogs = htseqCountsFlow.out.logs
      chVersions = chVersions.mix(htseqCountsFlow.out.versions)
     } else if (params.counts == 'star'){
       starCountsFlow (
         mappingStarFlow.out.counts,
         mappingStarFlow.out.logs,
         chGtf.collect(),
         strandnessFlow.out.strandnessResults
       )
      chCounts = starCountsFlow.out.counts
      chCountsTpm = starCountsFlow.out.tpm
      chCountsLogs = starCountsFlow.out.logs
      chVersions = chVersions.mix(starCountsFlow.out.versions)
     }

    // SUBWORKFLOW: gene counts qc
    geneCountsAnalysisFlow(
      chCounts,
      chCountsTpm,
      chGtf,
      chPcaHeader,
      chHeatmapHeader
    )
    chVersions = chVersions.mix(geneCountsAnalysisFlow.out.versions)

    // MultiQC
    if (!params.skipMultiQC){

      if (!params.skipSoftVersions){
        getSoftwareVersions(
          chVersions.unique().collectFile()
        )
      }

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
      fastqc.out.results.collect().ifEmpty([]),
      rRNAMapping.out.logs.collect().ifEmpty([]),
      chAlignedLogs.collect().ifEmpty([]),
      strandnessFlow.out.strandnessOutputFiles.collect().ifEmpty([]),
      qualimap.out.results.collect().ifEmpty([]),
      preseq.out.results.collect().ifEmpty([]),
      identitoFlow.out.results.collect().ifEmpty([]),
      markdupFlow.out.picardMetrics.collect().ifEmpty([]),
      markdupFlow.out.dupradarResults.collect().ifEmpty([]),
      chCountsLogs.collect().ifEmpty([]),
      geneCountsAnalysisFlow.out.geneSaturationResults.collect().ifEmpty([]),
      geneCountsAnalysisFlow.out.countsPerGenetype.collect().ifEmpty([]),
      geneCountsAnalysisFlow.out.expAnalysisResults.collect().ifEmpty([]),
      getSoftwareVersions.out.versionsYaml.collect().ifEmpty([]),
      workflowSummaryMqc.out.chWorkflowSummaryYaml.collect().ifEmpty([]),
      chWarn.collect().ifEmpty([])
      //skippedPoorAlignment
    )
  }

  //workflow.onComplete {}
}