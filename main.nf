#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2019-2021
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.
*/


/*
========================================================================================
                         RNA-seq DSL2
========================================================================================
 RNA-seq Analysis Pipeline.
  https://gitlab.curie.fr/data-analysis/rnaseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

// Initialize lintedParams and paramsWithUsage
NFTools.welcome(workflow, params)

// Use lintedParams as default params object
paramsWithUsage = NFTools.readParamsFromJsonSettings("${projectDir}/parameters.settings.json")
params.putAll(NFTools.lint(params, paramsWithUsage))

// Run name
customRunName = NFTools.checkRunName(workflow.runName, params.name)

// Custom functions
skippedPoorAlignment = []
include {checkAlignmentPercent } from './lib/functions'

/*
===================================
  SET UP CONFIGURATION VARIABLES
===================================
*/

// Genome-based variables
params.bowtie2Index = NFTools.getGenomeAttribute(params, 'bowtie2')
params.starIndex = NFTools.getGenomeAttribute(params, 'star')
params.hisat2Index = NFTools.getGenomeAttribute(params, 'hisat2')
params.rrna = NFTools.getGenomeAttribute(params, 'rrna')
params.gtf = NFTools.getGenomeAttribute(params, 'gtf')
params.bed12 = NFTools.getGenomeAttribute(params, 'bed12')
params.fasta = NFTools.getGenomeAttribute(params, 'fasta')
params.fastaFai = NFTools.getGenomeAttribute(params, 'fastafai')
params.polym = NFTools.getGenomeAttribute(params, 'polym')
params.starOptions = NFTools.getGenomeAttribute(params, 'starOpts')

// Stage config files
chMultiqcConfig = Channel.fromPath(params.multiqcConfig)
chOutputDocs = Channel.fromPath("$baseDir/docs/output.md")
chOutputDocsImages = file("$baseDir/docs/images/", checkIfExists: true)
chPcaHeader = Channel.fromPath("$baseDir/assets/pcaHeader.txt")
chHeatmapHeader = Channel.fromPath("$baseDir/assets/heatmapHeader.txt")

/*
===========================
   SUMMARY
===========================
*/

summary = [
  'Pipeline Release': workflow.revision ?: null,
  'Run Name': customRunName,
  'Inputs' : params.samplePlan ?: params.reads ?: null,
  'Genome' : params.genome,
  'GTF Annotation' : params.gtf ?: null,
  'BED Annotation' : params.bed12 ?: null,
  'Identito' : params.polym ?: null,
  'Strandedness' : params.stranded,
  'Aligner' : params.aligner,
  'Counts' : params.counts,
  'Max Resources': "${params.maxMemory} memory, ${params.maxCpus} cpus, ${params.maxTime} time per job",
  'Container': workflow.containerEngine && workflow.container ? "${workflow.containerEngine} - ${workflow.container}" : null,
  'Profile' : workflow.profile,
  'OutDir' : params.outDir,
  'WorkDir': workflow.workDir
].findAll{ it.value != null }

workflowSummaryCh = NFTools.summarize(summary, workflow, params)

/*
==========================
 VALIDATE INPUTS
==========================
*/

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
  exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

if (params.counts == 'star' && params.aligner != 'star'){
  exit 1, "Cannot run STAR counts without STAR aligner. Please check the '--aligner' and '--counts' parameters."
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
    exit 1, "No genome index specified!"
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
  log.warn "No polymorphisms available - identito monitoring will be skipped !"
  chPolymBed = Channel.empty()
}

/*
==============================
  LOAD INPUT DATA
==============================
*/

chRawReads = NFTools.getInputData(params.samplePlan, params.reads, params.readPaths, params.singleEnd, params)

//if(params.samplePlan){
//  if(params.singleEnd){
//    Channel
//      .from(file("${params.samplePlan}"))
//      .splitCsv(header: false)
//      .map{ row -> [ row[0], [file(row[2])]] }
//      .set { chRawReads }
//  }else{
//     Channel
//       .from(file("${params.samplePlan}"))
//       .splitCsv(header: false)
//       .map{ row -> [ row[0], [file(row[2]), file(row[3])]] }
//      .set { chRawReads } 
//   }
//   params.reads=false
//}
//else if(params.readPaths){
//  if(params.singleEnd){
//    Channel
//      .from(params.readPaths)
//      .map { row -> [ row[0], [file(row[1][0])]] }
//      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
//      .set { chRawReads }
//    } else {
//    Channel
//      .from(params.readPaths)
//      .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
//      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
//      .set { chRawReads } 
//  }
//} else {
//  Channel
//    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
//    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify //--singleEnd on the command line." }
//.set { chRawReads }
//}

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


/*
==================================
           INCLUDE
==================================
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
=====================================
            WORKFLOW 
=====================================
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
          chFasta.collect(),
          chFastaFai.collect(),
          chPolymBed.collect()
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
    if (!params.skipGeneCountsAnalysis){
      geneCountsAnalysisFlow(
        chCounts,
        chCountsTpm,
        chGtf,
        chPcaHeader,
        chHeatmapHeader
      )
      chGeneSatResults=geneCountsAnalysisFlow.out.geneSaturationResults
      chGeneTypeResults=geneCountsAnalysisFlow.out.countsPerGenetype
      chGeneExpAnResults=geneCountsAnalysisFlow.out.expAnalysisResults
      chVersions = chVersions.mix(geneCountsAnalysisFlow.out.versions)
    }else{
      chGeneSatResults=Channel.empty()
      chGeneTypeResults=Channel.empty()
      chGeneExpAnResults=Channel.empty()
    }

    // MultiQC
    if (!params.skipMultiQC){

      getSoftwareVersions(
        chVersions.unique().collectFile()
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
        chGeneSatResults.collect().ifEmpty([]),
        chGeneTypeResults.collect().ifEmpty([]),
        chGeneExpAnResults.collect().ifEmpty([]),
        getSoftwareVersions.out.versionsYaml.collect().ifEmpty([]),
        workflowSummaryCh.collectFile(name: "workflow_summary_mqc.yaml"),
        chWarn.collect().ifEmpty([])
      )
      mqcReport = multiqc.out.report.toList()
    }else{
      mqcReport = []
    }
}

workflow.onComplete {
  NFTools.makeReports(workflow, params, summary, customRunName, mqcReport)
}
