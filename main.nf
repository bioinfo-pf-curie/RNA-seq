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

// Custom functions/variables
mqcReport = []
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
params.transcriptsFasta = NFTools.getGenomeAttribute(params, 'transcriptsFasta')
params.bed12 = NFTools.getGenomeAttribute(params, 'bed12')
params.fasta = NFTools.getGenomeAttribute(params, 'fasta')
params.fastaFai = NFTools.getGenomeAttribute(params, 'fastaFai')
params.polym = NFTools.getGenomeAttribute(params, 'polym')
params.starOptions = params.genomes[ params.genome ].starOpts ? NFTools.getGenomeAttribute(params, 'starOpts') : params.starOpts
params.salmonIndex = NFTools.getGenomeAttribute(params, 'salmon')
params.gencode = NFTools.getGenomeAttribute(params, 'gencode')

// Stage config files
chMultiqcConfig = Channel.fromPath(params.multiqcConfig)
chOutputDocs = Channel.fromPath("$baseDir/docs/output.md")
chOutputDocsImages = file("$baseDir/docs/images/", checkIfExists: true)
chPcaHeader = Channel.fromPath("$baseDir/assets/pcaHeader.txt")
chHeatmapHeader = Channel.fromPath("$baseDir/assets/heatmapHeader.txt")

// Tools
starTwoPassOpts = params.starTwoPass ? '--twopassMode Basic' : ''
starCountsOpts = params.counts == 'star' ? '--quantMode GeneCounts' : ''
starBAMOpts = params.counts == 'salmon' ? '--quantMode TranscriptomeSAM' : ''
params.starAlignOptions = "${params.starOptions} ${starTwoPassOpts} ${starCountsOpts} ${starBAMOpts}"

denovoTools = params.denovo ? params.denovo.split(',').collect{it.trim().toLowerCase()} : []

/*
==========================
 VALIDATE INPUTS
==========================
*/

if (!params.genome){
  exit 1, "No genome provided. The --genome option is mandatory"
}

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
  exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

if (!params.pseudoAligner && !params.aligner){
  exit 1, "Please provide a pseudo-aligner and an aligner, using either of the '--aligner' and '--pseudoAligner' parameters."
}

if (params.pseudoAligner && params.aligner){
  exit 1, "Cannot use both a pseudo-aligner and an aligner. Please use either of the '--aligner' and '--pseudoAligner' parameters." 
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
  chSalmonIndex = Channel.empty()
}
else if ( params.hisat2Index && params.aligner == 'hisat2' ){
  Channel
    .fromPath("${params.hisat2Index}")
    .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2Index}" }
    .set{chHisat2Index}
  chStarIndex = Channel.empty()
  chSalmonIndex = Channel.empty()
}
else if ( params.salmonIndex && params.pseudoAligner == "salmon" ){
  Channel
    .fromPath("${params.salmonIndex}")
    .ifEmpty { exit 1, "Salmon index not found: ${params.salmonIndex}" }
    .set{chSalmonIndex}
  chStarIndex = Channel.empty()
  chHisat2Index = Channel.empty()
}
else {
    exit 1, "No genome index specified!"
}

if (params.bowtie2Index){
  Channel
    .fromPath("${params.bowtie2Index}")
    .ifEmpty { exit 1, "Bowtie2 index not found: ${params.bowtie2Index}" }
    .set{chBowtie2Index}
}else{
  chBowtie2Index = Channel.empty()
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

if( params.transcriptsFasta ){
  Channel
    .fromPath(params.transcriptsFasta)
    .ifEmpty { exit 1, "Transcripts fasta file not found: ${params.transcriptsFasta}" }
    .set { chTranscriptsFasta }
}else {
  chTranscriptsFasta = Channel.empty()
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
  'Gencode' : params.gencode ? 'yes' : 'no',
  'BED Annotation' : params.bed12 ?: null,
  'Identito' : params.polym ?: null,
  'Strandedness' : params.stranded,
  'Aligner' : params.aligner ?: null,
  'Star TwoPass' : params.starTwoPass ?:null,
  'PseudoAligner' : params.pseudoAligner ?: null,
  'Guided Assembly' : params.denovo ?: null,
  'Counts' : params.counts,
  'Max Resources': "${params.maxMemory} memory, ${params.maxCpus} cpus, ${params.maxTime} time per job",
  'Container': workflow.containerEngine && workflow.container ? "${workflow.containerEngine} - ${workflow.container}" : null,
  'Profile' : workflow.profile,
  'OutDir' : params.outDir,
  'WorkDir': workflow.workDir
].findAll{ it.value != null }

workflowSummaryCh = NFTools.summarize(summary, workflow, params)

/*
==============================
  LOAD INPUT DATA
==============================
*/

// Load raw reads
chRawReads = NFTools.getInputData(params.samplePlan, params.reads, params.readPaths, params.singleEnd, params)

// Make samplePlan if not available
chSplan = NFTools.getSamplePlan(params.samplePlan, params.reads, params.readPaths, params.singleEnd)

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
include { salmonQuantFromBamFlow } from './nf-modules/local/subworkflow/salmonQuantFromBam'
include { salmonQuantFromFastqFlow } from './nf-modules/local/subworkflow/salmonQuantFromFastq'
include { geneCountsAnalysisFlow } from './nf-modules/local/subworkflow/geneCountsAnalysis'
include { stringtieFlow } from './nf-modules/local/subworkflow/stringtie'
include { scallopFlow } from './nf-modules/local/subworkflow/scallop'

// Processes
include { getSoftwareVersions } from './nf-modules/local/process/getSoftwareVersions'
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

    // Init MultiQC Channels
    chFastqcMqc = Channel.empty()
    chrRNAMappingMqc = Channel.empty()
    chAlignedMqc = Channel.empty()
    chQualimapMqc = Channel.empty()
    chPreseqMqc = Channel.empty()
    chMarkDupMqc = Channel.empty()
    chDupradarMqc = Channel.empty()
    chIndentitoMqc = Channel.empty() 
    chCountsMqc = Channel.empty()
    chGeneSatResults=Channel.empty()
    chGeneTypeResults=Channel.empty()
    chGeneExpAnResults=Channel.empty()

    // Init Channels
    chCounts = Channel.empty()
    chCountsTpm = Channel.empty()

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
      chFastqMqc = fastqc.out.results.collect()
      chVersions = chVersions.mix(fastqc.out.versions)
    }

    // SUBWORKFLOW: Strandness rseqc
    strandnessFlow(
      chRawReads,
      chBedRseqc,
      chBowtie2Index
    )
    chVersions = chVersions.mix(strandnessFlow.out.versions)

    //*****************************************
    // ALIGNMENT-BASED ANALYSIS

    if (!params.pseudoAligner && params.aligner){

      // PROCESS: rRNA mapping 
      if (!params.skipRrna && params.rrna){
        rRNAMapping(
          chRawReads,
          chRrnaAnnot.collect()
        )
        chFilteredReads = rRNAMapping.out.filteredReads
	chrRNAMappingMqc = rRNAMapping.out.logs.collect()
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
        chAlignedMqc = mappingStarFlow.out.logs
        chAlignedFlagstat = mappingStarFlow.out.flagstat
        chVersions = chVersions.mix(mappingStarFlow.out.versions)
      }

      // SUBWORKFLOW: HISAT2 mapping      
      if (params.aligner == "hisat2"){
        mappingHisat2Flow(
          chFilteredReads,
	  strandnessFlow.out.strandnessResults,
          chHisat2Index,
          chGtf
        )
        chAlignedBam = mappingHisat2Flow.out.bam
        chAlignedBai = mappingHisat2Flow.out.bai
        chAlignedMqc = mappingHisat2Flow.out.logs
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
          chBamPassed.join(strandnessFlow.out.strandnessResults)
        )
        chVersions = chVersions.mix(bigWig.out.versions)
      }

      // PROCESS : Qualimap
      if (!params.skipQC && !params.skipQualimap){
        qualimap(
          chBamPassed.join(strandnessFlow.out.strandnessResults),
          chGtf.collect()
        )
	chQualimapMqc = qualimap.out.results.collect()
        chVersions = chVersions.mix(qualimap.out.versions)
      }

      // PROCESS : Saturation curves
      if (!params.skipQC && !params.skipSaturation){ 
        preseq(
          chBamPassed
        )
	chPreseqMqc = preseq.out.results.collect()
        chVersions = chVersions.mix(preseq.out.versions)
      }
      
      // SUBWORKFLOW: Duplicates
      markdupFlow(
        chBamPassed,
	strandnessFlow.out.strandnessResults,
        chGtf.collect()
      )
      chMarkDupMqc = markdupFlow.out.picardMetrics.collect()
      chDupradarMqc = markdupFlow.out.dupradarResults.collect()
      chVersions = chVersions.mix(markdupFlow.out.versions)

      // SUBWORKFLOW: Identito - polym and Monitoring
      if (!params.skipIdentito){
        identitoFlow(
          chBamPassed,
          chFasta.collect(),
          chFastaFai.collect(),
          chPolymBed.collect()
        )
	chIdentitoMqc = identitoFlow.out.results.collect()
        chVersions = chVersions.mix(identitoFlow.out.versions)
      }

      // SUBWORKFLOW: Counts
      if(params.counts == 'featureCounts'){
        featureCountsFlow(
          chBamPassed,
	  strandnessFlow.out.strandnessResults,
          chGtf.collect()
        )
        chCounts = featureCountsFlow.out.counts
        chCountsTpm = featureCountsFlow.out.tpm
        chCountsMqc = featureCountsFlow.out.logs
        chVersions = chVersions.mix(featureCountsFlow.out.versions)
      } else if (params.counts == 'HTseqCounts'){
        htseqCountsFlow (
          chBamPassed,
	  strandnessFlow.out.strandnessResults,
          chGtf.collect()
        )
        chCounts = htseqCountsFlow.out.counts
        chCountsTpm = htseqCountsFlow.out.tpm
        chCountsMqc = htseqCountsFlow.out.logs
        chVersions = chVersions.mix(htseqCountsFlow.out.versions)
      } else if (params.counts == 'star'){
        starCountsFlow (
          mappingStarFlow.out.counts,
	  mappingStarFlow.out.logs,
	  strandnessFlow.out.strandnessResults,
          chGtf.collect()
        )
        chCounts = starCountsFlow.out.counts
        chCountsTpm = starCountsFlow.out.tpm
        chCountsMqc = starCountsFlow.out.logs
        chVersions = chVersions.mix(starCountsFlow.out.versions)
      } else if (params.counts == 'salmon'){
        salmonQuantFromBamFlow (
          mappingStarFlow.out.transcriptsBam,
	  strandnessFlow.out.strandnessResults, 
	  chTranscriptsFasta,
	  chGtf
        )
        chCounts = salmonQuantFromBamFlow.out.countsGene
        chCountsTpm = salmonQuantFromBamFlow.out.tpmGene
        chCountsMqc = salmonQuantFromBamFlow.out.results
        chVersions = chVersions.mix(salmonQuantFromBamFlow.out.versions)
      }
    }


    //*****************************************
    // PSEUDO-ALIGNMENT-BASED ANALYSIS

    if (params.pseudoAligner == "salmon"){
      salmonQuantFromFastqFlow (
        chRawReads,
	strandnessFlow.out.strandnessResults,
	chSalmonIndex,
	chGtf
      )
      chCounts = salmonQuantFromFastqFlow.out.countsGene
      chCountsTpm = salmonQuantFromFastqFlow.out.tpmGene
      chCountsMqc = salmonQuantFromFastqFlow.out.results
      chVersions = chVersions.mix(salmonQuantFromFastqFlow.out.versions)
    }


    //******************************************
    // COUNTS-BASED QC
 
    // SUBWORKFLOW: gene counts qc
    if (!params.skipGeneCountsAnalysis && params.counts && (params.aligner || params.pseudoAligner)){
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
    }

    //*******************************************
    // GUIDED DE NOVO ASSEMBLY
    chgffCompareMqc = Channel.empty()
    if ("stringtie" in denovoTools){
      stringtieFlow(
        chBamPassed,
	strandnessFlow.out.strandnessResults,
        chGtf.collect()
      )
      chVersions = chVersions.mix(stringtieFlow.out.versions)
      chgffCompareMqc = chgffCompareMqc.mix(stringtieFlow.out.mqc)
    }

    if ("scallop" in denovoTools){
      scallopFlow(
        chBamPassed,
	strandnessFlow.out.strandnessResults,
	chGtf.collect()
      )
      chVersions = chVersions.mix(scallopFlow.out.versions)
      chgffCompareMqc = chgffCompareMqc.mix(scallopFlow.out.mqc)
    }

    //*******************************************
    // MULTIQC

    if (!params.skipMultiQC){

      getSoftwareVersions(
        chVersions.unique().collectFile()
      )

      // Warnings
      chWarnMapping = Channel.empty()
      if (skippedPoorAlignment.size() > 0){
        Channel.fromList(skippedPoorAlignment)
             .flatMap{ it -> it + ": Poor alignment rate. Sample discarded"}
             .set{chWarnMapping}
      }
      
      strandnessFlow.out.strandnessResults
        .map{it[1]}
        .unique()
	.count()
	.filter{it > 1}
	.flatMap{"Samples with different strandness detected !"}
	.set{chWarnStrand}

     chWarnStrand
       .mix(chWarnMapping)
       .collectFile(name: 'warnings.txt', newLine: true)
       .set{chWarn}

    
      multiqc(
        customRunName,
        chSplan.collect(),
        chMetadata.ifEmpty([]),
        chMultiqcConfig.ifEmpty([]),
        chFastqMqc.ifEmpty([]),
        chrRNAMappingMqc.ifEmpty([]),
        chAlignedMqc.collect().ifEmpty([]),
        strandnessFlow.out.strandnessOutputFiles.collect().ifEmpty([]),
        chQualimapMqc.ifEmpty([]),
        chPreseqMqc.ifEmpty([]),
        chIdentitoMqc.ifEmpty([]),
        chMarkDupMqc.ifEmpty([]),
        chDupradarMqc.ifEmpty([]),
        chCountsMqc.collect().ifEmpty([]),
        chGeneSatResults.collect().ifEmpty([]),
        chGeneTypeResults.collect().ifEmpty([]),
        chGeneExpAnResults.collect().ifEmpty([]),
	chgffCompareMqc.collect().ifEmpty([]),
        getSoftwareVersions.out.versionsYaml.collect().ifEmpty([]),
        workflowSummaryCh.collectFile(name: "workflow_summary_mqc.yaml"),
        chWarn.collect().ifEmpty([])
      )
      mqcReport = multiqc.out.report.toList()
    }
}

workflow.onComplete {
  NFTools.makeReports(workflow, params, summary, customRunName, mqcReport)
}
