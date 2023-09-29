#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2019-2022
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
include {checkAlignmentPercent} from './lib/functions'
include {combineStrandness} from './lib/functions'

/*
===================================
  SET UP CONFIGURATION VARIABLES
===================================
*/

// Genome-based variables
if (!params.genome){
  exit 1, "No genome provided. The --genome option is mandatory"
}

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
  exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

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
params.salmonIndex = NFTools.getGenomeAttribute(params, 'salmon')
params.gencode = NFTools.getGenomeAttribute(params, 'gencode')
params.pdxIndex = NFTools.getGenomeAttribute(params, 'xengsort', genome='pdx')

// Stage config files
chMultiqcConfig = Channel.fromPath(params.multiqcConfig)
chOutputDocs = Channel.fromPath("$projectDir/docs/output.md")
chOutputDocsImages = file("$projectDir/docs/images/", checkIfExists: true)
chPcaHeader = Channel.fromPath("$projectDir/assets/pcaHeader.txt")
chHeatmapHeader = Channel.fromPath("$projectDir/assets/heatmapHeader.txt")

// Tools
denovoTools = params.denovo ? params.denovo.split(',').collect{it.trim().toLowerCase()} : []

/*
==========================
 VALIDATE INPUTS
==========================
*/

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

if (!params.rrna){
  log.warn "No rRNA fasta file available - rRNA mapping - will be skipped !"
}

if (params.pdx && (params.genome != "hg19" && params.genome != "hg38" && params.genome != "mm9" && params.genome != "mm10" && params.genome != "mm39" )){
  exit 1, "Unexpected reference genome for PDX analysis. Please, specify a Human or Mouse reference genome ('--genome') for PDX analysis."
}

/*
==========================
 BUILD CHANNELS
==========================
*/

chStarIndex          = params.starIndex             ? Channel.fromPath(params.starIndex, checkIfExists: true).collect()              : Channel.empty()
chHisat2Index        = params.hisat2Index           ? Channel.fromPath(params.hisat2Index, checkIfExists: true).collect()            : Channel.empty()
chSalmonIndex        = params.salmonIndex           ? Channel.fromPath(params.salmonIndex, checkIfExists: true).collect()            : Channel.empty()
chBowtie2Index       = params.bowtie2Index          ? Channel.fromPath(params.bowtie2Index, checkIfExists: true).collect()           : Channel.empty()
chPdxIndex           = params.pdxIndex              ? Channel.fromPath(params.pdxIndex, checkIfExists: true).collect()               : Channel.empty()
chFasta              = params.fasta                 ? Channel.fromPath(params.fasta, checkIfExists: true).collect()                  : Channel.empty()
chFastaFai           = params.fastaFai              ? Channel.fromPath(params.fastaFai, checkIfExists: true).collect()               : Channel.empty()
chGtf                = params.gtf                   ? Channel.fromPath(params.gtf, checkIfExists: true).collect()                    : Channel.empty()
chTranscriptsFasta   = params.transcriptsFasta      ? Channel.fromPath(params.transcriptsFasta, checkIfExists: true).collect()       : Channel.empty()
chBedRseqc           = params.bed12                 ? Channel.fromPath(params.bed12, checkIfExists: true).collect()                  : Channel.empty()
chRrnaAnnot          = params.rrna                  ? Channel.fromPath(params.rrna, checkIfExists: true).collect()                   : Channel.empty()
chPolymBed           = params.polym                 ? Channel.fromPath(params.polym, checkIfExists: true).collect()                  : Channel.empty()
chMetadata           = params.metadata              ? Channel.fromPath(params.metadata, checkIfExists: true).collect()               : Channel.empty()

/*
===========================
   SUMMARY
===========================
*/

summary = [
  'Pipeline' : workflow.manifest.name ?: null,
  'Version': workflow.manifest.version ?: null,
  'DOI': workflow.manifest.doi ?: null,
  'Run Name': customRunName,
  'Inputs' : params.samplePlan ?: params.reads ?: null,
  'PDX' : params.pdx ?: null,
  'Genome' : params.genome,
  'GTF Annotation' : params.gtf ?: null,
  'BED Annotation' : params.bed12 ?: null,
  'Strandedness' : params.stranded,
  'Aligner' : params.aligner ?: null,
  'PseudoAligner' : params.pseudoAligner ?: null,
  'Guided Assembly' : params.denovo ?: null,
  'Counts' : params.counts,
  'Max Resources': "${params.maxMemory} memory, ${params.maxCpus} cpus, ${params.maxTime} time per job",
  'Profile' : workflow.profile,
  'OutDir' : params.outDir,
  'WorkDir': workflow.workDir,
  'CommandLine': workflow.commandLine
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
include { identitoFlow } from './nf-modules/common/subworkflow/identito'

include { strandnessFlow } from './nf-modules/local/subworkflow/strandness'
include { mappingStarFlow } from './nf-modules/local/subworkflow/mappingStar'
include { mappingHisat2Flow } from './nf-modules/local/subworkflow/mappingHisat2'
include { markdupFlow } from './nf-modules/local/subworkflow/markdup'
include { featureCountsFlow } from './nf-modules/local/subworkflow/featureCounts'
include { htseqCountsFlow } from './nf-modules/local/subworkflow/htseqCounts'
include { starCountsFlow } from './nf-modules/local/subworkflow/starCounts'
include { salmonQuantFromBamFlow } from './nf-modules/local/subworkflow/salmonQuantFromBam'
include { salmonQuantFromFastqFlow } from './nf-modules/local/subworkflow/salmonQuantFromFastq'
include { geneCountsAnalysisFlow } from './nf-modules/local/subworkflow/geneCountsAnalysis'
include { stringtieFlow } from './nf-modules/local/subworkflow/stringtie'
include { scallopFlow } from './nf-modules/local/subworkflow/scallop'

// Processes
include { getSoftwareVersions } from './nf-modules/common/process/utils/getSoftwareVersions'
include { outputDocumentation } from './nf-modules/common/process/utils/outputDocumentation'
include { trimGalore } from './nf-modules/common/process/trimGalore/trimGalore'
include { xengsort } from './nf-modules/common/process/xengsort/xengsort'
include { deeptoolsBamCoverage } from './nf-modules/common/process/deeptools/deeptoolsBamCoverage'
include { qualimapRNAseq } from './nf-modules/common/process/qualimap/qualimapRNAseq'
include { preseq } from './nf-modules/common/process/preseq/preseq'
include { fastqc } from './nf-modules/common/process/fastqc/fastqc'

include { rRNAMapping } from './nf-modules/local/process/rRNAMapping'
include { multiqc } from './nf-modules/local/process/multiqc'

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
    chIdentitoMqc = Channel.empty() 
    chCountsMqc = Channel.empty()
    chGeneSatResults = Channel.empty()
    chGeneTypeResults = Channel.empty()
    chGeneExpAnResults = Channel.empty()
    chXengsortMqc = Channel.empty()
    chTrimmingMqc = Channel.empty()

    // Init Channels
    chCounts = Channel.empty()
    chCountsTpm = Channel.empty()

    // subroutines
    outputDocumentation(
      chOutputDocs,
      chOutputDocsImages
    )

    // PROCESS: trimming
    if (params.trimming){
      trimGalore(
        chRawReads
      )
      chRawReads=trimGalore.out.fastq
      chTrimmingMqc=trimGalore.out.logs.map{it->it[1]}
      chVersions = chVersions.mix(trimGalore.out.versions)
    }

    // PROCESS: fastqc
    fastqc(
      chRawReads
    )
    chFastqcMqc = fastqc.out.results.collect()
    chVersions = chVersions.mix(fastqc.out.versions)

    // SUBWORKFLOW: Strandness rseqc
    strandnessFlow(
      chRawReads,
      chBedRseqc,
      chBowtie2Index
    )
    chVersions = chVersions.mix(strandnessFlow.out.versions)

    // Combine reads and strandness information
    chRawReads = combineStrandness(chRawReads, strandnessFlow.out.strand)

    //***************************************
    // PDX

    // PROCESS: xengsort
    if (params.pdx){
      xengsort(
        chRawReads,
        chPdxIndex
      )
      chRawReads = params.genome == "hg19" || params.genome == "hg38" ? xengsort.out.fastqHuman : xengsort.out.fastqMouse
      chXengsortMqc = xengsort.out.logs
      chVersions = chVersions.mix(xengsort.out.versions)
    }

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
        .filter { meta, logs, bam, bai -> checkAlignmentPercent(meta, logs) }
        .map { meta, logs, bam, bai -> [ meta, bam, bai ] }
        .set { chBamPassed }

      // PROCESS : bigwig file
      if (!params.skipBigwig){
        deeptoolsBamCoverage(
          chBamPassed.map{it->[it[0],it[1],it[2],null]},
	  Channel.empty().collect().ifEmpty([]),
	  Channel.empty().collect().ifEmpty([])
        )
        chVersions = chVersions.mix(deeptoolsBamCoverage.out.versions)
      }

      // PROCESS : Qualimap
      if (!params.skipQC && !params.skipQualimap){
        qualimapRNAseq(
          chBamPassed,
          chGtf.collect()
        )
	chQualimapMqc = qualimapRNAseq.out.results.collect()
        chVersions = chVersions.mix(qualimapRNAseq.out.versions)
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
      if (!params.skipMarkDup){
        markdupFlow(
          chBamPassed,
          chGtf.collect()
        )
        chMarkDupMqc = markdupFlow.out.picardMetrics.collect()
        chDupradarMqc = markdupFlow.out.dupradarResults.collect()
        chVersions = chVersions.mix(markdupFlow.out.versions)
      }

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
          chGtf.collect()
        )
        chCounts = featureCountsFlow.out.counts
        chCountsTpm = featureCountsFlow.out.tpm
        chCountsMqc = featureCountsFlow.out.logs
        chVersions = chVersions.mix(featureCountsFlow.out.versions)
      } else if (params.counts == 'HTseqCounts'){
        htseqCountsFlow (
          chBamPassed,
          chGtf.collect()
        )
        chCounts = htseqCountsFlow.out.counts
        chCountsTpm = htseqCountsFlow.out.tpm
        chCountsMqc = htseqCountsFlow.out.logs
        chVersions = chVersions.mix(htseqCountsFlow.out.versions)
      } else if (params.counts == 'star'){
        starCountsFlow (
	  chBamPassed,
	  mappingStarFlow.out.counts,
	  mappingStarFlow.out.countsLogs,
          chGtf.collect()
        )
        chCounts = starCountsFlow.out.counts
        chCountsTpm = starCountsFlow.out.tpm
        chCountsMqc = starCountsFlow.out.logs
        chVersions = chVersions.mix(starCountsFlow.out.versions)
      } else if (params.counts == 'salmon'){
        salmonQuantFromBamFlow (
          mappingStarFlow.out.transcriptsBam,
	  chTranscriptsFasta,
	  chGtf
        )
        chCounts = salmonQuantFromBamFlow.out.countsGeneLengthScaled
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
	chSalmonIndex,
	chGtf
      )
      chCounts = salmonQuantFromFastqFlow.out.countsGeneLengthScaled
      chCountsTpm = salmonQuantFromFastqFlow.out.tpmGene
      chCountsMqc = salmonQuantFromFastqFlow.out.results
      chVersions = chVersions.mix(salmonQuantFromFastqFlow.out.versions)
      chAlignedBam = channel.empty()
      chBamPassed = channel.empty()
    }


    //******************************************
    // COUNTS-BASED QC
 
    // SUBWORKFLOW: gene counts qc
    if (!params.skipGeneCountsAnalysis && (params.counts || params.pseudoAligner)){
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
        chGtf.collect()
      )
      chVersions = chVersions.mix(stringtieFlow.out.versions)
      chgffCompareMqc = chgffCompareMqc.mix(stringtieFlow.out.mqc)
    }

    if ("scallop" in denovoTools){
      scallopFlow(
        chBamPassed,
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
      chAlignedBam
        .join(chBamPassed, remainder: true)
        .filter{it -> it[2] == null}
        .flatMap{ it -> it[0] + ": Poor alignment rate. Sample discarded !"}
        .set{chWarnMapping}

      strandnessFlow.out.strand
        .map{it[1]}
        .unique()
	.count()
	.filter{it > 1}
	.flatMap{"Samples with different strandness detected !"}
	.set{chWarnStrand}

     chWarnStrand
       .concat(chWarnMapping)
       .collectFile(name: 'warnings.txt', newLine: true)
       .set{chWarn}
    
      multiqc(
        customRunName,
        chSplan.collect(),
        chMetadata.ifEmpty([]),
        chMultiqcConfig.ifEmpty([]),
        chTrimmingMqc.collect().ifEmpty([]),
        chXengsortMqc.collect().ifEmpty([]),
        chFastqcMqc.ifEmpty([]),
        chrRNAMappingMqc.ifEmpty([]),
        chAlignedMqc.collect().ifEmpty([]),
        strandnessFlow.out.logs.collect().ifEmpty([]),
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
