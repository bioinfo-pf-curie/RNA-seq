/* 
 * define the data analysis workflow 
 */

/* 
 * include requires tasks 
 */
include { checkStarLog } from './functions'
include { rRNAMapping } from '../process/rRNAMapping'
include { star } from '../process/star'
include { starSort } from '../process/starSort'
include { makeHisatSplicesites } from '../process/makeHisatSplicesites'
include { hisat2Align } from '../process/hisat2Align'
include { hisat2Sort } from '../process/hisat2Sort'

workflow mappingFlow {
    // required inputs
    take:
      chRawReads
      chRrnaAnnot
      chStarIndex
      chGtf
      chHisat2Index
      chStrandnessResults

    // workflow implementation
    main:
    // rRNA mapping 
    rRNAMapping(
      chRawReads,
      chRrnaAnnot.collect()
    )

    // Reads mapping
    // From nf-core
    // Function that checks the alignment rate of the STAR output
    // and returns true if the alignment passed and otherwise false
    skippedPoorAlignment = []

    // Update input channel
    if( params.rrna && !params.skipRrna){
      chRawReads = rRNAMapping.out.chRrnaMappingRes
    } 

    chAlignmentLogs = Channel.empty()
    chBam = Channel.empty()
    // STAR
    if(params.aligner == 'star'){
      chHisat2Version = Channel.empty()
      star(
        chRawReads,
        chStarIndex.collect(),
        chGtf.collect().ifEmpty([])
      )
      chAlignmentLogs = star.out.alignmentLogs
      chStarLog       = star.out.starLog
      chStarVersion   = star.out.version

      starSort(
        star.out.chStarSam
      )
      chSamtoolsVersionSort = starSort.out.version

      // Filter removes all 'aligned' channels that fail the check
      starSort.out.starAligned
       .filter { logs, bams -> checkStarLog(logs) }
       .map { logs, bams -> bams }
       .dump (tag:'starbams')
       .set { chBam }
      // chBam = starSort.out.starAligned
    }

    // HiSat2 
    if(params.aligner == 'hisat2'){
      chStarLog = Channel.empty()
      chStarVersion = Channel.empty() 
      makeHisatSplicesites(
        chGtf
      )

      hisat2Align(
        chRawReads,
        chHisat2Index.collect(),
        makeHisatSplicesites.out.alignmentSplicesites,
        chStrandnessResults
      )
      chAlignmentLogs = hisat2Align.out.alignmentLogs
      chHisat2Version = hisat2Align.out.version

      hisat2Sort(
        hisat2Align.out.hisat2Bam
      )
      chBam = hisat2Sort.out.bam
      chSamtoolsVersionSort = hisat2Sort.out.version
    } 

    emit:
      chRrnaLogs            = rRNAMapping.out.logs
      chBowtieVersion       = rRNAMapping.out.version
      chStarLogCounts       = star.out.starLogCounts
      chStarCounts          = star.out.starCounts
      chHisat2Version
      chBam
      chAlignmentLogs
      chStarLog
      chStarVersion
      chSamtoolsVersionSort
      //skippedPoorAlignment
}
