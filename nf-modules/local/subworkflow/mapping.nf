/* 
 * define the data analysis workflow 
 */

/* 
 * include requires tasks 
 */
include { checkStarLog1 } from './functions'
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
      chSamtoolsVersionSort = starSort.out.samtoolsVersionSort

      // Filter removes all 'aligned' channels that fail the check
      starSort.out.starAligned
        .filter { logs, bams -> checkStarLog(logs) }
        .map { logs, bams -> bams }
        .dump (tag:'starbams')
        .set { chBam }
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
        makeHisatSplicesites.out.alignmentSplicesites
        chStrandnessResults
      )
      chAlignmentLogs = hisat2Align.out.alignmentLogs
      hisat2Sort(
        hisat2Align.out.hisat2Bam
      )
      chBam = hisat2Sort.out.chBam
      chSamtoolsVersionSort = hisat2Sort.out.samtoolsVersionSort
    } 

    emit:
      chRrnaLogs            = rRNAMapping.out.logs
      chBowtieVersion       = rRNAMapping.out.version 
      chStarLogCounts       = star.out.starLogCounts
      chStarCounts          = star.out.starCounts        
      chHisat2Version       = hisat2Align.out.version
      chBam
      chAlignmentLogs
      chStarLog
      chStarVersion 
      chSamtoolsVersionSort
}
