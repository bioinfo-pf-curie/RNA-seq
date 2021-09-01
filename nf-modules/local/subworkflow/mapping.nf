/* 
 * define the data analysis workflow 
 */

/* 
 * include requires tasks 
 */
include { rRNAMapping } from '../process/rRNAMapping'
include { star } from '../process/star'
include { starSort } from '../process/starSort'
include { makeHisatSplicesites } from '../process/makeHisatSplicesites'
include { hisat2Align } from '../process/hisat2Align'
include { hisat2Sort } from '../process/hisat2Sort'

workflow rseqFlow {
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
      chAlignmentLogs
      chStarLogCounts       = star.out.starLogCounts
      chStarLog
      chStarCountsTo        = star.out.starCountsTo
      chStarVersion         
      chSamtoolsVersionSort
      chBam
      chHisat2Version       = hisat2Align.out.version
}