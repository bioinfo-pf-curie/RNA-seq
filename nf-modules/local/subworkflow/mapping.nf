/* 
 * define the data analysis workflow 
 */

/* 
 * include requires tasks 
 */
include { saveStrandness } from '../process/saveStrandness'
include { prepRseqc } from '../process/prepRseqc'
include { rseqc } from '../process/rseqc'

workflow rseqFlow {
    // required inputs
    take:
      reads 
      chRrnaAnnot
    // workflow implementation
    main:

    rRNAMapping(
      reads,
      chRrnaAnnot.collect()
    )

    // Update input channel
    chStarRawReads = Channel.empty()
    if( params.rrna && !params.skipRrna){
      chStarRawReads = rRNAMapping.out.chRrnaMappingRes
    } else {  
      chStarRawReads = chRawReadsStar
    }
     
    emit:
      chRrnaLogs      = rRNAMapping.out.logs
      chBowtieVersion = rRNAMapping.out.version 
}