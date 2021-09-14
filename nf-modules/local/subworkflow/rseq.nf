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
      chRawReads 
      chBedRseqc
    // workflow implementation
    main:
      // Strandness
      saveStrandness(
        chRawReads
      )
      chStrandedResults = Channel.empty()
      // User defined
      if (params.stranded == 'reverse' || params.stranded == 'forward' || params.stranded == 'no'){
        chRawReadsStrandness
        .map { file ->
           def key = params.stranded
           return tuple(key)
        }
        .set { chStrandedResults }
        chBowtie2Version = Channel.empty()
        chRseqcVersionInferExperiment = Channel.empty()  
      }else if (params.stranded == 'auto' && params.bed12){
        // auto     
        prepRseqc(
          chRawReads
        )
        chBowtie2Version  = prepRseqc.out.bowtie2Version
        rseqc(
          prepRseqc.out.bamRseqc,
          chBedRseqc.collect()
        )
        chRseqcVersionInferExperiment = rseqc.out.version
        chStrandedResults = rseqc.out.strandedResults
      }
      chStrandnessResults = Channel.empty()
      if (params.stranded == 'auto' && params.bed12){
        chStrandnessResults = rseqc.out.rseqcResults
      }else{
        chStrandnessResults = saveStrandness.out.savedStrandness
      }
      
    emit:
      chBowtie2Version
      chStrandedResults
      chRseqcVersionInferExperiment    
      chStrandnessResults
}
