/* 
 * define the data analysis workflow 
 */

/* 
 * include requires tasks 
 */
include { saveStrandness } from '../process/saveStrandness'
include { rseqcPrep } from '../process/rseqcPrep'
include { rseqc } from '../process/rseqc'

workflow strandnessFlow {
  // required inputs
  take:
  reads
  bed12

  // workflow implementation
  main:
  // if --stranded option is specified by the user
  if (params.stranded == 'reverse' || params.stranded == 'no' || params.stranded == 'yes'){

    // save it in a text file for MultiQC report
    saveStrandness(
      reads
    )

    // emit a channel with strandness information
    reads
      .map { file ->
         def key = params.stranded
         return tuple(key)
      }
      .set { strandnessResults }

    bowtie2Version = Channel.empty()
    rseqcVersionInferExperiment = Channel.empty()
    strandnessOutputFiles = saveStrandness.out.savedStrandness
  }

  // auto detection of strandness status
  if (params.stranded == 'auto' && params.bed12){

    rseqcPrep(
      reads
    )

    rseqc(
      rseqcPrep.out.bamRseqc,
      bed12.collect()
    )

    bowtie2Version  = rseqcPrep.out.bowtie2Version 
    rseqcVersionInferExperiment = rseqc.out.version
    strandnessResults = rseqc.out.strandnessResults
    strandnessOutputFiles = rseqc.out.rseqcResults
  } else {
    exit 1, "Cannot detect strandness without a bed12 annotation file. Please use the --stranded option."
  }
      
  emit:
  strandnessResults
  strandnessOutputFiles
  bowtie2Version
  rseqcVersionInferExperiment    
}
