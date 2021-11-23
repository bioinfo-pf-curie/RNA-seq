/* 
 * Strandness workflow
 */

include { saveStrandness } from '../process/saveStrandness'
include { rseqcPrep } from '../process/rseqcPrep'
include { rseqc } from '../process/rseqc'

workflow strandnessFlow {

  take:
  reads // Channel [prefix, [reads]]
  bed12 // Channel path(bed12)

  main:
  chVersions = Channel.empty()

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
    strandnessOutputFiles = saveStrandness.out.savedStrandness
  }

  // auto detection of strandness status
  if (params.stranded == 'auto' && params.bed12){
    rseqcPrep(
      reads
    )
    chVersions = chVersions.mix(rseqcPrep.out.versions)

    rseqc(
      rseqcPrep.out.bamRseqc,
      bed12.collect()
    )
    chVersions = chVersions.mix(rseqc.out.versions)
    strandnessResults = rseqc.out.strandnessResults.splitCsv()
    strandnessOutputFiles = rseqc.out.rseqcResults
  } else {
    exit 1, "Cannot detect strandness without a bed12 annotation file. Please use the --stranded option."
  }
      
  emit:
  strandnessResults
  strandnessOutputFiles
  versions = chVersions
}
