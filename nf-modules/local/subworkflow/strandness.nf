/* 
 * Strandness workflow
 */

include { saveStrandness } from '../process/saveStrandness'
include { bowtie2 } from '../../common/process/bowtie2/bowtie2'
include { rseqc } from '../../common/process/rseqc/rseqc'

workflow strandnessFlow {

  take:
  reads // Channel [meta, [reads]]
  bed12 // Channel path(bed12)
  bowtie2Index // Channel path(bowtie2Index)

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
      .map { meta, reads ->
         def key = params.stranded
         return [meta.id, key]
      }
      .set { strandnessResults }
    strandnessOutputFiles = saveStrandness.out.savedStrandness
  }

  // auto detection of strandness status
  if (params.stranded == 'auto' && params.bed12){
    bowtie2(
      reads,
      bowtie2Index.collect()
    )
    chVersions = chVersions.mix(bowtie2.out.versions)

    rseqc(
      bowtie2.out.bam,
      bed12.collect()
    )
    chVersions = chVersions.mix(rseqc.out.versions)
    strandnessResults = rseqc.out.strandnessResults.splitCsv()
    strandnessOutputFiles = rseqc.out.rseqcResults
  }
      
  emit:
  strand = strandnessResults
  logs = strandnessOutputFiles
  versions = chVersions
}
