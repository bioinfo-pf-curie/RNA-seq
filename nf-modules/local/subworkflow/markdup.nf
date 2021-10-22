/* 
 * define the data analysis workflow 
 */

/* 
 * include requires tasks 
 */
include { markDuplicates } from '../process/markDuplicates'
include { samtoolsIndex } from '../process/samtoolsIndex'
include { dupradar } from '../process/dupradar'

workflow markdupFlow {
  // required inputs
  take:
  bam
  gtf
  strandness

  // workflow implementation
  main:
  // Duplicates
  markDuplicates(
    bam
  )

  samtoolsIndex(
    markDuplicates.out.bam
  )

  if (!params.skipDupradar){
    dupradar (
      markDuplicates.out.bam,
      gtf.collect(),
      strandness
    )
  }

  emit:
  bam = markDuplicates.out.bam
  bai = samtoolsIndex.out.bai
  metrics = markDuplicates.out.metrics
  chPicardVersion = markDuplicates.out.version
}
