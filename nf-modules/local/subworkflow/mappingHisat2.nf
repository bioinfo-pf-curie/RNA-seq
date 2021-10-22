/* 
 * define the data analysis workflow 
 */

/* 
 * include requires tasks 
 */
include { makeHisatSplicesites } from '../process/makeHisatSplicesites'
include { hisat2Align } from '../process/hisat2Align'
include { samtoolsSort } from '../process/samtoolsSort'
include { samtoolsIndex } from '../process/samtoolsIndex'

workflow mappingHisat2Flow {

  // required inputs
  take:
  reads
  index
  gtf
  strandness

  // workflow implementation
  main:

  makeHisatSplicesites(
    gtf
  )

  hisat2Align(
    reads,
    index.collect(),
    makeHisatSplicesites.out.alignmentSplicesites,
    strandness
  )

  samtoolsSort(
    hisat2Align.out.bam
  )

  samtoolsIndex(
    samtoolsSort.out.bam
  )

  samtoolsFlagstat(
   samtoolsSort.out.bam
  )

  emit:
  bam = samtoolsSort.out.bam
  bai = samtoolsIndex.out.bai
  logs = hisat2Align.out.logs
  flagstat = samtoolsFlagstat.out.stats
  chHisat2Version = hisat2Align.out.version
  chSamtoolsVersionSort = samtoolsSort.out.version
}
