/* 
 * define the data analysis workflow 
 */

/* 
 * include requires tasks 
 */
include { checkStarLog } from './functions'
include { starAlign } from '../process/starAlign'
include { samtoolsSort } from '../process/samtoolsSort'
include { samtoolsIndex } from '../process/samtoolsIndex'
include { samtoolsFlagstat } from '../process/samtoolsFlagstat'

workflow mappingStarFlow {

  // required inputs
  take:
  reads
  index
  gtf

  // workflow implementation
  main:

  starAlign(
    reads,
    index.collect(),
    gtf.collect().ifEmpty([])
  )

  samtoolsSort(
    starAlign.out.bam
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
  flagstat = samtoolsFlagstat.out.stats
  logs = starAlign.out.logs
  counts = starAlign.out.counts
  countsLogs = starAlign.out.countsLogs
  chStarVersion = starAlign.out.version
  chSamtoolsVersionSort = samtoolsSort.out.version
}
