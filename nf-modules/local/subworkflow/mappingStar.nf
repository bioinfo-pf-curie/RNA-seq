/* 
 * STAR Workflow
 */

include { starAlign } from '../process/starAlign'
include { samtoolsSort } from '../process/samtoolsSort'
include { samtoolsIndex } from '../process/samtoolsIndex'
include { samtoolsFlagstat } from '../process/samtoolsFlagstat'

workflow mappingStarFlow {

  take:
  reads // Channel [val(prefix), [reads]]
  index // Channel path(index)
  gtf // Channel path(gtf)

  main:
  
  chVersions = Channel.empty()

  starAlign(
    reads,
    index.collect(),
    gtf.collect().ifEmpty([])
  )
  chVersions = chVersions.mix(starAlign.out.versions)

  samtoolsSort(
    starAlign.out.bam
  )
  chVersions = chVersions.mix(samtoolsSort.out.versions)

  samtoolsIndex(
    samtoolsSort.out.bam
  )
  chVersions = chVersions.mix(samtoolsIndex.out.versions)

  samtoolsFlagstat(
    samtoolsSort.out.bam
  )
  chVersions = chVersions.mix(samtoolsFlagstat.out.versions)

  emit:
  bam = samtoolsSort.out.bam
  bai = samtoolsIndex.out.bai
  flagstat = samtoolsFlagstat.out.stats
  logs = starAlign.out.logs
  counts = starAlign.out.counts
  countsLogs = starAlign.out.countsLogs
  versions = chVersions
}
