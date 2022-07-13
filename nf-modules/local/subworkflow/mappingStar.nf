/* 
 * STAR Workflow
 */

include { starAlign } from '../../common/process/star/starAlign'
include { samtoolsSort } from '../../common/process/samtools/samtoolsSort'
include { samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsFlagstat } from '../../common/process/samtools/samtoolsFlagstat'

workflow mappingStarFlow {

  take:
  reads // Channel [val(meta), [reads]]
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
  transcriptsBam = starAlign.out.transcriptsBam
  bai = samtoolsIndex.out.bai
  flagstat = samtoolsFlagstat.out.stats
  logs = starAlign.out.logs
  counts = starAlign.out.counts
  countsLogs = starAlign.out.countsLogs
  versions = chVersions
}
