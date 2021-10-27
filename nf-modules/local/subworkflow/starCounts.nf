/* 
 * Star counts workflow
 */

include { mergeCounts} from '../process/mergeCounts'

workflow starCountsFlow {

  take:
  starCounts // Channel path(starCounts)
  starCountsLogs // Channel path(starCountsLogs)
  gtf // Channel path(gtf)
  strandness

  main:
  tool = Channel.value("star")
  chVersions = Channel.empty()

  mergeCounts(
    starCounts.collect(),
    gtf.collect(),
    strandness.collect(),
    tool
  )
  chVersions = chVersions.mix(mergeCounts.out.versions)

  emit:
  counts = mergeCounts.out.countsTable
  tpm = mergeCounts.out.tpmTable
  logs = starCountsLogs
  versions = chVersions
}
