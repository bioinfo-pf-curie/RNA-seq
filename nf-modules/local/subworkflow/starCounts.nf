/* 
 * Star counts workflow
 */

include { mergeCounts} from '../process/mergeCounts'

workflow starCountsFlow {

  take:
  starCounts // Channel path(starCounts)
  starCountsLogs // Channel path(starCountsLogs)
  strandness // Channel val(strandness)
  gtf // Channel path(gtf)

  main:
  tool = Channel.value("star")
  chVersions = Channel.empty()

  // Merge and order counts
  starCounts
    .join(strandness)
    .multiMap { it ->
      counts: it[1]
      strand: it[2]
   }.set{chCountsAndStrand}

  mergeCounts(
    chCountsAndStrand.counts.collect(),
    chCountsAndStrand.strand.collect(),
    gtf.collect(),
    tool
  )
  chVersions = chVersions.mix(mergeCounts.out.versions)

  emit:
  counts = mergeCounts.out.countsTable
  tpm = mergeCounts.out.tpmTable
  logs = starCountsLogs
  versions = chVersions
}
