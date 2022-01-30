/* 
 * Star counts workflow
 */

include { mergeCounts} from '../process/mergeCounts'

workflow starCountsFlow {

  take:
  starBams
  starCounts // Channel path(starCounts)
  starCountsLogs // Channel path(starCountsLogs)
  strandness // Channel val(strandness)
  gtf // Channel path(gtf)

  main:
  tool = Channel.value("star")
  chVersions = Channel.empty()

  // Phase counts with aligned data if any has been skipped
  starBams
    .join(starCounts)
    .map(it -> [it[0], it[3]])
    .set{chStarCounts}

  starBams
    .join(starCountsLogs)
    .map(it -> [it[0], it[3]])
    .set{chStarCountsLogs}                                                                                                                                                                             

  // Merge and order counts
  chStarCounts
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
  logs = chStarCountsLogs
  versions = chVersions
}
