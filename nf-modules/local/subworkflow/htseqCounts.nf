/* 
 * HTSeq Workflow
 */

include { htseqCounts } from '../process/htseqCounts'
include { mergeCounts} from '../process/mergeCounts'

workflow htseqCountsFlow {

  take:
  bam // Channel [val(prefix), path(bam), path(bai)]
  strandness // Channel(strandness)
  gtf // Channel path(gtf)

  main:
  tool = Channel.value("HTseqCounts")
  chVersions = Channel.empty()
 
  htseqCounts(
    bam.join(strandness),
    gtf.collect()
  )
  chVersions = chVersions.mix(htseqCounts.out.versions)

  // Merge and order counts
  htseqCounts.out.counts
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
  logs = htseqCounts.out.logs
  versions = chVersions
}
