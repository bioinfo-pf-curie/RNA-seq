/* 
 * HTSeq Workflow
 */

include { htseqCounts } from '../process/htseqCounts'
include { mergeCounts} from '../process/mergeCounts'

workflow htseqCountsFlow {

  take:
  bam // Channel [val(prefix), path(bam), path(bai)]
  gtf // Channel path(gtf)
  strandness

  main:
  tool = Channel.value("HTseqCounts")
  chVersions = Channel.empty()
 
  htseqCounts(
    bam,
    gtf.collect(),
    strandness
  )
  chVersions = chVersions.mix(htseqCounts.out.versions)

  mergeCounts(
    htseqCounts.out.counts.collect(),
    gtf.collect(),
    strandness.collect(),
    tool
  )
  chVersions = chVersions.mix(mergeCounts.out.versions)

  emit:
  counts = mergeCounts.out.countsTable
  tpm = mergeCounts.out.tpmTable
  logs = htseqCounts.out.logs
  versions = chVersions
}
