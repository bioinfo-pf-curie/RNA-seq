/* 
 * HTSeq Workflow
 */

include { htseqCounts } from '../../common/process/htseq/htseqCounts'
include { mergeCounts} from '../process/mergeCounts'

workflow htseqCountsFlow {

  take:
  bam // Channel [val(meta), path(bam), path(bai)]
  gtf // Channel path(gtf)

  main:
  tool = Channel.value("HTseqCounts")
  chVersions = Channel.empty()
 
  htseqCounts(
    bam,
    gtf.collect()
  )
  chVersions = chVersions.mix(htseqCounts.out.versions)

  // Merge and order counts
  htseqCounts.out.counts
    .multiMap { it ->
      counts: it[1]
      strand: it[0].strandness
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
