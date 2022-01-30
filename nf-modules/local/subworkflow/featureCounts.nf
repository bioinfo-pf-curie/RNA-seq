/* 
 * FeatureCounts workflow
 */

include { featureCounts } from '../process/featureCounts'
include { mergeCounts} from '../process/mergeCounts'

workflow featureCountsFlow {

  take:
  bam // Channel [val(prefix), path(bam), path(bai)]
  strandness // Channel val(strandness)
  gtf // Channel path(gtf)

  main:
  tool = Channel.value("featureCounts")
  chVersions = Channel.empty()

  featureCounts(
    bam.join(strandness),
    gtf.collect(),
  )
  chVersions = chVersions.mix(featureCounts.out.versions)

  // merge and order counts and strand
  featureCounts.out.counts
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
  logs = featureCounts.out.logs
  versions = chVersions
}
