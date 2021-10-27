/* 
 * FeatureCounts workflow
 */

include { featureCounts } from '../process/featureCounts'
include { mergeCounts} from '../process/mergeCounts'

workflow featureCountsFlow {

  take:
  bam // Channel [val(prefix), path(bam), path(bai)]
  gtf // Channel path(gtf)
  strandness

  main:
  tool = Channel.value("featureCounts")
  chVersions = Channel.empty()

  featureCounts(
    bam,
    gtf.collect(),
    strandness
  )
  chVersions = chVersions.mix(featureCounts.out.versions)

  mergeCounts(
    featureCounts.out.counts.collect(),
    gtf.collect(),
    strandness.collect(),
    tool
  )
  chVersions = chVersions.mix(mergeCounts.out.versions)

  emit:
  counts = mergeCounts.out.countsTable
  tpm = mergeCounts.out.tpmTable
  logs = featureCounts.out.logs
  versions = chVersions
}
