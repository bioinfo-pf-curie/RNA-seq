/* 
 * FeatureCounts workflow
 */

include { featureCounts } from '../../common/process/featureCounts/featureCounts'
include { mergeCounts} from '../process/mergeCounts'

workflow featureCountsFlow {

  take:
  bam // Channel [val(prefix), path(bam), path(bai)]
  gtf // Channel path(gtf)

  main:
  tool = Channel.value("featureCounts")
  chVersions = Channel.empty()

  featureCounts(
    bam.combine(gtf)
  )
  chVersions = chVersions.mix(featureCounts.out.versions)

  // merge and order counts and strand
  featureCounts.out.counts
    .multiMap { it ->
      counts: it[1]
      strand: it[0].id + "\t" + it[1] + "\t" + it[0].strandness
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
