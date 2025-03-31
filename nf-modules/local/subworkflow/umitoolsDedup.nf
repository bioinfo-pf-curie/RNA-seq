/* 
 * UMItools dedup
 */

include { umitoolsDedup } from '../../common/process/umitools/umitoolsDedup'
include { samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'

workflow umitoolsDedupFlow {

  take:
  bam // Channel [val(meta), path(bam), path(bai)]

  main:
  chVersions = Channel.empty()

  umitoolsDedup(
    bam,
    Channel.of(true).collect()
  )
  chVersions = chVersions.mix(umitoolsDedup.out.versions)

  samtoolsIndex(
    umitoolsDedup.out.bam
  )
  chVersions = chVersions.mix(samtoolsIndex.out.versions)

  emit:
  bam = umitoolsDedup.out.bam.join(samtoolsIndex.out.bai)
  logs = umitoolsDedup.out.log
  versions = chVersions
}
