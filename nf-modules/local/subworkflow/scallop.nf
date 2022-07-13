/* 
 * Scallop workflow
 */

include { scallop } from '../../common/process/scallop/scallop'
include { gffcompare } from '../../common/process/gffcompare/gffcompare'

workflow scallopFlow {

  take:
  bam // Channel [val(meta), path(bam), path(bai)]
  gtf // Channel path(gtf)

  main:
  chVersions = Channel.empty()

  scallop(
    bam
  )
  chVersions = chVersions.mix(scallop.out.versions)

  gffcompare(
    scallop.out.transcriptGtf, 
    gtf.collect()
  )
  chVersions = chVersions.mix(gffcompare.out.versions)

  emit:
  gffCompareResults = gffcompare.out.results
  mqc = gffcompare.out.mqc
  versions = chVersions
}
