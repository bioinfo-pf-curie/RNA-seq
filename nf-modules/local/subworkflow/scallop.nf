/* 
 * Scallop workflow
 */

include { scallop } from '../process/scallop'
include { gffcompare } from '../process/gffcompare'

workflow scallopFlow {

  take:
  bam // Channel [val(prefix), path(bam), path(bai)]
  strandness // Channel val(strandness)
  gtf // Channel path(gtf)

  main:
  chVersions = Channel.empty()

  scallop(
    bam.join(strandness)
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
