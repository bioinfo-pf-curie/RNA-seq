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
    scallop.out.transcriptGtf.map{it[1]}.collect(), 
    gtf.collect(),
    Channel.value("scallop")
  )
  chVersions = chVersions.mix(gffcompare.out.versions)

  emit:
  combinedGtf = gffcompare.out.combinedGtf
  gffCompareResults = gffcompare.out.results
  versions = chVersions
}
