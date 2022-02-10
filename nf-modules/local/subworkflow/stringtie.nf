/* 
 * Stringtie workflow
 */

include { stringtie } from '../process/stringtie'
include { stringtieMerge } from '../process/stringtieMerge'
include { gffcompare } from '../process/gffcompare'

workflow stringtieFlow {

  take:
  bam // Channel [val(prefix), path(bam), path(bai)]
  strandness // Channel val(strandness)
  gtf // Channel path(gtf)

  main:
  chVersions = Channel.empty()

  stringtie(
    bam.join(strandness),
    gtf.collect(),
  )
  chVersions = chVersions.mix(stringtie.out.versions)

  stringtieMerge(
    stringtie.out.transcriptGtf.map{it[1]}.collect(),
    gtf.collect()
  )
  chVersions = chVersions.mix(stringtieMerge.out.versions)

  gffcompare(
    stringtieMerge.out.mergedGtf, 
    gtf.collect()
  )
  chVersions = chVersions.mix(gffcompare.out.versions)

  emit:
  combinedGtf = stringtieMerge.out.mergedGtf
  gffCompareResults = gffcompare.out.results
  mqc = gffcompare.out.mqc
  versions = chVersions
}
