/* 
 * MarkDuplicates
 */

include { markDuplicates } from '../../common/process/markDuplicates'
include { samtoolsIndex } from '../../common/process/samtoolsIndex'
include { dupradar } from '../process/dupradar'

workflow markdupFlow {

  take:
  bam
  gtf

  main:
  chVersions = Channel.empty()

  markDuplicates(
    bam
  )
  chVersions = chVersions.mix(markDuplicates.out.versions)

  samtoolsIndex(
    markDuplicates.out.bam
  )
  chVersions = chVersions.mix(samtoolsIndex.out.versions)

  if (!params.skipDupradar){
    dupradar (
      markDuplicates.out.bam,
      gtf.collect(),
    )
    chDupradarResults = dupradar.out.results
    chVersions = chVersions.mix(dupradar.out.versions)
  }else{
    chDupradarResults = Channel.empty()
  }

  emit:
  bam = markDuplicates.out.bam.join(samtoolsIndex.out.bai)
  picardMetrics = markDuplicates.out.metrics
  dupradarResults = chDupradarResults
  versions = chVersions
}