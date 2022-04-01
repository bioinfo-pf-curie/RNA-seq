/* 
 * HiSat2 Worflow
 */

include { makeHisatSplicesites } from '../process/makeHisatSplicesites'
include { hisat2Align } from '../process/hisat2Align'
include { samtoolsSort } from '../../common/process/samtoolsSort'
include { samtoolsIndex } from '../../common/process/samtoolsIndex'
include { samtoolsFlagstat } from '../../common/process/samtoolsFlagstat'

workflow mappingHisat2Flow {

  take:
  reads
  index
  gtf

  main:
  chVersions = Channel.empty()

  makeHisatSplicesites(
    gtf
  )
  chVersions = chVersions.mix(makeHisatSplicesites.out.versions)

  hisat2Align(
    reads,
    index.collect(),
    makeHisatSplicesites.out.alignmentSplicesites.collect(),
  )
  chVersions = chVersions.mix(hisat2Align.out.versions)

  samtoolsSort(
    hisat2Align.out.bam
  )
  chVersions = chVersions.mix(samtoolsSort.out.versions)

  samtoolsIndex(
    samtoolsSort.out.bam
  )
  chVersions = chVersions.mix(samtoolsIndex.out.versions)

  samtoolsFlagstat(
   samtoolsSort.out.bam
  )
  chVersions = chVersions.mix(samtoolsFlagstat.out.versions)

  emit:
  bam = samtoolsSort.out.bam
  bai = samtoolsIndex.out.bai
  logs = hisat2Align.out.logs
  flagstat = samtoolsFlagstat.out.stats
  versions = chVersions
}