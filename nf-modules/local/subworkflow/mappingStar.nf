/* 
 * STAR Workflow
 */

include { starAlign } from '../../common/process/star/starAlign'
include { samtoolsSort as samtoolsSortGenome } from '../../common/process/samtools/samtoolsSort'
include { samtoolsSort as samtoolsSortTranscriptome } from '../../common/process/samtools/samtoolsSort'
include { samtoolsIndex as samtoolsIndexGenome } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsIndex as samtoolsIndexTranscriptome } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsFlagstat } from '../../common/process/samtools/samtoolsFlagstat'

workflow mappingStarFlow {

  take:
  reads // Channel [val(meta), [reads]]
  index // Channel path(index)
  gtf // Channel path(gtf)

  main:
  
  chVersions = Channel.empty()

  starAlign(
    reads,
    index.collect(),
    gtf.collect().ifEmpty([])
  )
  chVersions = chVersions.mix(starAlign.out.versions)

  samtoolsSortGenome(
    starAlign.out.bam
  )
  chVersions = chVersions.mix(samtoolsSortGenome.out.versions)

  samtoolsIndexGenome(
    samtoolsSortGenome.out.bam
  )
  chVersions = chVersions.mix(samtoolsIndexGenome.out.versions)

  samtoolsSortTranscriptome(
    starAlign.out.transcriptomeBam
  )
  chVersions = chVersions.mix(samtoolsSortTranscriptome.out.versions)

  samtoolsIndexTranscriptome(
    samtoolsSortTranscriptome.out.bam
  )
  chVersions = chVersions.mix(samtoolsIndexTranscriptome.out.versions)

  samtoolsFlagstat(
    samtoolsSortGenome.out.bam
  )
  chVersions = chVersions.mix(samtoolsFlagstat.out.versions)

  emit:
  genomeBam = samtoolsSortGenome.out.bam.join(samtoolsIndexGenome.out.bai)
  transcriptomeBam = samtoolsSortTranscriptome.out.bam.join(samtoolsIndexTranscriptome.out.bai)
  flagstat = samtoolsFlagstat.out.stats
  logs = starAlign.out.logs
  counts = starAlign.out.counts
  countsLogs = starAlign.out.countsLogs
  versions = chVersions
}
