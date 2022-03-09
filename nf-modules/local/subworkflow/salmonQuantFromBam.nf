/* 
 * FeatureCounts workflow
 */

include { salmonQuantFromBam } from '../process/salmonQuantFromBam'
include { salmonTx2gene } from '../process/salmonTx2gene'
include { salmonTxImport } from '../process/salmonTxImport'
include { mergeCounts} from '../process/mergeCounts'

workflow salmonQuantFromBamFlow {

  take:
  bam // Channel [val(meta), path(bam), path(bai)]
  trsFasta // Channel path(transcriptsFasta)
  gtf // Channel path(gtf)

  main:
  chVersions = Channel.empty()

  salmonQuantFromBam(
    bam,
    trsFasta.collect(),
    gtf.collect(),
  )
  chVersions = chVersions.mix(salmonQuantFromBam.out.versions)

  salmonTx2gene(
    salmonQuantFromBam.out.results.collect(),
    gtf.collect()
  )
  chVersions = chVersions.mix(salmonTx2gene.out.versions)

  salmonTxImport(
    salmonQuantFromBam.out.results.collect(), 
    salmonTx2gene.out.results.collect()
  )
  chVersions = chVersions.mix(salmonTxImport.out.versions)

  emit:
  results = salmonQuantFromBam.out.results
  tpmGene = salmonTxImport.out.tpmGene
  countsGene = salmonTxImport.out.countsGene
  countsGeneLengthScaled = salmonTxImport.out.countsGeneLengthScaled
  countsGeneScaled = salmonTxImport.out.countsGeneScaled
  tpmTranscript = salmonTxImport.out.tpmTranscript
  countsTranscript = salmonTxImport.out.countsTranscript
  versions = chVersions
}
