/* 
 * FeatureCounts workflow
 */

include { salmonQuantFromBam } from '../process/salmonQuantFromBam'
include { salmonTx2gene } from '../process/salmonTx2gene'
include { salmonTxImport } from '../process/salmonTxImport'
include { mergeCounts} from '../process/mergeCounts'

workflow salmonQuantFromBamFlow {

  take:
  bam // Channel [val(prefix), path(bam), path(bai)]
  strandness // Channel val(strandness)
  trsFasta // Channel path(transcriptsFasta)
  gtf // Channel path(gtf)

  main:
  //tool = Channel.value("featureCounts")
  chVersions = Channel.empty()

  salmonQuantFromBam(
    bam.join(strandness),
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
