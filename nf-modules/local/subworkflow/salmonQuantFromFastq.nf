/* 
 * FeatureCounts workflow
 */

include { salmonQuant as salmonQuantFromFastq } from '../../common/process/salmon/salmonQuant'
include { salmonTx2gene } from '../process/salmonTx2gene'
include { salmonTxImport } from '../process/salmonTxImport'
include { mergeCounts} from '../process/mergeCounts'

workflow salmonQuantFromFastqFlow {

  take:
  reads // Channel [val(meta), reads]
  index // Channel path(salmonIndex)
  gtf // Channel path(gtf)

  main:
  chVersions = Channel.empty()

  salmonQuantFromFastq(
    reads,
    index.collect(),
    Channel.value([]),
    gtf.collect(),
  )
  chVersions = chVersions.mix(salmonQuantFromFastq.out.versions)

  salmonTx2gene(
    salmonQuantFromFastq.out.results.collect(),
    gtf.collect()
  )
  chVersions = chVersions.mix(salmonTx2gene.out.versions)

  salmonTxImport(
    salmonQuantFromFastq.out.results.collect(), 
    salmonTx2gene.out.results.collect()
  )
  chVersions = chVersions.mix(salmonTxImport.out.versions)

  emit:
  results = salmonQuantFromFastq.out.results
  tpmGene = salmonTxImport.out.tpmGene
  countsGene = salmonTxImport.out.countsGene
  countsGeneLengthScaled = salmonTxImport.out.countsGeneLengthScaled
  countsGeneScaled = salmonTxImport.out.countsGeneScaled
  tpmTranscript = salmonTxImport.out.tpmTranscript
  countsTranscript = salmonTxImport.out.countsTranscript
  versions = chVersions
}
