/* 
 * FeatureCounts workflow
 */

include { salmonQuantFromBam } from '../process/salmonQuantFromBam'
include { salmonTx2gene } from '../process/salmonTx2gene'
include { salmonTxImport } from '../process/salmonTxImport'
include { mergeCounts} from '../process/mergeCounts'

workflow salmonCountsFlow {

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
    salmonQuantFromBam.out.results,
    gtf.collect()
  )
  chVersions = chVersions.mix(salmonTx2gene.out.versions)

  salmonTxImport(
    salmonQuantFromBam.out.results.join(salmonTx2gene.out.results)
  )
  chVersions = chVersions.mix(salmonTxImport.out.versions)


  // merge and order counts and strand
  //featureCounts.out.counts
  //  .join(strandness)
  //  .multiMap { it ->
  //    counts: it[1]
  //    strand: it[2]
  // }.set{chCountsAndStrand}

  //mergeCounts(
  //  chCountsAndStrand.counts.collect(),
  //  chCountsAndStrand.strand.collect(),
  //  gtf.collect(),
  //  tool
  //)
  //chVersions = chVersions.mix(mergeCounts.out.versions)

  emit:
  counts = Channel.empty()
  tpm = Channel.empty()
  logs = Channel.empty()
  versions = chVersions
}
