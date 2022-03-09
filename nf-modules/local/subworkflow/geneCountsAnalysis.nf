/* 
 * Gene counts analysis workflow
 */

include { geneSaturation } from '../process/geneSaturation'
include { getCountsPerGeneType} from '../process/getCountsPerGeneType'
include { exploratoryAnalysis } from '../process/exploratoryAnalysis'

workflow geneCountsAnalysisFlow {

  take:
  counts
  tpm
  gtf
  pcaHeader
  heatmapHeader


  main:
  chVersions = Channel.empty()

  geneSaturation(
    counts.collect()
  )
  chVersions = chVersions.mix(geneSaturation.out.versions)

  getCountsPerGeneType(
    tpm.collect(),
    gtf.collect()
  )
  chVersions = chVersions.mix(getCountsPerGeneType.out.versions)

  exploratoryAnalysis(
    counts.collect(),
    tpm.collect(),
    pcaHeader,
    heatmapHeader 
  )
  chVersions = chVersions.mix(exploratoryAnalysis.out.versions)

  emit:
  geneSaturationResults = geneSaturation.out.results
  countsPerGenetype = getCountsPerGeneType.out.results
  expAnalysisResults = exploratoryAnalysis.out.results
  versions = chVersions
}
