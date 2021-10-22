/* 
 * define the data analysis workflow 
 */

/* 
 * include requires tasks 
 */
include { featureCounts } from '../process/featureCounts'
include { HTseqCounts } from '../process/HTseqCounts'
include { mergeCounts} from '../process/mergeCounts'
include { geneSaturation } from '../process/geneSaturation'
include { getCountsPerGeneType} from '../process/getCountsPerGeneType'
include { exploratoryAnalysis } from '../process/exploratoryAnalysis'

workflow countsFlow {
    // required inputs
    take:
      chBam
      chGtf
      chStrandedResults
      chStarCounts
      chStarLogCounts
      chPcaHeader
      chHeatmapHeader

    // workflow implementation
    main:
      // Counts

     chCounts = Channel.empty()
     chCountsLogs = Channel.empty()
     chFeaturecountsVersion = Channel.empty() 
     chHtseqVersion = Channel.empty()
     if( params.counts == 'featureCounts' ){
       featureCounts(
         chBam,
         chGtf.collect(),
         chStrandedResults
       )
       chCounts = featureCounts.out.counts
       chCountsLogs = featureCounts.out.logs
       chFeaturecountsVersion = featureCounts.out.version 
     } else if (params.counts == 'HTseqCounts'){
         HTseqCounts (
           chBam,
           chGtf.collect(),
           chStrandedResults
         )
       chCounts = HTseqCounts.out.counts
       chCountsLogs = HTseqCounts.out.counts
       chHtseqVersion = HTseqCounts.out.version
     } else if (params.counts == 'star'){
         chCounts = chStarCounts
         chCountsLogs = chStarLogCounts
     }

      mergeCounts(
        chCounts.collect(),
        chGtf.collect(),
        chStrandedResults.collect()
      )

      geneSaturation(
        mergeCounts.out.counts.collect()
      )

      getCountsPerGeneType(
        mergeCounts.out.tpmCounts,
        chGtf.collect()
      )

      exploratoryAnalysis(
        mergeCounts.out.counts.collect(),
        mergeCounts.out.tpmCounts.collect(),
        chCounts.count(),
        chPcaHeader,
        chHeatmapHeader 
      )

    emit:
      chMergeCountsVersion         = mergeCounts.out.version
      chGenesatResults             = geneSaturation.out.genesatResults
      chGeneSaturationVersion      = geneSaturation.out.version
      chCountsPerGenetype          = getCountsPerGeneType.out.countsPerGenetype
      chGeneTypeVersion            = getCountsPerGeneType.out.version
      chExploratoryAnalysisResults = exploratoryAnalysis.out.exploratoryAnalysisResults 
      chAnaExpVersion              = exploratoryAnalysis.out.version 
      chCountsLogs
      chFeaturecountsVersion
      chHtseqVersion
}
