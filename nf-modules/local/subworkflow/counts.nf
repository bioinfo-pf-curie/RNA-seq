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

workflow countFlow {
    // required inputs
    take:
      chBam
      chGtf
      chStrandedResults
      chStarCounts
      chStarLogCounts

    // workflow implementation
    main:
      // Counts
      featureCounts(
        chBam,
        chGtf.collect(),
        chStrandedResults
      )

      HTseqCounts (
        chBam,
        chGtf.collect(),
        chStrandedResults
      )

      chCounts = Channel.empty()
     if( params.counts == 'featureCounts' ){
       chCounts = featureCounts.out.counts
     } else if (params.counts == 'HTseqCounts'){
       chCounts = HTseqCounts.out.counts
     }else if (params.counts == 'star'){
       chCounts = chStarCounts
     }

      mergeCounts(
        chCounts.collect(),
        chGtf.collect(),
        chStrandedResults
      )

      chCountsLogs = Channel.empty()
      if( params.counts == 'featureCounts' ){
        chCountsLogs = featureCounts.out.logs
      } else if (params.counts == 'HTseqCounts'){
        chCountsLogs = HTseqCounts.out.counts
      }else if (params.counts == 'star'){
        chCountsLogs = chStarLogCounts
      }

      geneSaturation(

      )

      getCountsPerGeneType(

      )

      exploratoryAnalysis(

      )

    emit:
      chFeaturecountsVersion = featureCounts.out.version
      chHtseqCounts          = HTseqCounts.out.htseqCounts
      chHtseqVersion         = HTseqCounts.out.version
      chMergeCountsVersion   = mergeCounts.out.version
      chCountsLogs

}
