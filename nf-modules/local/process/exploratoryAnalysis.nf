/*
 * Exploratory analysis
 */

// Stage config files
//pcaHeader = Channel.fromPath("$baseDir/assets/pcaHeader.txt")
//heatmapHeader = Channel.fromPath("$baseDir/assets/heatmapHeader.txt")

process exploratoryAnalysis {
  label 'r'
  label 'minCpu'
  label 'lowMem'
  publishDir "${params.outDir}/exploratoryAnalysis", mode: 'copy'

  when:
  !params.skipExpan && numSample > 1

  input:
  path tableRaw 
  path tableTpm
  val numSample
  path pcaHeader
  path heatmapHeader

  output:
  path "*.{txt,pdf,csv}", emit: exploratoryAnalysisResults
  path("v_R.txt")       , emit: version

  script:
  """
  R --version &> v_R.txt
  exploratory_analysis.r ${tableRaw}
  cat $pcaHeader deseq2_pca_coords_mqc.csv >> tmp_file
  mv tmp_file deseq2_pca_coords_mqc.csv
  cat $heatmapHeader vst_sample_cor_mqc.csv >> tmp_file
  mv tmp_file vst_sample_cor_mqc.csv
  """
}

