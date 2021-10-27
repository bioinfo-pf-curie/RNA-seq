/*
 * Exploratory analysis
 */

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
  path "*.{txt,pdf,csv}", emit: results
  path("versions.txt"), emit: versions

  script:
  """
  echo \$(R --version | awk 'NR==1{print \$1,\$3}') > versions.txt
  exploratory_analysis.r ${tableRaw}
  cat $pcaHeader deseq2_pca_coords_mqc.csv >> tmp_file
  mv tmp_file deseq2_pca_coords_mqc.csv
  cat $heatmapHeader vst_sample_cor_mqc.csv >> tmp_file
  mv tmp_file vst_sample_cor_mqc.csv
  """
}

