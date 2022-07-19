/*
 * Exploratory analysis from gene counts table
 */

process exploratoryAnalysis {
  label 'r'
  label 'minCpu'
  label 'lowMem'

  input:
  path counts
  path pcaHeader
  path heatmapHeader

  output:
  path "*.{pdf,csv}", optional:true, emit: results
  path("versions.txt"), emit: versions

  script:
  """
  echo \$(R --version | awk 'NR==1{print \$1,\$3}') > versions.txt
  numSample=\$(awk -F, 'NR==2{print NF-1}' ${counts})
  if [[ \$numSample > 1 ]]; then
    exploratory_analysis.r ${counts}
    cat $pcaHeader deseq2_pca_coords_mqc.csv >> tmp_file
    mv tmp_file deseq2_pca_coords_mqc.csv
    cat $heatmapHeader vst_sample_cor_mqc.csv >> tmp_file
    mv tmp_file vst_sample_cor_mqc.csv
  fi
  """
}

