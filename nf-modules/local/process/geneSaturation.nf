/*
 * Gene-based saturation
 */

process geneSaturation {
  label 'r'
  label 'minCpu'
  label 'medMem'
  publishDir "${params.outDir}/geneSaturation" , mode: 'copy'

  input:
  path inputCounts

  output:
  path "*gcurve.txt", emit: results
  path("versions.txt"), emit: versions

  script:
  """
  echo \$(R --version | awk 'NR==1{print \$1,\$3}') > versions.txt
  gene_saturation.r $inputCounts counts.gcurve.txt
  """
}

