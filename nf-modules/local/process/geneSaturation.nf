/*
 * Gene-based saturation
 */

process geneSaturation {
  label 'r'
  label 'minCpu'
  label 'medMem'
  publishDir "${params.outDir}/geneSaturation" , mode: 'copy'

  when:
  !params.skipQC && !params.skipSaturation

  input:
  path inputCounts

  output:
  path "*gcurve.txt", emit: genesatResults
  path("v_R.txt")   , emit: version

  script:
  """
  R --version &> v_R.txt
  gene_saturation.r $inputCounts counts.gcurve.txt
  """
}

