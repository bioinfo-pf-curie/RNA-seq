/*
 * Reads distribution
 */

process getCountsPerGeneType {
  label 'r'
  label 'minCpu'
  label 'lowMem'
  publishDir "${params.outDir}/readDistribution", mode: 'copy'

  when:
  !params.skipReaddist

  input:
  path tpmGenetype
  path gtf

  output:
  path "*genetype.txt", emit: countsPerGenetype
  path("v_R.txt")     , emit: version

  script:
  """
  R --version &> v_R.txt
  gene_type_expression.r ${tpmGenetype} ${gtf} counts_genetype.txt
  """
}

