/*
 * Reads distribution - gene type distribution
 */

process getCountsPerGeneType {
  label 'r'
  label 'minCpu'
  label 'lowMem'

  input:
  path counts
  path gtf

  output:
  path "*genetype.txt", emit: results
  path("versions.txt"), emit: versions

  script:
  """
  echo \$(R --version | awk 'NR==1{print \$1,\$3}') > versions.txt
  gene_type_expression.r ${counts} ${gtf} counts_genetype.txt
  """
}

