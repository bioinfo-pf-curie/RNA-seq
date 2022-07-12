/*
 * Identito Monitoring - combine per sample polym information
 */

process identitoCombine {
  label 'lowCpu'
  label 'lowMem'
  label 'identito'

  input:
  path(matrix)

  output:
  path("*.{tsv,csv,png}"), emit: results
  path("versions.txt"), emit: versions

  script:
  """
  echo \$(R --version | awk 'NR==1{print \$1,\$3}') > versions.txt
  (head -1 "${matrix[0]}"; tail -n +2 -q *matrix.tsv) > identito_polym.tsv
  """
}
