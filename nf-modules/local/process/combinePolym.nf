/* 
 * Identito - polym
 */

process combinePolym {
  label 'minCpu'
  label 'lowMem'
  label 'identito'

  publishDir "${params.outDir}/identito", mode: 'copy'

  when:
  !params.skipQC && !params.skipIdentito

  input:
  path(matrix)

  output:
  path("*.csv")  , emit: clustPolymResults
  path("v_R.txt"), emit: version

  script:
  """
  R --version &> v_R.txt
  (head -1 "${matrix[0]}"; tail -n +2 -q *matrix.tsv) > clust_mat.tsv
  computeClust.R clust_mat.tsv ./ 10
  """
}