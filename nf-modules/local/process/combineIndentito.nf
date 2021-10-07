/*
 * IDENTITO MONITORING
 */
process combineIndentito {
  label 'lowCpu'
  label 'lowMem'
  label 'identito'

  publishDir "${params.outDir}/preprocessing/metrics/identito", mode: 'copy'

  when:
  !params.skipIdentito && !params.skipQC

  input:
  path(matrix)
  path(splan)

  output:
  path("*.{tsv,csv,png}"), emit: clustPolymResults
  path("*.png")          , optional:true, emit: clustPolymPlot

  script:
  """
  (head -1 "${matrix[0]}"; tail -n +2 -q *matrix.tsv) > identito_polym.tsv
  apComputeClust.R --input identito_polym.tsv --dist ejaccard
  """
}

