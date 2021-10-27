/*
 * Counts
 */

process mergeCounts {
  publishDir "${params.outDir}/counts", mode: 'copy'
  label 'r'
  label 'minCpu'
  label 'medMem'

  input:
  path inputCounts
  path gtf
  val strandness
  val tool

  output:
  path('tablecounts_raw.csv'), emit: countsTable
  path('tablecounts_tpm.csv'), emit: tpmTable
  path('tableannot.csv')     , emit: genesAnnot
  path('versions.txt')       , emit: versions

  script:
  """
  echo \$(R --version | awk 'NR==1{print \$1,\$3}') > versions.txt
  echo -e ${inputCounts} | tr " " "\n" > listofcounts.tsv
  echo -n "${strandness}" | sed -e "s/\\[//" -e "s/\\]//" -e "s/,//g" | tr " " "\n" > listofstrandness.tsv
  makeCountTable.r listofcounts.tsv ${gtf} ${tool} listofstrandness.tsv
  """
}
