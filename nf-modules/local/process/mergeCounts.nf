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
  val parseRes

  output:
  path 'tablecounts_raw.csv', emit: counts
  path 'tablecounts_tpm.csv', emit: tpmCounts
  path 'tableannot.csv'     , emit: genesAnnot
  path("v_R.txt")           , emit: version

  script:
  """
  R --version &> v_R.txt
  echo ${inputCounts} | tr " " "\\n" > listofcounts.tsv
  echo ${parseRes} | sed -e "s/\\[//" -e "s/\\]//" -e "s/,//g" | tr " " "\\n" > listofstrandness.tsv
  makeCountTable.r listofcounts.tsv ${gtf} ${params.counts} listofstrandness.tsv
  """
}
