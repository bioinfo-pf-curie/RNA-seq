/*
 * Strandness - save value in Txt file
 */

process saveStrandness {
  label 'unix'
  label 'minCpu'
  label 'minMem'
  
  input:
  tuple val(prefix), path(reads)

  output:
  path "*.txt", emit: savedStrandness

  script:
  """
  echo ${params.stranded} > ${prefix}_strandness.txt
  """
}

