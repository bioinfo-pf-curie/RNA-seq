/*
 * Strandness - save value in Txt file
 */

process saveStrandness {
  label 'unix'
  label 'minCpu'
  label 'minMem'
  
  input:
  tuple val(meta), path(reads)

  output:
  path "*.txt", emit: savedStrandness

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  echo ${params.stranded} > ${prefix}_strandness.txt
  """
}

