/*
 * Strandness - save value in Txt file
 */

process saveStrandness {
  label 'unix'
  label 'minCpu'
  label 'minMem'
  publishDir "${params.outDir}/strandness" , mode: 'copy',
    saveAs: {filename ->
      if (filename.indexOf(".txt") > 0) "$filename"
      else null
    }
  
  input:
  tuple val(prefix), path(reads)

  output:
  path "*.txt", emit: savedStrandness

  script:
  """
  echo ${params.stranded} > ${prefix}_strandness.txt
  """
}

