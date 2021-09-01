/*
 * Strandness
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
  
  when:
  params.stranded == 'reverse' || params.stranded == 'no' || params.stranded == 'yes' || (params.stranded == 'auto' && !params.bed12)

  input:
  tuple val(prefix), path(reads)

  output:
  file "*.txt", emit: savedStrandness

  script:
  """
  echo ${params.stranded} > ${prefix}_strandness.txt
  """
}
