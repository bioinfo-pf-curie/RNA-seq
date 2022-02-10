/*
 * DupRadar for duplicates assessment
 * External parameters :
 * @ params.singleEnd : is data single-end sequencing ?
 */

process dupradar {
  tag "${prefix}"
  label 'dupradar'
  label 'minCpu'
  label 'lowMem'

  input:
  tuple val(prefix), path(bam), val(strandness)
  path gtf

  output:
  path "${prefix}*.{pdf,txt}", emit : results
  path("versions.txt"), emit: versions

  script: 
  def dupradarDirection = 0
  if (strandness == 'forward'){
      dupradarDirection = 1
  } else if ((strandness == 'reverse')){
      dupradarDirection = 2
  }
  def paired = params.singleEnd ? 'single' : 'paired'
  """
  echo \$(R --version | awk 'NR==1{print \$1,\$3}') > versions.txt
  dupRadar.r ${bam} ${gtf} ${dupradarDirection} ${paired} ${task.cpus}
  """
}
