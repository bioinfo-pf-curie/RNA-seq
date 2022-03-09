/*
 * DupRadar for duplicates assessment
 */

process dupradar {
  tag "${meta.id}"
  label 'dupradar'
  label 'minCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(bam)
  path gtf

  output:
  path "${bam.baseName}*.{pdf,txt}", emit : results
  path("versions.txt"), emit: versions

  script: 
  def dupradarDirection = 0
  if (meta.strandness == 'forward'){
      dupradarDirection = 1
  } else if ((meta.strandness == 'reverse')){
      dupradarDirection = 2
  }
  def paired = meta.singleEnd ? 'single' : 'paired'
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  echo \$(R --version | awk 'NR==1{print \$1,\$3}') > versions.txt
  dupRadar.r ${bam} ${gtf} ${dupradarDirection} ${paired} ${task.cpus}
  """
}
