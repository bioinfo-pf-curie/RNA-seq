/*
 * HTSeqCounts - gene counting
 */

process htseqCounts {
  tag "${meta.id}"
  label 'htseq'
  label 'medCpu'
  label 'medMem'
  
  input:
  tuple val(meta), path(bam), path(bai), val(strandness)
  path gtf

  output: 
  tuple val(meta), path("*_counts.csv"), emit: counts
  path("*_counts.csv"), emit: logs
  path("versions.txt"), emit: versions 

  script:
  def args   = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def strandedOpt = '-s no' 
  if (strandness == 'forward'){
      strandedOpt= '-s yes'
  } else if ((strandness == 'reverse')){
      strandedOpt= '-s reverse'
  }
  """
  echo \$(htseq-count --version | awk '{print "HTSeq "\$1}') > versions.txt
  htseq-count ${args} ${strandedOpt} ${bam} ${gtf} > ${prefix}_counts.csv
  """
}
