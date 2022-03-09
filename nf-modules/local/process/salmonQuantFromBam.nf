/*
 * Salmon quant from BAM file
 */

process salmonQuantFromBam {
  tag "${meta.id}"
  label "salmon"
  label "medCpu"
  label "medMem"

  input:
  tuple val(meta), path(bam)
  path(transcriptsFasta)
  path(gtf)

  output:
  path("${prefix}"), emit: results
  path("versions.txt"), emit: versions

  script:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"
  def strandOpts = meta.singleEnd ? 'U' : 'IU'
  if (meta.strandness == 'forward') {
    strandOpts = meta.singleEnd ? 'SF' : 'ISF'
  } else if (meta.strandness == 'reverse') {
    strandOpts = meta.singleEnd ? 'SR' : 'ISR'
  }    
  """
  echo \$(salmon --version 2>&1) > versions.txt
  salmon quant \\
    --libType=$strandOpts \\
    -a ${bam} \\
    --threads ${task.cpus} \\
    -t ${transcriptsFasta} \\
    --geneMap ${gtf} \\
    ${args} \\
    -o ${prefix}
  """
}
