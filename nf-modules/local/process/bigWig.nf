/*
 * Generate bigwig file
 */

process bigWig {
  tag "${meta.id}"
  label 'deeptools'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bam), path(bai)

  output:
  path('*.bigwig') , emit: bigWig
  path("versions.txt"), emit: versions

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  strandOpts = meta.strandness == 'forward' ? '--filterRNAstrand forward' : meta.strandness == 'reverse' ? '--filterRNAstrand reverse' : ''
  """
  echo \$(bamCoverage --version) > versions.txt
  bamCoverage -b ${bam} \\
              -o ${prefix}.bigwig \\
              -p ${task.cpus} \\
              ${strandOpts} \\
	      ${args}
  """
}