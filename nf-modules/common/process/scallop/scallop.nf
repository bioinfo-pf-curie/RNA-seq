/*
 * Scallop - Reference-guided de novo isoform assembly
 * https://github.com/Kingsford-Group/scallop
 */

process scallop {
  tag "${meta.id}"
  label 'scallop'
  label 'lowCpu'
  label 'highMem'

  input:
  tuple val(meta), path(bam), path(bai)

  output:
  tuple val(meta), path("*Transcripts.gtf"), emit: transcriptGtf
  path  "versions.txt"                     , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def strandOpts = '--library_type unstranded'
  if (meta.strandness == 'forward') {
    strandOpts = '--library_type second'
  } else if (meta.strandness == 'reverse') {
    strandOpts = '--library_type first'
  }
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  scallop \\
    -i ${bam} \\
    ${strandOpts} \\
    ${args} \\
    -o ${prefix}.scallopTranscripts.gtf
  echo "scallop "\$(scallop --version 2>&1) > versions.txt
  """
}
