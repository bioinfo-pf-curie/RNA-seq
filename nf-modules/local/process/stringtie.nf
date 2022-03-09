/*
 * Stringtie - Reference-guided de novo isoform assembly
 * https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
 */

process stringtie {
  tag "${meta.id}"
  label 'stringtie'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bam), path(bai)
  path(gtf)

  output:
  tuple val(meta), path("*.coverage.gtf")   , emit: coverageGtf
  tuple val(meta), path("*.transcripts.gtf"), emit: transcriptGtf
  tuple val(meta), path("*.abundance.txt")  , emit: abundance
  tuple val(meta), path("*.ballgown")       , emit: ballgown
  path  "versions.txt"                      , emit: versions

  script:
  def args = task.ext.args ?: ''
  def strandOpts = ''
  if (meta.strandness == 'forward') {
      strandOpts = '--fr'
  } else if (meta.strandness == 'reverse') {
      strandOpts = '--rf'
  }
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  stringtie \\
    $bam \\
    $strandOpts \\
    -G $gtf \\
    -o ${prefix}.transcripts.gtf \\
    -A ${prefix}.gene.abundance.txt \\
    -C ${prefix}.coverage.gtf \\
    -b ${prefix}.ballgown \\
    -p $task.cpus \\
    ${args}
  echo "stringtie "\$(stringtie --version 2>&1) > versions.txt
  """
}
