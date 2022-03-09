/*
 * Samtools - Index
 */

process samtoolsIndex {
  tag "${meta.id}"
  label 'samtools'
  label 'minCpu'
  label 'lowMem'
 
  input:
  tuple val(meta), path (bam)

  output:
  tuple val(meta), path("*bam.bai"), emit: bai
  path("versions.txt") , emit: versions

  script:
  """
  echo \$(samtools --version | head -1) > versions.txt
  samtools index ${bam}
  """
}
    