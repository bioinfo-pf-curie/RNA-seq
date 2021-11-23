/*
 * Samtools - Flagstat
 */

process samtoolsFlagstat {
  tag "$prefix"
  label 'samtools'
  label 'minCpu'
  label 'lowMem'
  publishDir "${params.outDir}/mapping", mode: 'copy'
 
  input:
  tuple val(prefix), path (bam)

  output:
  tuple val(prefix), path("*flagstats"), emit: stats
  path("versions.txt") , emit: versions

  script:
  """
  echo \$(samtools --version | head -1) > versions.txt
  samtools flagstat ${bam} > ${prefix}.flagstats
  """
}
    