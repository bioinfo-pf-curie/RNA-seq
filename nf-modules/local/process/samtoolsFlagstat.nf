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
  path("v_samtools.txt") , emit: version

  script:
  """
  samtools --version &> v_samtools.txt
  samtools flagstat ${bam} > ${prefix}.flagstats
  """
}
    