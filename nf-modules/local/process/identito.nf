/*
 * IDENTITO MONITORING
 */

process identito {
  label 'lowCpu'
  label 'medMem'
  label 'identito'

  when:
  !params.skipIdentito && !params.skipQC && params.polym

  input:
  path(fasta)
  path(fastaFai)
  path(polyms)
  path(bam)
  //tuple val(sampleId), val(sampleName), path("${sampleId}.md.bam"), path("${sampleId}.md.bam.bai")

  output:
  path("${prefix}_matrix.tsv"), emit: clustPolym
  path("v_bcftools.txt")      , emit: version 

  script:
  prefix = bam[0].toString() - ~/(_sorted)?(.markDups)?(.bam)?$/
  """
  bcftools --version &> v_bcftools.txt 2>&1 || true
  bcftools mpileup -R ${polyms} -f ${fasta} -x -A -B -q 20 -I -Q 0 -d 1000 --annotate FORMAT/DP,FORMAT/AD ${bam[0]} > ${prefix}_bcftools.vcf
  SnpSift extractFields -e "."  -s ";" ${prefix}_bcftools.vcf CHROM POS REF ALT GEN[*].DP GEN[*].AD > ${prefix}_bcftools.tsv
  apComputePolym.R ${prefix}_bcftools.tsv ${prefix}_matrix.tsv ${prefix} ${polyms}
  """
}

