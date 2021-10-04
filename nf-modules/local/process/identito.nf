/*
 * IDENTITO MONITORING
 */

process identito {
  label 'lowCpu'
  label 'medMem'
  label 'identito'

  when:
  !params.skipIdentito && !params.skipQC && params.polyms

  input:
  path(fasta)
  path(fastaFai)
  path(polyms)
  tuple val(sampleId), val(sampleName), path("${sampleId}.md.bam"), path("${sampleId}.md.bam.bai")

  output:
  path("${sampleId}_matrix.tsv"), emit: clustPolym
  path("v_bcftools.txt")        , emit: version 

  script:
  """
  bcftools --version &> v_bcftools.txt 2>&1 || true
  bcftools mpileup -R ${polyms} -f ${fasta} -x -A -B -q 20 -I -Q 0 -d 1000 --annotate FORMAT/DP,FORMAT/AD ${sampleId}.md.bam > ${sampleId}_bcftools.vcf
  SnpSift extractFields -e "."  -s ";" ${sampleId}_bcftools.vcf CHROM POS REF ALT GEN[*].DP GEN[*].AD > ${sampleId}_bcftools.tsv
  apComputePolym.R ${sampleId}_bcftools.tsv ${sampleId}_matrix.tsv ${sampleId} ${polyms}
  """
}

