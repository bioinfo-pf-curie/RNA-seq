/* 
 * Identito - polym
 */

process getPolym {
  label 'minCpu'
  label 'medMem'
  label 'identito'

  publishDir "${params.outDir}/identito", mode: 'copy'

  when:
  !params.skipQC && !params.skipIdentito

  input:
  path(fasta)
  path(polyms)
  path(bam)

  output:
  path("v_bcftools.txt"), emit: version
  file("*matrix.tsv")   , emit: clustPolym

  script:
  """
  bcftools --version &> v_bcftools.txt 2>&1 || true
  bcftools mpileup -R ${polyms} -f ${fasta} -x -A -B -q 20 -I -Q 0 -d 1000 --annotate FORMAT/DP,FORMAT/AD ${bam[0]} > ${bam[0].baseName}_bcftools.vcf
  SnpSift extractFields -e "."  -s ";" ${bam[0].baseName}_bcftools.vcf CHROM POS REF ALT GEN[*].DP GEN[*].AD > ${bam[0].baseName}_bcftools.tsv
  computePolym.R ${bam[0].baseName}_bcftools.tsv ${bam[0].baseName}_matrix.tsv ${bam[0].baseName} ${polyms}
  """
}
