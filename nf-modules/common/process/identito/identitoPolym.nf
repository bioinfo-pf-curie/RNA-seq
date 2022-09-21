/*
 * BQSR - detects systematic errors made by the sequencing machine
 */

process identitoPolym {
  tag "${meta.id}"
  label 'lowCpu'
  label 'medMem'
  label 'identito'

  input:
  tuple val(meta), path(bam), path(bai)
  path(fasta)
  path(fastaFai)
  path(polyms)

  output:
  path("*_matrix.tsv"), emit: polyms
  path("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  echo \$(bcftools --version | head -1) > versions.txt
  echo \$(SnpSift 2>&1| awk 'NR==1{print \$1,\$3}') >> versions.txt
  echo \$(R --version | awk 'NR==1{print \$1,\$3}') >> versions.txt

  bcftools mpileup -R ${polyms} -f ${fasta} -x -A -B -q 20 -I -Q 0 -d 1000 --annotate FORMAT/DP,FORMAT/AD ${bam[0]} > ${prefix}_bcftools.vcf
  SnpSift extractFields -e "."  -s ";" ${prefix}_bcftools.vcf CHROM POS REF ALT GEN[*].DP GEN[*].AD > ${prefix}_bcftools.tsv
  apComputePolym.R ${prefix}_bcftools.tsv ${prefix}_matrix.tsv ${prefix} ${polyms}
  """
}
