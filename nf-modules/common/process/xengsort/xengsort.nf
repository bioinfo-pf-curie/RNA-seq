/*
 * Xengsort - deconvolute Human/Mouse reads from PDX samples
 */

process xengsort {
  tag "${meta.id}"
  label 'medCpu'
  label 'highMem'
  label 'xengsort'

  input:
  tuple val(meta),path(reads)
  path (index)

  output :
  tuple val(meta),path("*graft*.fastq.gz"), emit: fastqHuman
  tuple val(meta),path("*host*.fastq.gz"), emit: fastqMouse
  path("*.log"), emit: logs
  path("versions.txt") , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def inputs = meta.singleEnd ? "--fastq <(zcat ${reads})" : "--fastq <(zcat ${reads[0]})  --pairs <(zcat ${reads[1]})"
  """
  echo "xengsort "\$(xengsort --version) > versions.txt
  xengsort classify -T ${task.cpus} --index ${index} ${inputs} --prefix ${prefix} ${args} > ${prefix}_xengsort.log
  for f in *.fq; do out=\$(echo "\$f" | sed -e 's/.1.fq/_R1.fastq/g' -e 's/.2.fq/_R2.fastq/' -e 's/.fq/.fastq/'); mv -- "\$f" "\$out"; done
  gzip *.fastq
  """
}
