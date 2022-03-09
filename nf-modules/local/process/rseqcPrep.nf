/*
 * RSeQC - Quick mapping on a subset of reads for RSeQC
 */

process rseqcPrep {
  tag "${meta.id}"
  label 'bowtie2'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(reads)
  path(index)

  output:
  tuple val(meta), path("*_subsample.bam"), emit: bamRseqc
  path("versions.txt"), emit: versions

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def inputOpts = meta.singleEnd ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
  """
  localIndex=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
  echo \$(bowtie2 --version | awk 'NR==1{print "bowtie2 "\$3}') > versions.txt
  bowtie2 --fast --end-to-end --reorder \\
          -p ${task.cpus} \\
          ${args} \\
          -x \${localIndex} \\
          ${inputOpts} > ${prefix}_subsample.bam 
   """
}