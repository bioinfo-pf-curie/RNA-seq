/*
 * rRNA mapping using bowtie
 */

process rRNAMapping {
  tag "${meta.id}"
  label 'bowtie'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(reads)
  path(index)

  output:
  tuple val(meta), path("*fastq.gz"), emit: filteredReads
  path "*.log"                        , emit: logs
  path("versions.txt")                , emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  inputOpts = meta.singleEnd ? "${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
  """
  localIndex=`find -L ./ -name "*.rev.1.ebwt" | sed 's/.rev.1.ebwt//'`
  echo \$(bowtie --version | awk 'NR==1{print "bowtie "\$3}') > versions.txt
  bowtie ${args} \\
         -p ${task.cpus} \\
         --un ${prefix}_norRNA.fastq \\
         --sam \${localIndex} \\
         ${inputOpts} \\
         ${prefix}.sam  2> ${prefix}.log

  ## fix bug in MultiQC
  sed -i -e s'/one alignment/one reported alignment/' ${prefix}.log

  gzip -f ${prefix}_norRNA*.fastq 
  rm -f ${prefix}.sam
  """
}


