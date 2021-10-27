
process rseqcPrep {
  tag "${prefix}"
  label 'bowtie2'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(prefix), path(reads)

  output:
  tuple val("${prefix}"), path("${prefix}_subsample.bam"), emit: bamRseqc
  path("versions.txt")                                   , emit: versions

  script:
  inputOpts = params.singleEnd ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
  """
  echo \$(bowtie2 --version | awk 'NR==1{print "bowtie2 "\$3}') > versions.txt
  bowtie2 --fast --end-to-end --reorder \\
          -p ${task.cpus} \\
          -u ${params.nCheck} \\
          -x ${params.bowtie2Index} \\
          ${inputOpts} > ${prefix}_subsample.bam 
   """
}