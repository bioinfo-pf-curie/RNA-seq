/*
 * RSeQC - Quick mapping on a subset of reads for RSeQC
 * External parameters :
 * @ params.singleEnd :	is data	single-end sequencing ?
 * @ params.nCheck : number of reads to subset
 */

process rseqcPrep {
  tag "${prefix}"
  label 'bowtie2'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(prefix), path(reads)
  path(index)

  output:
  tuple val("${prefix}"), path("${prefix}_subsample.bam"), emit: bamRseqc
  path("versions.txt")                                   , emit: versions

  script:
  inputOpts = params.singleEnd ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
  subset = params.nCheck ?: 200000
  """
  localIndex=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
  echo \$(bowtie2 --version | awk 'NR==1{print "bowtie2 "\$3}') > versions.txt
  bowtie2 --fast --end-to-end --reorder \\
          -p ${task.cpus} \\
          -u ${subset} \\
          -x \${localIndex} \\
          ${inputOpts} > ${prefix}_subsample.bam 
   """
}