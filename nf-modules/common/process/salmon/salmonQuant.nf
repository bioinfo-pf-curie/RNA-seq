/*
 * Salmon quant from Fastq or BAM file
 */

process salmonQuant {
  tag "${meta.id}"
  label "salmon"
  label "medCpu"
  label "medMem"

  input:
  tuple val(meta), path(reads) // bam or fastq files
  path(index)
  path(transcriptsFasta)
  path(gtf)

  output:
  path("${prefix}"), emit: results
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"
  
  def reference = "--index $index"
  def inputReads = meta.singleEnd ? "-r $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
  if (reads.getExtension() == "bam"){
    inputReads = "-a $reads"
    reference = "-t ${transcriptsFasta}"
  }

  def strandOpts = meta.singleEnd ? 'U' : 'IU'
  if (meta.strandness == 'forward') {
    strandOpts = meta.singleEnd ? 'SF' : 'ISF'
  } else if (meta.strandness == 'reverse') {
    strandOpts = meta.singleEnd ? 'SR' : 'ISR'
  }    
  """
  echo \$(salmon --version 2>&1) > versions.txt
  salmon quant \\
    ${inputReads} \\
    ${reference} \\
    --libType=$strandOpts \\
    --threads ${task.cpus} \\
    --geneMap ${gtf} \\
    ${args} \\
    -o ${prefix}
  """
}
