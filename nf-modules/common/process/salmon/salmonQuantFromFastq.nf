/*
 * Salmon quant from Fastq file
 */

process salmonQuantFromFastq {
  tag "${meta.id}"
  label "salmon"
  label "medCpu"
  label "extraMem"

  input:
  tuple val(meta), path(reads)
  path  index
  path  gtf

  output:
  path("${prefix}"), emit: results
  path("versions.txt"), emit: versions

  script:
  def args = task.ext.args ?: ''
  def strandOpts = meta.singleEnd ? 'U' : 'IU'
  if (meta.strandness == 'forward') {
    strandOpts = meta.singleEnd ? 'SF' : 'ISF'
  } else if (meta.strandness == 'reverse') {
    strandOpts = meta.singleEnd ? 'SR' : 'ISR'
  }
  def inputReads = meta.singleEnd ? "-r $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
  prefix = task.ext.prefix ?: "${meta.id}"
  """
  echo \$(salmon --version 2>&1) > versions.txt
  salmon quant \\
    --threads $task.cpus \\
    --libType=$strandOpts \\
    --validateMappings \\
    --geneMap ${gtf} \\
    ${args} \\
    $inputReads \\
    --index $index \\
    -o ${prefix}
  """
}
