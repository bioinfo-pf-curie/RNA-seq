/*
 * Salmon quant from Fastq file
 * External parameters :
 * @ params.singleEnd : is the data single-end ?
 */

process salmonQuantFromFastq {
  tag "$prefix"
  label "salmon"
  label "medCpu"
  label "medMem"

  input:
  tuple val(prefix), path(reads), val(strandness)
  path  index
  path  gtf

  output:
  path("${prefix}"), emit: results
  path("versions.txt"), emit: versions

  script:
  def args = task.ext.args ?: ''
  strandOpts = params.singleEnd ? 'U' : 'IU'
  if (strandness == 'forward') {
    strandOpts = params.singleEnd ? 'SF' : 'ISF'
  } else if (strandness == 'reverse') {
    strandOpts = params.singleEnd ? 'SR' : 'ISR'
  }
  inputReads = params.singleEnd ? "-r $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
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
