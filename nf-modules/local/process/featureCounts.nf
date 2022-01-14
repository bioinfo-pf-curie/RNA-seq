/* 
 * FeatureCounts - gene counting
 */

process featureCounts {
  tag "${prefix}"
  label 'featurecounts'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(prefix), path(bam), path(bai), val(strandness) // Channel [prefix, bam, bai, strandness]
  path gtf

  output:
  tuple val(prefix), path("${prefix}_counts.csv"), emit: counts
  path "${prefix}_counts.csv.summary", emit: logs
  path "versions.txt"                , emit: versions

  script:
  def args   = task.ext.args ?: ''
  def featureCountsDirection = 0
  if (strandness == 'forward'){
      featureCountsDirection = 1
  } else if ((strandness == 'reverse')){
      featureCountsDirection = 2
  }
  """
  echo \$(featureCounts -v 2>&1 | sed '/^\$/d') > versions.txt
  featureCounts ${args} -T ${task.cpus} -a ${gtf} -o ${prefix}_counts.csv -p -s ${featureCountsDirection} ${bam}
  """
}

