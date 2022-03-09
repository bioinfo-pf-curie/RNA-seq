/* 
 * FeatureCounts - gene counting
 */

process featureCounts {
  tag "${meta.id}"
  label 'featurecounts'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bam), path(bai)
  path gtf

  output:
  tuple val(meta), path("*_counts.csv"), emit: counts
  path "*counts.csv.summary", emit: logs
  path "versions.txt", emit: versions

  script:
  def args   = task.ext.args ?: ''
  def featureCountsDirection = 0
  if (meta.strandness == 'forward'){
      featureCountsDirection = 1
  } else if ((meta.strandness == 'reverse')){
      featureCountsDirection = 2
  }
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  echo \$(featureCounts -v 2>&1 | sed '/^\$/d') > versions.txt
  featureCounts ${args} -T ${task.cpus} -a ${gtf} -o ${prefix}_counts.csv -p -s ${featureCountsDirection} ${bam}
  """
}

