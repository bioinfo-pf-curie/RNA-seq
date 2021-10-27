/* 
 * FeatureCounts
 */

process featureCounts {
  tag "${prefix}"
  label 'featurecounts'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outDir}/counts", mode: 'copy',
    saveAs: {filename ->
      if (filename.indexOf("_counts.csv.summary") > 0) "gene_count_summaries/$filename"
      else if (filename.indexOf("_counts.csv") > 0) "gene_counts/$filename"
      else "$filename"
   }

  input:
  tuple val(prefix), path(bam), path(bai) // Channel [prefix, bam, bai]
  path gtf
  val strandness

  output:
  path "${prefix}_counts.csv"        , emit: counts
  path "${prefix}_counts.csv.summary", emit: logs
  path "versions.txt"                , emit: versions

  script:
  def featureCountsDirection = 0
  if (strandness == 'forward'){
      featureCountsDirection = 1
  } else if ((strandness == 'reverse')){
      featureCountsDirection = 2
  }
  """
  echo \$(featureCounts -v 2>&1 | sed '/^\$/d') > versions.txt
  featureCounts ${params.featurecountsOpts} -T ${task.cpus} -a ${gtf} -o ${prefix}_counts.csv -p -s ${featureCountsDirection} ${bam}
  """
}

