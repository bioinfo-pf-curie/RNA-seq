// Counts

process featureCounts {
  tag "${bam.baseName}"
  label 'featurecounts'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outDir}/counts", mode: 'copy',
    saveAs: {filename ->
      if (filename.indexOf("_counts.csv.summary") > 0) "gene_count_summaries/$filename"
      else if (filename.indexOf("_counts.csv") > 0) "gene_counts/$filename"
      else "$filename"
   }

  when:
  params.counts == 'featureCounts'

  input:
  path bam
  path gtf
  val parseRes

  output:
  path "${bam.baseName}_counts.csv"        , emit: counts
  path "${bam.baseName}_counts.csv.summary", emit: logs
  path "v_featurecounts.txt"               , emit: version

  script:
  def featureCountsDirection = 0
  if (parseRes == 'forward'){
      featureCountsDirection = 1
  } else if ((parseRes == 'reverse')){
      featureCountsDirection = 2
  }
  """
  featureCounts -v &> v_featurecounts.txt
  featureCounts ${params.featurecountsOpts} -T ${task.cpus} -a ${gtf} -o ${bam.baseName}_counts.csv -p -s ${featureCountsDirection} ${bam}
  """
}

