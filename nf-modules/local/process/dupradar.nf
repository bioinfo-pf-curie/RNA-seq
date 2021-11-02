process dupradar {
  tag "${prefix}"
  label 'dupradar'
  label 'minCpu'
  label 'lowMem'
  publishDir "${params.outDir}/dupradar", mode: 'copy',
    saveAs: {filename ->
      if (filename.indexOf("_duprateExpDens.pdf") > 0) "scatter_plots/$filename"
      else if (filename.indexOf("_duprateExpBoxplot.pdf") > 0) "box_plots/$filename"
      else if (filename.indexOf("_expressionHist.pdf") > 0) "histograms/$filename"
      else if (filename.indexOf("_dupMatrix.txt") > 0) "gene_data/$filename"
      else if (filename.indexOf("_duprateExpDensCurve.txt") > 0) "scatter_curve_data/$filename"
      else if (filename.indexOf("_intercept_slope.txt") > 0) "intercepts_slopes/$filename"
      else "$filename"
    }

  input:
  tuple val(prefix), path(bam), val(strandness)
  path gtf

  output:
  path "${prefix}*.{pdf,txt}", emit : results
  path("versions.txt"), emit: versions

  script: 
  def dupradarDirection = 0
  if (strandness == 'forward'){
      dupradarDirection = 1
  } else if ((strandness == 'reverse')){
      dupradarDirection = 2
  }
  def paired = params.singleEnd ? 'single' : 'paired'
  """
  echo \$(R --version | awk 'NR==1{print \$1,\$3}') > versions.txt
  dupRadar.r ${bam} ${gtf} ${dupradarDirection} ${paired} ${task.cpus}
  """
}
