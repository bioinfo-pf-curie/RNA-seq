process dupradar {
  tag "${bamMd[0]}"
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

  when:
  !params.skipQC && !params.skipDupradar

  input:
  path bamMd
  path gtf
  val parseRes

  output:
  path "*.{pdf,txt}", emit : dupradarResults

  script: 
  def dupradarDirection = 0
  if (parseRes == 'forward'){
      dupradarDirection = 1
  } else if ((parseRes == 'reverse')){
      dupradarDirection = 2
  }
  def paired = params.singleEnd ? 'single' :  'paired'
  """
  dupRadar.r ${bamMd[0]} ${gtf} ${dupradarDirection} ${paired} ${task.cpus}
  """
}
