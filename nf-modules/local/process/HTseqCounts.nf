/*
 * Counts
 */

process HTseqCounts {
  tag "${bam}"
  label 'htseq'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outDir}/counts", mode: 'copy',
    saveAs: {filename ->
      if (filename.indexOf("_gene.HTseqCounts.txt.summary") > 0) "gene_count_summaries/$filename"
      else if (filename.indexOf("_gene.HTseqCounts.txt") > 0) "gene_counts/$filename"
      else "$filename"
    }
  
  when:
  params.counts == 'HTseqCounts'

  input:
  path bam
  path gtf
  val parseRes

  output: 
  path "*_counts.csv", emit: counts
  path("v_htseq.txt"), emit: version 

  script:
  def strandedOpt = '-s no' 
  if (parseRes == 'forward'){
      strandedOpt= '-s yes'
  } else if ((parseRes == 'reverse')){
      strandedOpt= '-s reverse'
  }
  """
  htseq-count -h | grep version  &> v_htseq.txt
  htseq-count ${params.htseqOpts} ${strandedOpt} ${bam[0]} $gtf > ${bam[0].baseName}_counts.csv
  """
}
