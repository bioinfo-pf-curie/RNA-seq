/*
 * HTSeqCounts
 */

process htseqCounts {
  tag "${prefix}"
  label 'htseq'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outDir}/counts", mode: 'copy',
    saveAs: {filename ->
      if (filename.indexOf("_gene.HTseqCounts.txt.summary") > 0) "gene_count_summaries/$filename"
      else if (filename.indexOf("_gene.HTseqCounts.txt") > 0) "gene_counts/$filename"
      else "$filename"
    }
  
  input:
  tuple val(prefix), path(bam), path(bai), val(strandness)
  path gtf

  output: 
  tuple val(prefix), path("${prefix}_counts.csv"), emit: counts
  path("${prefix}_counts.csv"), emit: logs
  path("versions.txt"), emit: versions 

  script:
  def strandedOpt = '-s no' 
  if (strandness == 'forward'){
      strandedOpt= '-s yes'
  } else if ((strandness == 'reverse')){
      strandedOpt= '-s reverse'
  }
  """
  echo \$(htseq-count --version | awk '{print "HTSeq "\$1}') > versions.txt
  htseq-count ${params.htseqOpts} ${strandedOpt} ${bam} $gtf > ${prefix}_counts.csv
  """
}
