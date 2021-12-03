/*
 * MultiQC for RNA-seq report
 * External parameters :
 * @ params.singleEnd :	is data	single-end sequencing ?
 */

process multiqc {
  label 'multiqc'
  label 'minCpu'
  label 'lowMem'
  publishDir "${params.outDir}/MultiQC", mode: 'copy'

  when:
  !params.skipMultiQC

  input:
  val customRunName
  path splan
  path metadata
  path multiqcConfig
  path ('fastqc/*')
  path ('rrna/*')
  path ('alignment/*')
  path ('strandness/*')
  path ('qualimap/*')
  path ('preseq/*')
  path ('identito/*')
  path ('picard/*')
  path ('dupradar/*')
  path ('counts/*')
  path ('genesat/*')
  path ('genetype/*')
  path ('exploratoryAnalysis/*')
  path ('gffCompare/*')
  path ('softwareVersions/*')
  path ('workflowSummary/*')
  path warnings
  //val skippedPoorAlignment

  output:
  path splan, emit: splan
  path "*report.html", emit: report
  path "*_data", emit: data

  script:
  rtitle = customRunName ? "--title \"$customRunName\"" : ''
  rfilename = customRunName ? "--filename " + customRunName + "_rnaseq_report" : "--filename rnaseq_report"
  alignerOpts = params.aligner ?: params.pseudoAligner ?: ''
  metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
  splanOpts = params.samplePlan ? "--splan ${params.samplePlan}" : ""
  isPE = params.singleEnd ? 0 : 1
    
  modulesList = "-m custom_content -m preseq -m rseqc -m bowtie1 -m hisat2 -m star -m cutadapt -m fastqc -m qualimap -m salmon -m gffcompare"
  modulesList = params.counts == 'featureCounts' ? "${modulesList} -m featureCounts" : "${modulesList}"  
  modulesList = params.counts == 'HTseqCounts' ? "${modulesList} -m htseq" : "${modulesList}"  
 
  //warn=skippedPoorAlignment.size() > 0 ? "--warn warnings.txt" : ""
  warn = warnings.name == 'warnings.txt' ? "--warn warnings.txt" : ""
  """
  stats2multiqc.sh ${splan} ${alignerOpts} ${isPE}
  medianReadNb="\$(sort -t, -k3,3n mq.stats | awk -F, '{a[i++]=\$3;} END{x=int((i+1)/2); if (x<(i+1)/2) printf "%.0f", (a[x-1]+a[x])/2; else printf "%.0f",a[x-1];}')"
  mqc_header.py --name "RNA-seq" --version ${workflow.manifest.version} ${metadataOpts} ${splanOpts} --nbreads \${medianReadNb} ${warn} > multiqc-config-header.yaml
  multiqc . -f $rtitle $rfilename -c $multiqcConfig -c multiqc-config-header.yaml $modulesList
  """    
}
