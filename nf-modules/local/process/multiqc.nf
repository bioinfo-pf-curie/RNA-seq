process multiqc {
  label 'multiqc'
  label 'minCpu'
  label 'lowMem'
  publishDir "${params.outDir}/MultiQC", mode: 'copy'

  when:
  !params.skipMultiQC

  input:
  path splan
  path metadata
  path multiqcConfig
  path (fastqc:'fastqc/*')
  path ('rrna/*')
  path ('alignment/*')
  path ('strandness/*')
  path ('qualimap/*')
  path ('preseq/*')
  path ('genesat/*')
  path ('dupradar/*')
  path ('picard/*')
  path ('counts/*')
  path ('genetype/*')
  path ('identito/*')
  path ('exploratoryAnalysis_results/*')
  path ('softwareVersions/*')
  path ('workflowSummary/*')
  path ('workflowSummary/*')

  output:
  path splan
  path "*report.html"
  path "*_data"

  script:
  rtitle = customRunName ? "--title \"$customRunName\"" : ''
  rfilename = customRunName ? "--filename " + customRunName + "_rnaseq_report" : "--filename rnaseq_report"
  metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
  splanOpts = params.samplePlan ? "--splan ${params.samplePlan}" : ""
  isPE = params.singleEnd ? 0 : 1
    
  modulesList = "-m custom_content -m preseq -m rseqc -m bowtie1 -m hisat2 -m star -m cutadapt -m fastqc -m qualimap"
  modulesList = params.counts == 'featureCounts' ? "${modulesList} -m featureCounts" : "${modulesList}"  
  modulesList = params.counts == 'HTseqCounts' ? "${modulesList} -m htseq" : "${modulesList}"  
 
  warn=skippedPoorAlignment.size() > 0 ? "--warn workflowSummary/warnings.txt" : ""
  """
  stats2multiqc.sh ${splan} ${params.aligner} ${isPE}
  medianReadNb="\$(sort -t, -k3,3n mq.stats | awk -F, '{a[i++]=\$3;} END{x=int((i+1)/2); if (x<(i+1)/2) printf "%.0f", (a[x-1]+a[x])/2; else printf "%.0f",a[x-1];}')"
  mqc_header.py --name "RNA-seq" --version ${workflow.manifest.version} ${metadataOpts} ${splanOpts} --nbreads \${medianReadNb} ${warn} > multiqc-config-header.yaml
  multiqc . -f $rtitle $rfilename -c $multiqcConfig -c multiqc-config-header.yaml $modulesList
  """    
}
