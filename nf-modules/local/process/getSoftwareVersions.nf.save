process getSoftwareVersions{
  label 'python'
  label 'minCpu'
  label 'lowMem'
  publishDir path: "${params.outDir}/softwareVersions", mode: "copy"

  when:
  !params.skipSoftVersions

  input:
  path 'v_fastqc.txt' 
  path 'v_star.txt'
  path 'v_hisat2.txt' 
  path 'v_bowtie.txt' 
  path 'v_bowtie2.txt'
  path 'v_samtools.txt' 
  path 'v_picard.txt' 
  path 'v_preseq.txt' 
  path 'v_R.txt'
  path 'v_rseqc.txt' 
  path 'v_featurecounts.txt' 
  path 'v_deeptools.txt' 
  path 'v_bcftools.txt' 
  path 'v_htseq.txt' 
  path 'v_qualimap.txt' 

  output:
  path 'software_versions_mqc.yaml', emit: chSoftwareVersionsYaml

  script:
  """
  echo $workflow.manifest.version &> v_pipeline.txt
  echo $workflow.nextflow.version &> v_nextflow.txt
  scrape_software_versions.py &> software_versions_mqc.yaml
  """
}
