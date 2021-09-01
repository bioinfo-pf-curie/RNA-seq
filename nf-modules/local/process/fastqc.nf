process fastqc {
  tag "${prefix}"
  label 'fastqc'
  label 'medCpu'
  label 'lowMem'
  publishDir "${params.outDir}/fastqc", mode: 'copy',
    saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

  when:
  !params.skipQC && !params.skipFastqc

  input:
  tuple val(prefix), path(reads)

  output:
  path("*_fastqc.{zip,html}"), emit: mqc 
  path("v_fastqc.txt")       , emit: version 

  script:
  pbase = reads[0].toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
  """
  fastqc --version &> v_fastqc.txt
  fastqc -q $reads
  mv ${pbase}_fastqc.html ${prefix}_fastqc.html
  mv ${pbase}_fastqc.zip ${prefix}_fastqc.zip
  """
}
