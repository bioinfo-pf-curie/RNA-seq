/*
 * FastQC process
 */

process fastqc {
  tag "${prefix}"
  label 'fastqc'
  label 'medCpu'
  label 'lowMem'
  publishDir "${params.outDir}/fastqc", mode: 'copy',
    saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

  input:
  tuple val(prefix), path(reads)

  output:
  path("*_fastqc.{zip,html}"), emit: results
  path("versions.txt")       , emit: versions

  script:
  if (params.singleEnd){
    pbase = reads[0].toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
    """
    echo \$(fastqc --version) > versions.txt
    fastqc -q $reads --threads ${task.cpus}
    mv ${pbase}_fastqc.html ${prefix}_fastqc.html
    mv ${pbase}_fastqc.zip ${prefix}_fastqc.zip
    """
  }else{
    pbase_1 = reads[0].toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
    pbase_2 = reads[1].toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
    """
    echo \$(fastqc --version) > versions.txt
    fastqc -q $reads --threads ${task.cpus}
    mv ${pbase_1}_fastqc.html ${prefix}_R1_fastqc.html
    mv ${pbase_1}_fastqc.zip ${prefix}_R1_fastqc.zip
    mv ${pbase_2}_fastqc.html ${prefix}_R2_fastqc.html
    mv ${pbase_2}_fastqc.zip ${prefix}_R2_fastqc.zip
    """
  }
}
