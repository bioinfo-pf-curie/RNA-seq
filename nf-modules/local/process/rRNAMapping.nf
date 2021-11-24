/*
 * rRNA mapping using bowtie
 * External parameters :
 * @ params.singleEnd :	is data	single-end sequencing ?
 * @ params.bowtieOpts : additional bowtie parameters
 */

process rRNAMapping {
  tag "${prefix}"
  label 'bowtie'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outDir}/rRNAmapping", mode: 'copy',
    saveAs: {filename ->
      if (filename.indexOf("fastq.gz") > 0 &&  params.saveAlignedIntermediates) filename
      else if (filename.indexOf(".log") > 0) "logs/$filename"
      else null
    }

  input:
  tuple val(prefix), path(reads)
  path(index)

  output:
  tuple val(prefix), path("*fastq.gz"), emit: filteredReads
  path "*.log"                        , emit: logs
  path("versions.txt")                , emit: versions

  script:
  inputOpts = params.singleEnd ? "${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
  """
  localIndex=`find -L ./ -name "*.rev.1.ebwt" | sed 's/.rev.1.ebwt//'`
  echo \$(bowtie --version | awk 'NR==1{print "bowtie "\$3}') > versions.txt
  bowtie ${params.bowtieOpts} \\
         -p ${task.cpus} \\
         --un ${prefix}_norRNA.fastq \\
         --sam \${localIndex} \\
         ${inputOpts} \\
         ${prefix}.sam  2> ${prefix}.log && \
  gzip -f ${prefix}_norRNA*.fastq 
  rm -f ${prefix}.sam
  """
}


