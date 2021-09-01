/*
 * rRNA mapping 
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

  when:
  !params.skipRrna && params.rrna

  input:
  tuple val(prefix), file(reads)
  path annot

  output:
  tuple val(prefix), path("*fastq.gz"), emit: chRrnaMappingRes
  tuple val(prefix), path("*.sam")    , emit: rnaSam
  path "*.log"                        , emit: logs
  path("v_bowtie.txt")                , emit: version

  script:
  inputOpts = params.singleEnd ? "${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
  """
  bowtie --version &> v_bowtie.txt
  bowtie ${params.bowtieOpts} \\
         -p ${task.cpus} \\
         --un ${prefix}_norRNA.fastq \\
         --sam ${params.rrna} \\
         ${inputOpts} \\
         ${prefix}.sam  2> ${prefix}.log && \
  gzip -f ${prefix}_norRNA*.fastq 
  """
}

