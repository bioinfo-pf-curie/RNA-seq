/*
 * Duplicates
 */

process markDuplicates {
  tag "${bam[0]}"
  label 'picard'
  label 'minCpu'
  label 'medMem'
  publishDir "${params.outDir}/markDuplicates", mode: 'copy',
    saveAs: {filename -> 
      if (filename.indexOf("_metrics.txt") > 0) "metrics/$filename" 
      else if (params.saveAlignedIntermediates) filename
    }

  when:
  !params.skipQC && !params.skipDupradar

  input:
  path bam

  output:
  path('*markDups.{bam,bam.bai}'), emit: BamMd
  path('*markDups_metrics.txt')  , emit: PicardResults
  path('v_picard.txt')           , emit: PicardVersion

  script:
  markdupMemOption = "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
  """
  echo \$(picard MarkDuplicates --version 2>&1) &> v_picard.txt
  picard ${markdupMemOption} -Djava.io.tmpdir=${params.tmpDir} MarkDuplicates \\
      MAX_RECORDS_IN_RAM=50000 \\
      INPUT=${bam[0]} \\
      OUTPUT=${bam[0].baseName}.markDups.bam \\
      METRICS_FILE=${bam[0].baseName}.markDups_metrics.txt \\
      REMOVE_DUPLICATES=false \\
      ASSUME_SORTED=true \\
      PROGRAM_RECORD_ID='null' \\
      VALIDATION_STRINGENCY=LENIENT
  samtools index ${bam[0].baseName}.markDups.bam
  """
}