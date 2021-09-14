/* 
 * Qualimap
 */

process qualimap {
  tag "${bam[0].baseName}"
  label 'qualimap'
  label 'minCpu'
  label 'medMem'
  publishDir "${params.outDir}/qualimap/" , mode: 'copy'

  when:
  !params.skipQualimap

  input:
  path bam
  path gtf
  val stranded

  output:
  path ("${bam[0].baseName}"), emit: chQualimapResults
  path ("v_qualimap.txt")    , emit: chQualimapVersion

  script:
  peOpts = params.singleEnd ? '' : '-pe'
  memory     = task.memory.toGiga() + "G"
  strandnessOpts = 'non-strand-specific'
  if (stranded == 'forward') {
    strandnessOpts = 'strand-specific-forward'
  } else if (stranded == 'reverse') {
    strandnessOpts = 'strand-specific-reverse'
  }
  """
  mkdir tmp
  export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
  qualimap \\
    --java-mem-size=$memory \\
    rnaseq \\
    -bam $bam[0] \\
    -gtf $gtf \\
    -p $strandnessOpts \\
    $peOpts \\
    -outdir ${bam[0].baseName}
  echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//' > v_qualimap.txt
  """
}
