/*
 * Saturation Curves
 */

process preseq {
  tag "${bam[0]}"
  label 'preseq'
  label 'minCpu'
  label 'medMem'
  publishDir "${params.outDir}/preseq", mode: 'copy'

  when:
  !params.skipQC && !params.skipSaturation

  input:
  path bam

  output:
  path "*ccurve.txt"  , emit: chPreseqResults
  path("v_preseq.txt"), emit: chPreseqVersion

  script:
  peOpts = params.singleEnd ? '' : '-pe'
  """
  preseq &> v_preseq.txt
  preseq lc_extrap -seed 1 -v -B ${bam[0]} ${peOpts} -o ${bam[0].baseName}_extrap_ccurve.txt -e 200e+06 -seg_len 100000000
  """
}

