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
  tuple val(prefix), path(bam), path(bai)

  output:
  path("*ccurve.txt"), emit: results
  path("versions.txt"), emit: versions

  script:
  peOpts = params.singleEnd ? '' : '-pe'
  """
  echo \$(preseq 2>&1 | awk '\$0~"Version"{print "Preseq",\$2}') > versions.txt
  preseq lc_extrap -seed 1 -v -B ${bam} ${peOpts} -o ${prefix}_extrap_ccurve.txt -e 200e+06 -seg_len 100000000
  """
}

