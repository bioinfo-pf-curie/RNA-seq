/*
 * Generate bigwig file
 */

process bigWig {
  tag "${prefix}"
  label 'deeptools'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(prefix), path(bam), path(bai), val(strandness)

  output:
  path('*.bigwig') , emit: bigWig
  path("versions.txt"), emit: versions

  script:
  strandOpts = strandness == 'forward' ? '--filterRNAstrand forward' : strandness == 'reverse' ? '--filterRNAstrand reverse' : ''
  """
  echo \$(bamCoverage --version) > versions.txt
  bamCoverage -b ${bam} \\
              -o ${prefix}_cpm.bigwig \\
              -p ${task.cpus} \\
              ${strandOpts} \\
	      --normalizeUsing CPM \\
	      --skipNonCoveredRegions
  """
}