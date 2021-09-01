/*
 * Generate bigwig file
 */

process bigWig {
  tag "${prefix}"
  label 'deeptools'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outDir}/bigWig", mode: "copy",
    saveAs: {filename ->
    	     if ( filename.endsWith(".bigwig") ) "$filename"
             else null}

  when:
  !params.skipBigwig

  input:
  path(bam)
  val parseRes

  output:
  path('*.bigwig')       , emit: bigWig
  path("v_deeptools.txt"), emit: version

  script:
  prefix = bam[0].toString() - ~/(_sorted)?(.bam)?$/
  strandOpt = parseRes == 'forward' ? '--filterRNAstrand forward' : parseRes == 'reverse' ? '--filterRNAstrand reverse' : ''
  """
  bamCoverage --version &> v_deeptools.txt
  bamCoverage -b ${bam[0]} \\
              -o ${prefix}_cpm.bigwig \\
              -p ${task.cpus} \\
              ${strandOpt} \\
	      --normalizeUsing CPM \\
	      --skipNonCoveredRegions
  """
}