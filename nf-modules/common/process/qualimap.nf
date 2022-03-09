/* 
 * Qualimap - quality controls on RNA-seq BAM file
 */

process qualimap {
  tag "${meta.id}"
  label 'qualimap'
  label 'minCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bam), path(bai)
  path gtf

  output:
  path ("${prefix}"), emit: results
  path ("versions.txt"), emit: versions

  script:
  def peOpts = meta.singleEnd ? '' : '-pe'
  def memory = task.memory.toGiga() + "G"
  def strandnessOpts = 'non-strand-specific'
  if (meta.stranded == 'forward') {
    strandnessOpts = 'strand-specific-forward'
  } else if (meta.stranded == 'reverse') {
    strandnessOpts = 'strand-specific-reverse'
  }
  prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  mkdir tmp
  export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
  qualimap \\
    --java-mem-size=$memory \\
    rnaseq \\
    -bam ${bam} \\
    -gtf ${gtf} \\
    -p ${strandnessOpts} \\
    ${peOpts} \\
    -outdir ${prefix}
  echo \$(qualimap -h 2>&1 |  grep "QualiMap" | sed -e 's/v.//') > versions.txt
  """
}
