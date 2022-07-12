/*
 * gffCompare - compare reference and new GTF files
 */

process gffcompare {
  label 'gffcompare'
  label 'minCpu'
  label 'medMem'

  input:
  tuple val(meta), path(denovoGtf)
  path(refGtf)

  output:
  tuple val(meta), path("gffcompare_*"), emit: results
  path("*.stats"), emit: mqc
  path("*combined.gtf"), optional: true, emit: combinedGtf
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  gffcompare \\
    -r ${refGtf} \\
    -o gffcompare_${prefix} \\
    ${denovoGtf}

  echo \$(gffcompare --version 2>&1) > versions.txt
  """
}
