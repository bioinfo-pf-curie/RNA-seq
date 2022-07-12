/*
 * Hisat2 genome mapping
 */

process hisat2Align {
  tag "${meta.id}"
  label 'hisat2'
  label 'highCpu'
  label 'extraMem'

  input:
  tuple val(meta), path(reads)
  path hs2Index
  path alignmentSplicesites

  output:
  tuple val(meta), path("${prefix}.bam"), emit: bam
  path "*hisat2_summary.txt"     , emit: logs
  path ("versions.txt")                   , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  indexBase = hs2Index[0].toString() - ~/.\d.ht2/
  def strandOpts = ''
  if (meta.strandness=='forward'){
    strandOpts = meta.singleEnd ? '--rna-strandness F' : '--rna-strandness FR'
  } else if (meta.strandness=='reverse'){
    strandOpts = meta.singleEnd ? '--rna-strandness R' : '--rna-strandness RF'
  }
  inputOpts = meta.singleEnd ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  echo \$(hisat2 --version | awk 'NR==1{print "hisat2 "\$3}') > versions.txt
  localIndex=`find -L ./ -name "*.1.ht2" | sed 's/.1.ht2//'`
  hisat2 -x \${localIndex} \\
         ${inputOpts} \\
         ${strandOpts} \\
         --known-splicesite-infile $alignmentSplicesites \\
         -p ${task.cpus} \\
         --met-stderr \\
         --new-summary \\
         --rg-id ${prefix} --rg SM:${prefix} --rg PL:ILLUMINA \\
         --summary-file ${prefix}.hisat2_summary.txt \\
         ${args} \\
         | samtools view -bS -F 256 - > ${prefix}.bam
  """
}
  
  
