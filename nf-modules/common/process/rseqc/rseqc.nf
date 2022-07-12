/*
 * RSeQC - infer_experiment.py
 */

process rseqc {
  tag "${meta.id}"
  label 'rseqc'
  label 'medCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(bamRseqc)
  path bed12 

  output:
  path "${meta.id}*.{txt,pdf,r,xls}", emit: rseqcResults
  path("versions.txt"), emit: versions
  stdout emit: strandnessResults

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args   = task.ext.args ?: ''
  """
  echo \$(infer_experiment.py --version | awk '{print "rseqc "\$2}') > versions.txt    
  infer_experiment.py -i $bamRseqc -r $bed12 ${args} > ${prefix}.txt
  res=\$(parse_rseq_output.sh ${prefix}.txt)
  echo "${prefix},\$res"> ${prefix}_strandness.txt
  cat ${prefix}_strandness.txt
  """  
}
