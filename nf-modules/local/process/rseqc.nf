/*
 * RSeQC - infer_experiment.py
 */

process rseqc {
    tag "${prefix - '_subsample'}"
    label 'rseqc'
    label 'medCpu'
    label 'lowMem'
    publishDir "${params.outDir}/strandness" , mode: 'copy',
      saveAs: {filename ->
        if (filename.indexOf(".txt") > 0) "$filename"
        else null }

    input:
    tuple val(prefix), path(bamRseqc)
    path bed12 

    output:
    path "${prefix}*.{txt,pdf,r,xls}", emit: rseqcResults
    path("versions.txt")             , emit: versions
    stdout emit: strandnessResults

    script:
    """
    echo \$(infer_experiment.py --version | awk '{print "rseqc "\$2}') > versions.txt    
    infer_experiment.py -i $bamRseqc -r $bed12 > ${prefix}.txt
    res=\$(parse_rseq_output.sh ${prefix}.txt)
    echo "$prefix,\$res"> ${prefix}_strandness.txt
    cat ${prefix}_strandness.txt
    """  
  }
