
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
    path("v_rseqc.txt")              , emit: version
    stdout emit: strandedResults

    script:
    """
    infer_experiment.py --version &> v_rseqc.txt    
    infer_experiment.py -i $bamRseqc -r $bed12 > ${prefix}.txt
    parse_rseq_output.sh ${prefix}.txt > ${prefix}_strandness.txt
    cat ${prefix}_strandness.txt
    """  
  }
