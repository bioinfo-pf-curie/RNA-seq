process hisat2Align {
    tag "$prefix"
    label 'hisat2'
    label 'highCpu'
    label 'highMem'
    publishDir "${params.outDir}/mapping", mode: 'copy',
      saveAs: {filename ->
        if (filename.indexOf(".hisat2_summary.txt") > 0) "logs/$filename"
        else if (params.saveAlignedIntermediates) filename
        else null
      }

    input:
    tuple val(prefix), path(reads)
    path hs2Index
    path alignmentSplicesites
    val parseRes

    output:
    path "${prefix}.bam"                , emit: hisat2Bam
    path "${prefix}.hisat2_summary.txt" , emit: alignmentLogs
    path ("v_hisat2.txt")               , emit: version

    script:
    indexBase = hs2Index[0].toString() - ~/.\d.ht2/
    def rnastrandness = ''
    if (parseRes=='forward'){
        rnastrandness = params.singleEnd ? '--rna-strandness F' : '--rna-strandness FR'
    } else if (parseRes=='reverse'){
        rnastrandness = params.singleEnd ? '--rna-strandness R' : '--rna-strandness RF'
    }
    inputOpts = params.singleEnd ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    hisat2 --version &> v_hisat2.txt
    hisat2 -x $indexBase \\
           ${inputOpts} \\
           $rnastrandness \\
           --known-splicesite-infile $alignmentSplicesites \\
           -p ${task.cpus} \\
           --met-stderr \\
           --new-summary \\
	   --rg-id ${prefix} \\
           --summary-file ${prefix}.hisat2_summary.txt \\
           | samtools view -bS -F 4 -F 256 - > ${prefix}.bam
    """
  }
  
  
