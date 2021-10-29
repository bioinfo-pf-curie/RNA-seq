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
    val strandness

    output:
    tuple val(prefix), path("${prefix}.bam"), emit: bam
    path "${prefix}.hisat2_summary.txt"     , emit: logs
    path ("versions.txt")                   , emit: versions

    script:
    indexBase = hs2Index[0].toString() - ~/.\d.ht2/
    def strandOpts = ''
    if (strandness=='forward'){
      strandOpts = params.singleEnd ? '--rna-strandness F' : '--rna-strandness FR'
    } else if (strandness=='reverse'){
      strandOpts = params.singleEnd ? '--rna-strandness R' : '--rna-strandness RF'
    }
    inputOpts = params.singleEnd ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    echo \$(hisat2 --version | awk 'NR==1{print "hisat2 "\$3}') > versions.txt
    hisat2 -x $indexBase \\
           ${inputOpts} \\
           ${strandOpts} \\
           --known-splicesite-infile $alignmentSplicesites \\
           -p ${task.cpus} \\
           --met-stderr \\
           --new-summary \\
	   --rg-id ${prefix} --rg SM:${prefix} --rg PL:ILLUMINA \\
           --summary-file ${prefix}.hisat2_summary.txt \\
           | samtools view -bS -F 256 - > ${prefix}.bam
    """
  }
  
  
