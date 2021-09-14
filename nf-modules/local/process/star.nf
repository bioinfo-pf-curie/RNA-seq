  // STAR

  process star {
    tag "$prefix"
    label 'star'
    label 'highCpu'
    label 'highMem'
    publishDir "${params.outDir}/mapping", mode: 'copy',
      saveAs: {filename ->
        if (filename.indexOf(".bam") == -1) "logs/$filename"
        else if (params.saveAlignedIntermediates) filename
        else null
      }
    publishDir "${params.outDir}/counts", mode: 'copy',
      saveAs: {filename ->
        if (filename.indexOf("ReadsPerGene.out.tab") > 0) "$filename"
        else null
      }

    input:
    tuple val(prefix), path(reads)
    path index
    file gtf

    output:
    tuple val(prefix), path ("*Log.final.out"), path ('*.bam'), emit: chStarSam
    path "*.out"                                              , emit: alignmentLogs
    path "*.out.tab"                                          , emit: starLogCounts
    path "*Log.out"                                           , emit: starLog
    path "*ReadsPerGene.out.tab" optional true                , emit: starCounts
    path("v_star.txt")                                        , emit: version

    script:
    def starCountOpt = params.counts == 'star' && params.gtf ? params.starOptsCounts : ''
    def starGtfOpt = params.gtf ? "--sjdbGTFfile $gtf" : ''
    """
    STAR --version &> v_star.txt
    STAR --genomeDir $index \\
         ${starGtfOpt} \\
         --readFilesIn $reads  \\
         --runThreadN ${task.cpus} \\
         --runMode alignReads \\
         --outSAMtype BAM Unsorted  \\
         --readFilesCommand zcat \\
         --runDirPerm All_RWX \\
         --outTmpDir "${params.tmpDir}/rnaseq_\$(date +%d%s%S%N)"\\
         --outFileNamePrefix $prefix  \\
         --outSAMattrRGline ID:$prefix SM:$prefix LB:Illumina PL:Illumina  \\
         ${params.starOptions} \\
	 --limitOutSJcollapsed 5000000 \\
	 ${starCountOpt}
    """
  }
