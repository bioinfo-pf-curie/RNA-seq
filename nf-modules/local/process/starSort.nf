process starSort {
    tag "$prefix"
    label 'samtools'
    label 'medCpu'
    label 'medMem'
    publishDir "${params.outDir}/mapping", mode: 'copy'
 
    input:
    tuple val(prefix), path(LogFinalOut), path (starBam)

    output:
    tuple path("${prefix}Log.final.out"), path ("*.{bam,bam.bai}"), emit: starAligned
    path("v_samtools.txt")                                        , emit: version
    path "${prefix}_sorted.bam.bai"

    script:
    """
    samtools --version &> v_samtools.txt
    samtools sort  \\
        -@  ${task.cpus}  \\
        -m ${params.sortMaxMemory} \\
        -o ${prefix}_sorted.bam  \\
        ${starBam}
    samtools index ${prefix}_sorted.bam
    """
    }
    