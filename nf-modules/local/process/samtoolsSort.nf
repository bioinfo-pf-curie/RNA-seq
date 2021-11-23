/*
 * Samtools - Sort
 * External parameters :
 * @ params.sortMaxMemory : maximum sort memory
 */

process samtoolsSort {
    tag "$prefix"
    label 'samtools'
    label 'medCpu'
    label 'medMem'
    publishDir "${params.outDir}/mapping", mode: 'copy'
 
    input:
    tuple val(prefix), path (bam)

    output:
    tuple val(prefix), path ("*_sorted.bam"), emit: bam
    path("versions.txt") , emit: versions

    script:
    maxMem = params.sortMaxMemory ?: '900M' 
    """
    echo \$(samtools --version | head -1 ) > versions.txt
    samtools sort  \\
        -@  ${task.cpus}  \\
        -m ${maxMem} \\
        -o ${prefix}_sorted.bam  \\
        ${bam}
    """
}
    