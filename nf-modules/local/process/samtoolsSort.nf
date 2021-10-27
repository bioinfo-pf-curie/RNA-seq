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
    """
    echo \$(samtools --version | head -1 ) > versions.txt
    samtools sort  \\
        -@  ${task.cpus}  \\
        -m ${params.sortMaxMemory} \\
        -o ${prefix}_sorted.bam  \\
        ${bam}
    """
}
    