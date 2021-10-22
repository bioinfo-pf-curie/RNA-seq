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
    path("v_samtools.txt") , emit: version

    script:
    """
    samtools --version &> v_samtools.txt
    samtools sort  \\
        -@  ${task.cpus}  \\
        -m ${params.sortMaxMemory} \\
        -o ${prefix}_sorted.bam  \\
        ${bam}
    """
}
    