process hisat2Sort {
    tag "${hisat2Bam.baseName}"
    label 'samtools'
    label 'medCpu'
    label 'medMem'  
    publishDir "${params.outDir}/mapping", mode: 'copy'

    input:
    path hisat2Bam

    output:
    path ('*sorted.{bam,bam.bai}')             , emit: chBam
    file "${hisat2Bam.baseName}_sorted.bam.bai"
    file("v_samtools.txt")                     , emit: chSamtoolsVersionSort 

    script:
    def availMem = task.memory ? "-m ${task.memory.toBytes() / task.cpus}" : ''
    """
    samtools --version &> v_samtools.txt
    samtools sort \\
             ${hisat2Bam} \\
             -@ ${task.cpus} $availMem \\
             -m ${params.sortMaxMemory} \\
             -o ${hisat2Bam.baseName}_sorted.bam
    samtools index ${hisat2Bam.baseName}.sorted.bam
    """
  }