process umitoolsExtract {
    tag "$meta.id"
    label 'umitools'
    label 'minCpu'
    label 'medMem'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    tuple val(meta), path("*.log")     , emit: log
    path  "versions.txt"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        umi_tools \\
            extract \\
            -I $reads \\
            -S ${prefix}.umi_extract.fastq.gz \\
            $args \\
            > ${prefix}.umi_extract.log

        umi_tools --version | cut -f1,3 -d" " &> versions.txt
        """
    }  else {
        """
        umi_tools \\
            extract \\
            -I ${reads[0]} \\
            --read2-in=${reads[1]} \\
            -S ${prefix}.umi_extract_1.fastq.gz \\
            --read2-out=${prefix}.umi_extract_2.fastq.gz \\
            $args \\
            > ${prefix}.umi_extract.log

        umi_tools --version | cut -f1,3 -d" " &> versions.txt
        """
    }
}