/*
 * Scallop - Reference-guided de novo isoform assembly
 * https://github.com/Kingsford-Group/scallop
 * External parameters :
 * @ params.scallopOpts : Additional Scallop parameters
 */

process scallop {
    tag "$prefix"
    label 'scallop'
    label 'medCpu'
    label 'medMem'
    publishDir "${params.outDir}/scallop", mode: 'copy'

    input:
    tuple val(prefix), path(bam), path(bai), val(strandness) // Channel [prefix, bam, bai, strandness]

    output:
    tuple val(prefix), path("*Transcripts.gtf"), emit: transcriptGtf
    path  "versions.txt"                       , emit: versions

    script:
    def strandOpts = '--library_type unstranded'
    if (strandness == 'forward') {
        strandOpts = '--library_type second'
    } else if (strandness == 'reverse') {
        strandOpts = '--library_type first'
    }
    """
    scallop \\
        -i ${bam} \\
        ${strandOpts} \\
        ${params.scallopOpts} \\
        -o ${prefix}.scallopTranscripts.gtf

    echo "scallop "\$(scallop --version 2>&1) > versions.txt
    """
}
