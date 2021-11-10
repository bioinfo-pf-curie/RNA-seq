/*
 * Reference-guided de novo isoform assembly
 * https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
 */

process stringtie {
    tag "$prefix"
    label 'stringtie'
    label 'medCpu'
    label 'medMem'
    publishDir "${params.outDir}/stringtie", mode: 'copy'

    input:
    tuple val(prefix), path(bam), path(bai), val(strandness) // Channel [prefix, bam, bai, strandness]
    path(gtf)

    output:
    tuple val(prefix), path("*.coverage.gtf")   , emit: coverageGtf
    tuple val(prefix), path("*.transcripts.gtf"), emit: transcriptGtf
    tuple val(prefix), path("*.abundance.txt")  , emit: abundance
    tuple val(prefix), path("*.ballgown")       , emit: ballgown
    path  "versions.txt"                        , emit: versions

    script:
    def strandOpts = ''
    if (strandness == 'forward') {
        strandOpts = '--fr'
    } else if (strandness == 'reverse') {
        strandOpts = '--rf'
    }
    """
    stringtie \\
        $bam \\
        $strandOpts \\
        -G $gtf \\
        -o ${prefix}.transcripts.gtf \\
        -A ${prefix}.gene.abundance.txt \\
        -C ${prefix}.coverage.gtf \\
        -b ${prefix}.ballgown \\
        -p $task.cpus \\
        $params.stringtieOpts

    echo "stringtie "\$(stringtie --version 2>&1) > versions.txt
    """
}
