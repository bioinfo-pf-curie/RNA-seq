/*
 * Stringtie - Merge
 * https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
 */

process stringtieMerge {
    label 'stringtie'
    label 'medCpu'
    label 'medMem'
    publishDir "${params.outDir}/stringtie", mode: 'copy'

    input:
    path(stringtieGtf)
    path(gtf)

    output:
    tuple val("stringtieMerge"), path("mergedTranscripts.gtf"), emit: mergedGtf
    path("versions.txt"), emit: versions

    script:
    """
    echo -e ${stringtieGtf} | tr " " "\n" > listofgtf.tsv
    stringtie \\
        --merge \\
        -G ${gtf} \\
        -o mergedTranscripts.gtf \\
        -p ${task.cpus} \\
	listofgtf.tsv

    echo "stringtie "\$(stringtie --version 2>&1) > versions.txt
    """
}
