/*
 * gffCompare
 */

process gffcompare {
    label 'gffcompare'
    label 'minCpu'
    label 'medMem'

    publishDir "${params.outDir}/gffcompare", mode: 'copy'

    input:
    path(denovoGtf)
    path(refGtf)

    output:
    path("gffcompare*")   , emit: results
    path("*combined.gtf") , optional: true, emit: combinedGtf
    path("versions.txt")  , emit: versions

    script:
    """
    gffcompare \\
      -r ${refGtf} \\
      -o gffcompare \\
      ${denovoGtf}

    echo \$(gffcompare --version 2>&1) > versions.txt
    """
}
