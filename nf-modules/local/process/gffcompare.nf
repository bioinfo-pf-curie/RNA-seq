/*
 * gffCompare - compare reference and new GTF files
 */

process gffcompare {
    label 'gffcompare'
    label 'minCpu'
    label 'medMem'

    publishDir "${params.outDir}/gffcompare/${source}/", mode: 'copy'

    input:
    path(denovoGtf)
    path(refGtf)
    val(source)

    output:
    path("gffcompare*")   , emit: results
    path("*combined.gtf") , optional: true, emit: combinedGtf
    path("versions.txt")  , emit: versions

    script:
    """
    gffcompare \\
      -r ${refGtf} \\
      -o gffcompare_${source} \\
      ${denovoGtf}

    echo \$(gffcompare --version 2>&1) > versions.txt
    """
}
