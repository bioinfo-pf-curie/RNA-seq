/*
 * gffCompare - compare reference and new GTF files
 */

process gffcompare {
    label 'gffcompare'
    label 'minCpu'
    label 'medMem'

    publishDir "${params.outDir}/gffcompare/", mode: 'copy'

    input:
    tuple val(prefix), path(denovoGtf)
    path(refGtf)

    output:
    tuple val(prefix), path("gffcompare_${prefix}*"), emit: results
    path("*.stats"), emit: mqc
    path("*combined.gtf"), optional: true, emit: combinedGtf
    path("versions.txt"), emit: versions

    script:
    """
    gffcompare \\
      -r ${refGtf} \\
      -o gffcompare_${prefix} \\
      ${denovoGtf}

    echo \$(gffcompare --version 2>&1) > versions.txt
    """
}
