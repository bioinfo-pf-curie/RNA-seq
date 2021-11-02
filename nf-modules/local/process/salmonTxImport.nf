/*
 * Salmon Tx Import
 */

process salmonTxImport {
    tag "$prefix"
    label "R"
    label "medCpu"
    label "medMem"

    publishDir "${params.outDir}/counts", mode: 'copy'

    input:
    tuple val(prefix), path(salmonResults), path(salmonTx2gene)

    output:
    path("*.tsv"), emit: results
    path("versions.txt"), emit: versions

    script:
    """
    echo \$(R --version | awk 'NR==1{print \$1,\$3}') > versions.txt
    salmon_tximport.r \
      ${salmonTx2gene} \
      ${salmonResults} \
      ${prefix}
    """
}
