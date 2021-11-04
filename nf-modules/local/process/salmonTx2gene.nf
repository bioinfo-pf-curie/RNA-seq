/*
 * Salmon Tx2gene
 */

process salmonTx2gene {
    tag "$gtf"
    label "python"
    label "medCpu"
    label "medMem"

    publishDir "${params.outDir}/counts", mode: 'copy'

    input:
    path("salmon/*")
    path(gtf)

    output:
    path("salmon_tx2gene.tsv"), emit: results
    path("versions.txt"), emit: versions

    script:
    """
    echo \$(python --version 2>&1) > versions.txt
    salmon_tx2gene.py \
      --gtf ${gtf} \
      --salmon salmon \
      --id "gene_id" \
      --extra "gene_name" \
      -o salmon_tx2gene.tsv 

    """
}
