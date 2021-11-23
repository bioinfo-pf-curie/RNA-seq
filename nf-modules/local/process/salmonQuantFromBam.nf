/*
 * Salmon quant from BAM file
 * External parameters :
 * @ params.singleEnd : is the data single-end ?
 * @ params.gencode : is the annotation from Gencode ?
 * @ params.salmonQuantOpts : addition option for Salmon quantification
 */

process salmonQuantFromBam {
    tag "$prefix"
    label "salmon"
    label "medCpu"
    label "medMem"

    publishDir "${params.outDir}/counts", mode: 'copy'

    input:
    tuple val(prefix), path(bam), val(strandness)
    path(transcriptsFasta)
    path(gtf)

    output:
    path("${prefix}"), emit: results
    path("versions.txt"), emit: versions

    script:
    gencodeOpts = params.gencode ? '--gencode' : ''
    strandOpts = params.singleEnd ? 'U' : 'IU'
    if (strandness == 'forward') {
      strandOpts = params.singleEnd ? 'SF' : 'ISF'
    } else if (strandness == 'reverse') {
      strandOpts = params.singleEnd ? 'SR' : 'ISR'
    }    
    """
    echo \$(salmon --version 2>&1) > versions.txt
    salmon quant \\
      --libType=$strandOpts \\
      -a ${bam} \\
      --threads ${task.cpus} \\
      -t ${transcriptsFasta} \\
      --geneMap ${gtf} \\
      ${params.salmonQuantOpts} \\
      ${gencodeOpts} \\
      -o ${prefix}
    """
}
