/*
 * Salmon quant from BAM file
 */

process salmonQuantFromBam {
    tag "$prefix"
    label "salmon"
    label "medCpu"
    label "medMem"

    publishDir "${params.outDir}/counts", mode: 'copy',

    input:
    tuple val(prefix), path(bam), val(strandness)
    path  transcriptFasta
    path  gtf

    output:
    tuple val(${prefix}), path("${prefix}"), emit: results
    path("versions.txt"), emit: versions

    script:


    def reference   = "--index $index"
    def input_reads = meta.single_end ? "-r $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    if (alignment_mode) {
        reference   = "-t $transcript_fasta"
        input_reads = "-a $reads"
    }

    def strandedness_opts = [
        'A', 'U', 'SF', 'SR',
        'IS', 'IU' , 'ISF', 'ISR',
        'OS', 'OU' , 'OSF', 'OSR',
        'MS', 'MU' , 'MSF', 'MSR'
    ]
    def strandedness =  'A'
    if (lib_type) {
        if (strandedness_opts.contains(lib_type)) {
            strandedness = lib_type
        } else {
            log.info "[Salmon Quant] Invalid library type specified '--libType=${lib_type}', defaulting to auto-detection with '--libType=A'."
        }
    } else {
        strandedness = meta.single_end ? 'U' : 'IU'
        if (meta.strandedness == 'forward') {
            strandedness = meta.single_end ? 'SF' : 'ISF'
        } else if (meta.strandedness == 'reverse') {
            strandedness = meta.single_end ? 'SR' : 'ISR'
        }
    }
    """
    salmon quant \\
        --geneMap $gtf \\
        --threads $task.cpus \\
        --libType=$strandedness \\
        $reference \\
        $input_reads \\
        $options.args \\
        -o $prefix
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
