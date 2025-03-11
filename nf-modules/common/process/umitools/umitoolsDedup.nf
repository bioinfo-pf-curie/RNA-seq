process umitoolsDedup {
    tag "$meta.id"
    label 'umitools'
    label 'minCpu'
    label 'extraMem'

    input:
    tuple val(meta), path(bam), path(bai)
    val stats

    output:
    tuple val(meta), path("${prefix}.bam")     , emit: bam
    tuple val(meta), path("*.log")             , emit: log
    tuple val(meta), path("*edit_distance.tsv"), optional:true, emit: tsvDistance
    tuple val(meta), path("*per_umi.tsv")      , optional:true, emit: tsvUmi
    tuple val(meta), path("*per_position.tsv") , optional:true, emit: tsvPosition
    path  "versions.txt"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def paired = meta.singleEnd ? "" : "--paired"
    statsOpt = stats ? "--output-stats ${prefix}" : ""
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    if (!(args ==~ /.*--random-seed.*/)) {args += " --random-seed=100"}
    """
    PYTHONHASHSEED=0 umi_tools \\
        dedup \\
        -I $bam \\
        -S ${prefix}.bam \\
        -L ${prefix}.log \\
        $statsOpt \\
        $paired \\
        $args

    umi_tools --version | cut -f1,3 -d" " &> versions.txt
    """
}