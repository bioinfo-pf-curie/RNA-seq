process outputDocumentation {
    label 'python'
    label 'minCpu'
    label 'minMem'
    publishDir "${params.summaryDir}/", mode: 'copy'

    input:
    path outputDocs 
    path images 

    output:
    path "results_description.html"

    script:
    """
    markdown_to_html.py $outputDocs -o results_description.html
    """
}
