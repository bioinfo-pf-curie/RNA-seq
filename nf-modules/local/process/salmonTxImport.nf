/*
 * Salmon
 * Import Salmon results and generate counts table using TxImport
 */

process salmonTxImport {
  label "r"
  label "medCpu"
  label "medMem"

  input:
  path("salmon/*")
  path(salmonTx2gene)

  output:
  path "*gene_tpm.csv"                 , emit: tpmGene
  path "*gene_counts.csv"              , emit: countsGene
  path "*gene_counts_length_scaled.csv", emit: countsGeneLengthScaled
  path "*gene_counts_scaled.csv"       , emit: countsGeneScaled
  path "*transcript_tpm.csv"           , emit: tpmTranscript
  path "*transcript_counts.csv"        , emit: countsTranscript
  path("versions.txt"), emit: versions

  script:
  """
  echo \$(R --version | awk 'NR==1{print \$1,\$3}') > versions.txt
  salmon_tximport.r \\
    NULL \\
    salmon \\
    salmon.merged
  """
}
