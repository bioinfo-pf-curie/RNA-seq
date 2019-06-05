#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2019
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data.
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

This script is based on the nf-core guidelines. See https://nf-co.re/ for more information
*/


/*
========================================================================================
                         mypipeline
========================================================================================
 RNA-seq Analysis Pipeline.
 #### Homepage / Documentation
 https://gitlab.curie.fr/rnaseq
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info """
    rnaseq v${workflow.manifest.version}
    =======================================================

    Usage:
    nextflow run rnaseq --reads '*_R{1,2}.fastq.gz' --genome hg19 -profile conda

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Configuration profile to use. test / curie / conda / docker / singularity

    Options:
      --genome                      Name of genomes reference
      --singleEnd                   Specifies that the input is single end reads

    Strandedness:
      --stranded                    Library strandness ['auto', 'yes', 'reverse', 'no']. Default: 'auto'

    Mapping:
      --aligner                     Tool for read alignments ['star', 'hisat2', 'tophat2']. Default: 'star'

    Counts:
      --counts                      Tool to use to estimate the raw counts per gene ['star', 'featureCounts', 'HTseqCounts']. Default: 'star'

    References:                     If not specified in the configuration file or you wish to overwrite any of the references.
      --star_index                  Path to STAR index
      --hisat2_index                Path to HiSAT2 index
      --tophat2_index		    Path to TopHat2 index
      --gtf                         Path to GTF file
      --bed12                       Path to gene bed12 file
      --saveAlignedIntermediates    Save the BAM files from the Aligment step  - not done by default

    Other options:
      --metadata                    Add metadata file for multiQC report
      --outdir                      The output directory where the results will be saved
      -w/--work-dir                 The temporary directory where intermediate data will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    QC options:
      --skip_qc                     Skip all QC steps apart from MultiQC
      --skip_rrna                   Skip rRNA mapping
      --skip_fastqc                 Skip FastQC
      --skip_preseq                 Skip Preseq
      --skip_dupradar               Skip dupRadar (and Picard MarkDups)
      --skip_read_dist              Skip read distribution step
      --skip_expan                  Skip exploratory analysis
      --skip_multiqc                Skip MultiQC

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */


// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Configurable reference genomes
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
  }

// Reference index path configuration
// Define these here - after the profiles are loaded with the iGenomes paths
params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
params.bowtie2_index = params.genome ? params.genomes[ params.genome ].bowtie2 ?: false : false
params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2 ?: false : false
params.rrna = params.genome ? params.genomes[ params.genome ].rrna ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.bed12 = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Stage config files
ch_multiqc_config = Channel.fromPath(params.multiqc_config)
//ch_multiqc_logo = Channel.fromPath("$baseDir/assets/institut_curie.jpg")
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")
ch_pca_header = Channel.fromPath("$baseDir/assets/pca_header.txt")
ch_heatmap_header = Channel.fromPath("$baseDir/assets/heatmap_header.txt")

/*
 * CHANNELS
 */

// Validate inputs
if (params.aligner != 'star' && params.aligner != 'hisat2' && params.aligner != 'tophat2'){
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star', 'hisat2', 'tophat2'"
}
if (params.counts != 'star' && params.counts != 'featureCounts' && params.counts != 'HTseqCounts'){
    exit 1, "Invalid counts option: ${params.counts}. Valid options: 'star', 'featureCounts', 'HTseqCounts'"
}
if (params.counts == 'star' && params.aligner != 'star'){
    exit 1, "Cannot run STAR counts without STAR aligner."
}
if (params.stranded != 'auto' && params.stranded != 'reverse' && params.stranded != 'yes' && params.stranded != 'no'){
    exit 1, "Invalid stranded option: ${params.stranded}. Valid options: 'auto', 'reverse', 'yes', 'no'"
}

if( params.star_index && params.aligner == 'star' ){
    star_index = Channel
        .fromPath(params.star_index)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
}
else if ( params.hisat2_index && params.aligner == 'hisat2' ){
    hs2_indices = Channel
        .fromPath("${params.hisat2_index}*")
        .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }
}
else if ( params.bowtie2_index && params.aligner == 'tophat2' ){
    Channel.fromPath("${params.bowtie2_index}*")
        .ifEmpty { exit 1, "TOPHAT2 index not found: ${params.bowtie2_index}" }
        .into { tophat2_indices}
}
else {
    exit 1, "No reference genome specified!"
}

if( params.gtf ){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .into { gtf_star; gtf_dupradar; gtf_featureCounts; gtf_HTseqCounts; gtf_tophat; gtf_table }
} else {
    exit 1, "No GTF annotation specified!"
}

if( params.bed12 ){
    bed12 = Channel
        .fromPath(params.bed12)
        .ifEmpty { exit 1, "BED12 annotation file not found: ${params.bed12}" }
        .into { bed_rseqc; bed_read_dist} 
}

if( params.rrna ){
    Channel
        .fromPath(params.rrna)
        .ifEmpty { exit 1, "rRNA annotation file not found: ${params.rrna}" }
        .set { rrna_annot }
}

if( params.aligner == 'hisat2' && params.splicesites ){
    Channel
        .fromPath(params.bed12)
        .ifEmpty { exit 1, "HISAT2 splice sites file not found: $alignment_splicesites" }
        .into { indexing_splicesites; alignment_splicesites }
}

if ( params.metadata ){
   Channel
       .fromPath( params.metadata )
       .ifEmpty { exit 1, "Metadata file not found: ${params.metadata}" }
       .set { ch_metadata }
}

/*
 * Create a channel for input read files
 */
if(params.readPaths){
    if(params.singleEnd){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { raw_reads_fastqc; raw_reads_star; raw_reads_hisat2; raw_reads_tophat2; raw_reads_rna_mapping; raw_reads_prep_rseqc; }
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { raw_reads_fastqc; raw_reads_star; raw_reads_hisat2; raw_reads_tophat2; raw_reads_rna_mapping; raw_reads_prep_rseqc; }
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { raw_reads_fastqc; raw_reads_star; raw_reads_hisat2; raw_reads_tophat2; raw_reads_rna_mapping; raw_reads_prep_rseqc; }
}

// Header log info
log.info """=======================================================

 rnaseq : RNA-Seq workflow v${workflow.manifest.version}
======================================================="""
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Metadata']	= params.metadata
summary['Reads']        = params.reads
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Genome']       = params.genome
summary['Strandedness'] = params.stranded
if(params.aligner == 'star'){
  summary['Aligner'] = "STAR"
  if(params.star_index) summary['STAR Index'] = params.star_index
} else if(params.aligner == 'tophat2') {
  summary['Aligner'] = "Tophat2"
  if(params.bowtie2_index) summary['Tophat2 Index'] = params.bowtie2_index
} else if(params.aligner == 'hisat2') {
  summary['Aligner'] = "HISAT2"
  if(params.hisat2_index) summary['HISAT2 Index'] = params.hisat2_index
  if(params.splicesites) summary['Splice Sites'] = params.splicesites
}
summary['Counts'] = params.counts
if(params.gtf)                 summary['GTF Annotation']  = params.gtf
if(params.bed12)               summary['BED Annotation']  = params.bed12
summary['Save Intermeds'] = params.saveAlignedIntermediates ? 'Yes' : 'No'
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Config Profile'] = workflow.profile

if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


/*
 * FastQC
 */
process fastqc {
    tag "${prefix}"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when:
    !params.skip_qc && !params.skip_fastqc

    input:
    set val(name), file(reads) from raw_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    prefix = reads[0].toString() - ~/(_1)?(_2)?(_R1)?(_R2)?(.R1)?(.R2)?(_val_1)?(_val_2)?(trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/

    """
    fastqc -q $reads
    """
}


/*
 * rRNA mapping 
 */
process rRNA_mapping {
  tag "${prefix}"
  publishDir "${params.outdir}/rRNA_mapping", mode: 'copy',
      saveAs: {filename ->
          if (filename.indexOf("sorted.bam") > 0 &&  params.saveAlignedIntermediates) filename
	  else if (filename.indexOf("fastq.gz") > 0 &&  params.saveAlignedIntermediates) filename
	  else if (filename.indexOf(".log") > 0) "logs/$filename"
          else null
      }

  when:
    !params.skip_rrna

  input:
    set val(name), file(reads) from raw_reads_rna_mapping
    file annot from rrna_annot.collect()

  output:
    set val(name), file("${prefix}_norRNA_{1,2}.fastq.gz") into rrna_mapping_res
    file "*.log" into rrna_logs
    file "*sorted.bam"


  script:
  prefix = reads[0].toString() - ~/(_1)?(_2)?(_R1)?(_R2)?(.R1)?(.R2)?(_val_1)?(_val_2)?(trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/

  if (params.singleEnd) {
     """
     bowtie ${params.bowtie_opts} \\
     -p ${task.cpus} \\
     --un ${prefix}_norRNA.fastq \\
     --sam ${params.rrna} \\
     -1 ${reads} \\
     ${prefix}.sam  2> ${prefix}.log && \
     gzip -f ${prefix}_norRNA*.fastq 
     samtools view -@  ${task.cpus} -bS ${prefix}.sam > ${prefix}.bam  && \
     samtools sort \\
          ${prefix}.bam \\
          -@ ${task.cpus} \\
          -o ${prefix}.sorted.bam
    """
  } else {
     """
     bowtie ${params.bowtie_opts} \\
     -p ${task.cpus} \\
     --un ${prefix}_norRNA.fastq \\
     --sam ${params.rrna} \\
     -1 ${reads[0]} \\
     -2 ${reads[1]} \\
     ${prefix}.sam  2> ${prefix}.log && \
     gzip -f ${prefix}_norRNA_*.fastq 
     samtools view -@  ${task.cpus} -bS ${prefix}.sam > ${prefix}.bam  && \
     samtools sort \\
          ${prefix}.bam \\
          -@ ${task.cpus} \\
          -o ${prefix}.sorted.bam
     """
  }  
}


/*
 * Strandness
 */

process prep_rseqc {
  tag "${prefix}"
  input:
  set val(name), file(reads) from raw_reads_prep_rseqc

  output:
  file("${prefix}_subsample.bam") into bam_rseqc

  when:
  params.stranded == 'auto'

  script:
  prefix = reads[0].toString() - ~/(_1)?(_2)?(_R1)?(_R2)?(.R1)?(.R2)?(_val_1)?(_val_2)?(trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
  if (params.singleEnd) {
     """
     bowtie2 --fast --end-to-end --reorder \\
     -p ${task.cpus} \\
     -u ${params.n_check} \\
     -x ${params.bowtie2_index} \\
     -U ${reads} > ${prefix}_subsample.bam 
     """
  } else {
     """
     bowtie2 --fast --end-to-end --reorder \\
     -p ${task.cpus} \\
     -u ${params.n_check} \\
     -x ${params.bowtie2_index} \\
     -1 ${reads[0]} \\
     -2 ${reads[1]} > ${prefix}_subsample.bam
     """
  }
}

process rseqc {
  tag "${bam_rseqc.baseName - '_subsample'}"
  publishDir "${params.outdir}/rseqc" , mode: 'copy',
      saveAs: {filename ->
               if (filename.indexOf("bam_stat.txt") > 0) "bam_stat/$filename"
          else if (filename.indexOf("infer_experiment.txt") > 0) "infer_experiment/$filename"
          else filename
      }

  when:
  params.stranded == 'auto'

  input:
  file bam_rseqc
  file bed12 from bed_rseqc.collect()

  output:
  file "*.{txt,pdf,r,xls}" into rseqc_results
  file "${bam_rseqc.baseName}.ret_parserseq_output.txt" into parse_rseqc
   
  script:
  pathworkdir = "$baseDir/results/tmp" 
  """
  infer_experiment.py -i $bam_rseqc -r $bed12 > ${bam_rseqc.baseName}.infer_experiment.txt
  parse_rseq_output.sh ${bam_rseqc.baseName}.infer_experiment.txt > ${bam_rseqc.baseName}.ret_parserseq_output.txt
  mkdir -p $pathworkdir
  cp ${bam_rseqc.baseName}.ret_parserseq_output.txt $pathworkdir
  """
}

process parse_infer_experiment {
  tag "${name}"

  when:
  params.stranded == 'auto'

  input:
  file parse_rseqc

  output:
  val parse_res into rseqc_results_featureCounts, rseqc_results_HTseqCounts, rseqc_results_dupradar, rseqc_results_tophat, rseqc_results_table

  script:
  name = parse_rseqc[0].toString()
  pathworkdir = "$baseDir/results/tmp"  
  lines = new File("${pathworkdir}/${name}").findAll { it.startsWith('') }
  parse_res = lines[0] 
  
  """
  echo '${parse_res}'  > res.stranded.txt
  """
}

if (params.stranded != 'auto'){
   Channel.from( params.stranded )
     .into { rseqc_results_featureCounts; rseqc_results_HTseqCounts; rseqc_results_dupradar; rseqc_results_tophat; rseqc_results_table }
}


/*
 * Reads mapping
 */

// From nf-core
// Function that checks the alignment rate of the STAR output
// and returns true if the alignment passed and otherwise false
skipped_poor_alignment = []
def check_log(logs) {
  def percent_aligned = 0;
  logs.eachLine { line ->
    if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
      percent_aligned = matcher[0][1]
    }
  }
  logname = logs.getBaseName() - 'Log.final'
  if(percent_aligned.toFloat() <= '2'.toFloat() ){
      log.info "#################### VERY POOR ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($logname)    >> ${percent_aligned}% <<"
      skipped_poor_alignment << logname
      return false
  } else {
      log.info "          Passed alignment > star ($logname)   >> ${percent_aligned}% <<"
      return true
  }
}

// Update input channel
star_raw_reads = Channel.create()
if( params.rrna && !params.skip_rrna){
  star_raw_reads = rrna_mapping_res
}
else {  
  star_raw_reads = raw_reads_star
}

if(params.aligner == 'star'){
  hisat_stdout = Channel.from(false)
  process star {
    tag "$prefix"
    publishDir "${params.outdir}/STAR", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".bam") == -1) "logs/$filename"
            else if (!params.saveAlignedIntermediates) filename
            else if (params.saveAlignedIntermediates) filename
            else null
        }

    input:
    set val(name), file(reads) from star_raw_reads
    file index from star_index.collect()
    file gtf from gtf_star.collect()

    output:
    set file("*Log.final.out"), file ('*.bam') into star_aligned
    file "*.out" into alignment_logs
    file "*.out.tab" into star_log_counts
    file "*Log.out" into star_log
    file "${prefix}Aligned.sortedByCoord.out.bam.bai" into bam_index_rseqc
    file "*ReadsPerGene.out.tab" optional true into star_counts_to_merge, star_counts_to_r

    script:
    prefix = reads[0].toString() - ~/(_1)?(_2)?(_R1)?(_R2)?(.R1)?(.R2)?(_val_1)?(_val_2)?(trimmed)?(_norRNA)?(\.fq)?(\.fastq)?(\.gz)?$/
    def star_mem = task.memory ?: params.star_memory ?: false
    def avail_mem = star_mem ? "--limitBAMsortRAM ${star_mem.toBytes() - 100000000}" : ''
    def star_opt_add = params.counts == 'star' ? params.star_opts_counts : ''
    seqCenter = params.seqCenter ? "--outSAMattrRGline ID:$prefix 'CN:$params.seqCenter'" : ''
    """
    STAR --genomeDir $index \\
         --sjdbGTFfile $gtf \\
         --readFilesIn $reads  \\
         --runThreadN ${task.cpus} \\
         --runMode alignReads \\
         --outSAMtype BAM SortedByCoordinate  \\
         --readFilesCommand zcat \\
         --runDirPerm All_RWX \\
         --outTmpDir /local/scratch/rnaseq_\$(date +%d%s%S) \\
         --outFileNamePrefix $prefix  \\
         --outSAMattrRGline ID:$prefix SM:$prefix LB:Illumina PL:Illumina  \\
         ${params.star_opts} ${star_opt_add} 
            
    samtools index ${prefix}Aligned.sortedByCoord.out.bam
    """
    }

    // Filter removes all 'aligned' channels that fail the check
    star_aligned
        .filter { logs, bams -> check_log(logs) }
        .flatMap {  logs, bams -> bams }
    .into { bam_count; bam_preseq; bam_markduplicates; bam_featurecounts; bam_HTseqCounts; bam_read_dist }
}


// Update HiSat2 channel
hisat2_raw_reads = Channel.create()
if( params.rrna && !params.skip_rrna ){
    hisat2_raw_reads = rrna_mapping_res
}
else {
    hisat2_raw_reads = raw_reads_hisat2 
}

if(params.aligner == 'hisat2'){
  star_log = Channel.from(false)
  process hisat2Align {
    tag "$prefix"
    publishDir "${params.outdir}/HISAT2", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".hisat2_summary.txt") > 0) "logs/$filename"
            else if (!params.saveAlignedIntermediates) filename
            else if (params.saveAlignedIntermediates) filename
            else null
        }

    input:
    set val(name), file(reads) from hisat2_raw_reads 
    file hs2_indices from hs2_indices.collect()
    file alignment_splicesites from alignment_splicesites.collect()

    output:
    file "${prefix}.bam" into hisat2_bam
    file "${prefix}.hisat2_summary.txt" into alignment_logs

    script:
    index_base = hs2_indices[0].toString() - ~/.\d.ht2/
    prefix = reads[0].toString() - ~/(_1)?(_2)?(_R1)?(_R2)?(.R1)?(.R2)?(_val_1)?(_val_2)?(trimmed)?(_norRNA)?(\.fq)?(\.fastq)?(\.gz)?$/
    seqCenter = params.seqCenter ? "--rg-id ${prefix} --rg CN:${params.seqCenter.replaceAll('\\s','_')}" : ''
    def rnastrandness = ''
    if (params.stranded=='yes'){
        rnastrandness = params.singleEnd ? '--rna-strandness F' : '--rna-strandness FR'
    } else if (params.stranded=='reverse'){
        rnastrandness = params.singleEnd ? '--rna-strandness R' : '--rna-strandness RF'
    }
    if (params.singleEnd) {
    """
    hisat2 -x $index_base \\
           -U $reads \\
           $rnastrandness \\
           --known-splicesite-infile $alignment_splicesites \\
           -p ${task.cpus} \\
           --met-stderr \\
           --new-summary \\
           --summary-file ${prefix}.hisat2_summary.txt $seqCenter \\
           | samtools view -bS -F 4 -F 256 - > ${prefix}.bam
    """
    } else {
    """
    hisat2 -x $index_base \\
           -1 ${reads[0]} \\
           -2 ${reads[1]} \\
           $rnastrandness \\
           --known-splicesite-infile $alignment_splicesites \\
           --no-mixed \\
           --no-discordant \\
           -p ${task.cpus} \\
           --met-stderr \\
           --new-summary \\
           --summary-file ${prefix}.hisat2_summary.txt $seqCenter \\
           | samtools view -bS -F 4 -F 8 -F 256 - > ${prefix}.bam
     """
    }
  }

  process hisat2_sortOutput {
      tag "${hisat2_bam.baseName}"
      publishDir "${params.outdir}/HISAT2", mode: 'copy',
          saveAs: { filename ->
              if (!params.saveAlignedIntermediates) filename
              else if (params.saveAlignedIntermediates) "aligned_sorted/$filename"
              else null
          }

      input:
      file hisat2_bam

      output:
      file "${hisat2_bam.baseName}.sorted.bam" into bam_count, bam_preseq, bam_markduplicates, bam_featurecounts, bam_HTseqCounts, bam_read_dist 
      file "${hisat2_bam.baseName}.sorted.bam.bai" into bam_index_rseqc
 
      script:
      def avail_mem = task.memory ? "-m ${task.memory.toBytes() / task.cpus}" : ''
      """
      samtools sort \\
          $hisat2_bam \\
          -@ ${task.cpus} $avail_mem \\
          -o ${hisat2_bam.baseName}.sorted.bam
      samtools index ${hisat2_bam.baseName}.sorted.bam
      """
  }
}

// Update channel for TopHat2
tophat2_raw_reads = Channel.create()
if( params.rrna && !params.skip_rrna ){
    tophat2_raw_reads = rrna_mapping_res
}
else {
    tophat2_raw_reads = raw_reads_tophat2
}

if(params.aligner == 'tophat2'){
 process tophat2 {
  tag "${name}"
  publishDir "${params.outdir}/tophat2", mode: 'copy'

  input:
    set val(name), file(reads) from tophat2_raw_reads 
    file "tophat2" from tophat2_indices.collect()
    file gtf from gtf_tophat.collect()
    val parse_res from rseqc_results_tophat

  output:
    file "${prefix}.bam" into bam_count, bam_preseq, bam_markduplicates, bam_featurecounts, bam_HTseqCounts, bam_read_dist
    file "${prefix}.align_summary.txt" into alignment_logs

  script:
    prefix = reads[0].toString() - ~/(_1)?(_2)?(_R1)?(_R2)?(.R1)?(.R2)?(_val_1)?(_val_2)?(trimmed)?(_norRNA)?(\.fq)?(\.fastq)?(\.gz)?$/
    def avail_mem = task.memory ? "-m ${task.memory.toBytes() / task.cpus}" : ''
    def stranded_opt = '--library-type fr-unstranded'
    if (parse_res == 'yes'){
        stranded_opt = '--library-type fr-secondstrand'
    } else if ((parse_res == 'reverse')){
        stranded_opt = '--library-type fr-firststrand'
    }
    def out = './mapping'
    def sample = "--rg-id ${name} --rg-sample ${name} --rg-library Illumina --rg-platform Illumina --rg-platform-unit ${name}"
    """
    mkdir -p ${out}
    tophat2 -p ${task.cpus} \\
    ${sample} \\
    ${params.tophat2_opts} \\
    --GTF $gtf \\
    ${stranded_opt} \\
    -o ${out} \\
    ${params.bowtie2_index} \\
    ${reads} && \\
    mv ${out}/accepted_hits.bam ./${prefix}.bam
    mv ${out}/align_summary.txt ./${prefix}.align_summary.txt
    """
 }
}


/*
 * Saturation Curves
 */
process preseq {
  tag "${bam_preseq.baseName - 'Aligned.sortedByCoord.out'}"
  publishDir "${params.outdir}/preseq", mode: 'copy'

  when:
  !params.skip_qc && !params.skip_preseq

  input:
  file bam_preseq

  output:
  file "${bam_preseq.baseName}.ccurve.txt" into preseq_results

  script:
  """
  preseq lc_extrap -v -B $bam_preseq -o ${bam_preseq.baseName}.ccurve.txt
  """
}


/*
 * Duplicates
 */
process markDuplicates {
  tag "${bam.baseName - 'Aligned.sortedByCoord.out'}"
  publishDir "${params.outdir}/markDuplicates", mode: 'copy',
      saveAs: {filename -> filename.indexOf("_metrics.txt") > 0 ? "metrics/$filename" : "$filename"}

  when:
  !params.skip_qc && !params.skip_dupradar

  input:
  file bam from bam_markduplicates

  output:
  file "${bam.baseName}.markDups.bam" into bam_md
  file "${bam.baseName}.markDups_metrics.txt" into picard_results
  file "${bam.baseName}.markDups.bam.bai"

  script:
  if( !task.memory ){
    log.info "[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this."
    avail_mem = 3
  } else {
    avail_mem = task.memory.toGiga()
  }
  """
  picard -Xmx${avail_mem}g MarkDuplicates \\
      INPUT=$bam \\
      OUTPUT=${bam.baseName}.markDups.bam \\
      METRICS_FILE=${bam.baseName}.markDups_metrics.txt \\
      REMOVE_DUPLICATES=false \\
      ASSUME_SORTED=true \\
      PROGRAM_RECORD_ID='null' \\
      VALIDATION_STRINGENCY=LENIENT
  samtools index ${bam.baseName}.markDups.bam
  """
}

process dupradar {
  tag "${bam_md.baseName - 'Aligned.sortedByCoord.out.markDups'}"
  publishDir "${params.outdir}/dupradar", mode: 'copy',
      saveAs: {filename ->
          if (filename.indexOf("_duprateExpDens.pdf") > 0) "scatter_plots/$filename"
          else if (filename.indexOf("_duprateExpBoxplot.pdf") > 0) "box_plots/$filename"
          else if (filename.indexOf("_expressionHist.pdf") > 0) "histograms/$filename"
          else if (filename.indexOf("_dupMatrix.txt") > 0) "gene_data/$filename"
          else if (filename.indexOf("_duprateExpDensCurve.txt") > 0) "scatter_curve_data/$filename"
          else if (filename.indexOf("_intercept_slope.txt") > 0) "intercepts_slopes/$filename"
          else "$filename"
      }

  when:
    !params.skip_qc && !params.skip_dupradar

  input:
  file bam_md
  file gtf from gtf_dupradar.collect()
  val parse_res from rseqc_results_dupradar

  output:
  file "*.{pdf,txt}" into dupradar_results

  script: // This script is bundled with the pipeline, in nfcore/rnaseq/bin/
  def dupradar_direction = 0
  if (parse_res == 'yes'){
      dupradar_direction = 1
  } else if ((parse_res == 'reverse')){
      dupradar_direction = 2
  }
  def paired = params.singleEnd ? 'single' :  'paired'
  """
  dupRadar.r $bam_md $gtf $dupradar_direction $paired ${task.cpus}
  """
}

/*
 * Counts
 */

process featureCounts {
  tag "${bam_featurecounts.baseName - 'Aligned.sortedByCoord.out'}"
  publishDir "${params.outdir}/counts", mode: 'copy',
      saveAs: {filename ->
          if (filename.indexOf("_gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$filename"
          else if (filename.indexOf("_gene.featureCounts.txt") > 0) "gene_counts/$filename"
          else "$filename"
      }

  when:
  params.counts == 'featureCounts'

  input:
  file bam_featurecounts
  file gtf from gtf_featureCounts.collect()
  val parse_res from rseqc_results_featureCounts

  output:
  file "${bam_featurecounts.baseName}_gene.featureCounts.txt" into featureCounts_counts_to_merge, featureCounts_counts_to_r
  file "${bam_featurecounts.baseName}_gene.featureCounts.txt.summary" into featureCounts_logs

  script:
  def featureCounts_direction = 0
  if (parse_res == 'yes'){
      featureCounts_direction = 1
  } else if ((parse_res == 'reverse')){
      featureCounts_direction = 2
  }

  // Try to get real sample name
  sample_name = bam_featurecounts.baseName - 'Aligned.sortedByCoord.out'
  """
  featureCounts ${params.featurecounts_opts} -T ${task.cpus} -a $gtf -o ${bam_featurecounts.baseName}_gene.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
  """
}

process HTseqCounts {
  tag "${bam_HTseqCounts.baseName - 'Aligned.sortedByCoord.out'}"
  publishDir "${params.outdir}/counts", mode: 'copy',
      saveAs: {filename ->
          if (filename.indexOf("_gene.HTseqCounts.txt.summary") > 0) "gene_count_summaries/$filename"
          else if (filename.indexOf("_gene.HTseqCounts.txt") > 0) "gene_counts/$filename"
          else "$filename"
      }
  when:
  params.counts == 'HTseqCounts'

  input:
  file bam_HTseqCounts
  file gtf from gtf_HTseqCounts.collect()
  val parse_res from  rseqc_results_HTseqCounts

  output: 
  file "*_counts.csv" into htseq_counts_to_merge, htseq_counts_to_r

  script:
  def stranded_opt = '-s no' 
  if (parse_res == 'yes'){
      stranded_opt= '-s yes'
  } else if ((parse_res == 'reverse')){
      stranded_opt= '-s reverse'
  }

  // Try to get real sample name
  sample_name = bam_HTseqCounts.baseName - 'Aligned.sortedByCoord.out'
  """
  htseq-count ${params.htseq_opts} $stranded_opt $bam_HTseqCounts $gtf > ${sample_name}_counts.csv
  """
}


counts_to_merge = Channel.create()
counts_to_r = Channel.create()
if( params.counts == 'featureCounts' ){
    counts_to_merge = featureCounts_counts_to_merge
    counts_to_r = featureCounts_counts_to_r
} else if (params.counts == 'HTseqCounts'){
    counts_to_merge = htseq_counts_to_merge
    counts_to_r = htseq_counts_to_r	
}else if (params.counts == 'star'){
    counts_to_merge = star_counts_to_merge
    counts_to_r = star_counts_to_r
}

process merge_counts {
  publishDir "${params.outdir}/counts", mode: 'copy'

  input:
  file input_counts from counts_to_merge.collect()
  file gtf from gtf_table.collect()
  val parse_res from rseqc_results_table.collect()

  output:
  file 'tablecounts_raw.csv' into raw_counts
  file 'tablecounts_tpm.csv' into tpm_counts

  script:
  """
  echo -e ${input_counts} | tr " " "\n" > listofcounts.tsv
  echo -e ${parse_res} | sed -e "s/\\[//" -e "s/\\]//" -e "s/,//" | tr " " "\n" > listofstrandness.tsv
  makeCountTable.r listofcounts.tsv ${gtf} ${params.counts} listofstrandness.tsv
  """
}

counts_logs = Channel.create()
if( params.counts == 'featureCounts' ){
    counts_logs = featureCounts_logs
} else if (params.counts == 'HTseqCounts'){
    counts_logs = HTSeqCounts_logs
}else if (params.counts == 'star'){
    counts_logs = star_log_counts
}

/*
 * Reads distribution
 */

process read_distribution {
  tag "${bam_read_dist.baseName - '_subsample'}"
  publishDir "${params.outdir}/read_distribution" , mode: 'copy'

  when:
  !params.skip_read_dist

  input:
  file bam_read_dist
  file bed12 from bed_read_dist.collect()

  output:
  file "*.txt" into read_dist_results

  script:
  """
  read_distribution.py -i ${bam_read_dist} -r ${bed12} > ${bam_read_dist.baseName}.read_distribution.txt
  """
}


/*
 * Exploratory analysis
 */

process exploratory_analysis {
  publishDir "${params.outdir}/exploratory_analysis", mode: 'copy'

  when:
  !params.skip_expan

  input:
  file table_raw from raw_counts.collect()
  file table_tpm from tpm_counts.collect()
  val num_sample from counts_to_r.count()
  file pca_header from ch_pca_header
  file heatmap_header from ch_heatmap_header

  output:
  file "*.{txt,pdf,csv}" into exploratory_analysis_results

  when:
  num_sample > 1

  script:
  """
  exploratory_analysis.r ${table_raw}
  cat $pca_header deseq2_pca_coords_mqc.csv >> tmp_file
  mv tmp_file deseq2_pca_coords_mqc.csv 
  cat $heatmap_header vst_sample_cor_mqc.csv >> tmp_file
  mv tmp_file vst_sample_cor_mqc.csv
  """
}

/*
 * MultiQC
 */

process get_software_versions {
  output:
  file 'software_versions_mqc.yaml' into software_versions_yaml

  script:
  """
  echo $workflow.manifest.version &> v_rnaseq.txt
  echo $workflow.nextflow.version &> v_nextflow.txt
  fastqc --version &> v_fastqc.txt
  STAR --version &> v_star.txt
  tophat2 --version &> v_tophat2.txt
  hisat2 --version &> v_hisat2.txt
  preseq &> v_preseq.txt
  infer_experiment.py --version &> v_rseqc.txt
  read_duplication.py --version &> v_read_duplication.txt
  featureCounts -v &> v_featurecounts.txt
  htseq-count -h | grep version  &> v_htseq.txt
  picard MarkDuplicates --version &> v_markduplicates.txt  || true
  samtools --version &> v_samtools.txt
  multiqc --version &> v_multiqc.txt
  scrape_software_versions.py &> software_versions_mqc.yaml
  """
}

process workflow_summary_mqc {
  when:
  !params.skip_multiqc

  output:
  file 'workflow_summary_mqc.yaml' into workflow_summary_yaml

  exec:
  def yaml_file = task.workDir.resolve('workflow_summary_mqc.yaml')
  yaml_file.text  = """
  id: 'summary'
  description: " - this information is collected when the pipeline is started."
  section_name: 'Workflow Summary'
  section_href: 'https://gitlab.curie.fr/rnaseq'
  plot_type: 'html'
  data: |
      <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
      </dl>
  """.stripIndent()
}

process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    file metadata from ch_metadata.ifEmpty([])
    file multiqc_config from ch_multiqc_config    
    //file multiqc_custom_logo from ch_multiqc_logo
    file (fastqc:'fastqc/*') from fastqc_results.collect().ifEmpty([])
    file ('rrna/*') from rrna_logs.collect().ifEmpty([])
    file ('alignment/*') from alignment_logs.collect()
    file ('rseqc/*') from rseqc_results.collect().ifEmpty([])
    file ('rseqc/*') from read_dist_results.collect().ifEmpty([])
    file ('preseq/*') from preseq_results.collect().ifEmpty([])
    file ('dupradar/*') from dupradar_results.collect().ifEmpty([])
    file ('picard/*') from picard_results.collect().ifEmpty([])	
    file ('counts/*') from counts_logs.collect()
    file ('exploratory_analysis_results/*') from exploratory_analysis_results.collect().ifEmpty([]) // If the Edge-R is not run create an Empty array
    file ('software_versions/*') from software_versions_yaml.collect()
    file ('workflow_summary/*') from workflow_summary_yaml.collect()

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    makemeta = params.metadata ? true : false
    isPE = params.singleEnd ? 0 : 1
    
    modules_list = "-m custom_content -m picard -m preseq -m rseqc -m bowtie1 -m hisat2 -m star -m tophat -m cutadapt -m fastqc"
    modules_list = params.counts == 'featureCounts' ? "${modules_list} -m featureCounts" : "${modules_list}"  
    modules_list = params.counts == 'HTseqCounts' ? "${modules_list} -m HTSeq" : "${modules_list}"  
 
    if ( makemeta ){
      """
      metadata2multiqc.py $metadata > multiqc-config-metadata.yaml
      stats2multiqc.sh ${params.aligner} ${isPE}
      multiqc . -f $rtitle $rfilename -c $multiqc_config -c multiqc-config-metadata.yaml $modules_list
      """    
    }else{
      """	
      stats2multiqc.sh ${params.aligner} ${isPE}
      multiqc . -f $rtitle $rfilename -c $multiqc_config $modules_list
      """
    }
}


/*
 * Sub-routine
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}

workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[rnaseq] Successful: $workflow.runName"
    if(skipped_poor_alignment.size() > 0){
        subject = "[rnaseq] Partially Successful (${skipped_poor_alignment.size()} skipped): $workflow.runName"
    }
    if(!workflow.success){
      subject = "[nfcore/rnaseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['skipped_poor_alignment'] = skipped_poor_alignment

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success && !params.skip_multiqc) {
            mqc_report = multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList){
                log.warn "[nfcore/rnaseq] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
                }
        }
    } catch (all) {
        log.warn "[nfcore/rnaseq] Could not attach MultiQC report to summary email"
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nfcore/rnaseq] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nfcore/rnaseq] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    if(skipped_poor_alignment.size() > 0){
        log.info "[nfcore/rnaseq] WARNING - ${skipped_poor_alignment.size()} samples skipped due to poor alignment scores!"
    }

    log.info "[nfcore/rnaseq] Pipeline Complete"

    if(!workflow.success){
        if( workflow.profile == 'standard'){
            if ( "hostname".execute().text.contains('.uppmax.uu.se') ) {
                log.error "====================================================\n" +
                        "  WARNING! You are running with the default 'standard'\n" +
                        "  pipeline config profile, which runs on the head node\n" +
                        "  and assumes all software is on the PATH.\n" +
                        "  This is probably why everything broke.\n" +
                        "  Please use `-profile uppmax` to run on UPPMAX clusters.\n" +
                        "============================================================"
            }
        }
    }
}
