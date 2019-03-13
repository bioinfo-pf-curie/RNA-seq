#!/usr/bin/env nextflow



/************************
 * CONFIG
 **************************/

/*
genomes.config have been updated with path to local annotation (Mouse + Human)
Note that these paths are not always available for all organisms
*/

/*
base.config contains the cluster parameters for each tool
*/

/*
curie.config contains internal configuration values
Note that this file is not expected to be in this repository. Will move to nf-config repository
*/

/*
How to set up the parameters for all tools ?
*/

/*
How to set up the list of containers if we are using singularity ?
*/


/************************
 * CURIE RNA-SEQ WORKFLOW
 **************************/


/*
 * STEP 1 - FastQC - nf-core
*/
process fastqc {
   when:
     !params.skip_fastqc

   input:
     set val(name), file(reads) from raw_reads_fastqc
     
   output:
     file "*_fastqc.{zip,html}" into fastqc_results
}


/*
 * STEP 2 - RSeQC analysis
 * infer_experiment.py

/*
IMPORTANT
=========
We will need to translate the rseqc results and to use this results during the mapping
see the parse_rseqc_output() function at https://gitlab.curie.fr/data-analysis/RNA-seq/blob/devel/scripts/rna.inc.sh
========
*/

process rseqc {
  when:
    !params.skip_rseqc && !params.strandness

  input:
    file bam_rseqc
    file index from bam_index_rseqc
    file bed12 from bed_rseqc.collect()

  output:
    file "*.{txt}" into rseqc_results
}


/*
 * STEP 3 - rRNA mapping - NEW !!!
*/
process rRNA_mapping {
  when:
    rrna annotation are available 

  input:
    raw reads

  output:
    fastq file
}



/*
 * STEP 4 - align with STAR - nf-core
*/
process star {
  input:
    file reads from raw_reads
    file index from star_index.collect()
    file gtf from gtf_star.collect()
    file wherearemyfiles from ch_where_star.collect()

  output:
    set file("*Log.final.out"), file ('*.bam') into star_aligned
    file "*.out" into alignment_logs
    file "*SJ.out.tab"
    file "*Log.out" into star_log
    file "where_are_my_files.txt"
    file "${prefix}Aligned.sortedByCoord.out.bam.bai" into bam_index_rseqc, bam_index_genebody
} 


/*
 * STEP 4 - align with HISAT2 - nf-core
*/
process hisat2 {
  input:
    file reads from raw_reads
    file hs2_indices from hs2_indices.collect()
    file alignment_splicesites from alignment_splicesites.collect()
    file wherearemyfiles from ch_where_hisat2.collect()

  output:
    file "${prefix}.bam" 
    file "${prefix}.hisat2_summary.txt" into alignment_logs
    file "where_are_my_files.txt"
}

/*
 * STEP 4 - align with TOPHAT2 - NEW !!!
 */

process tophat2 {
  input:
    file reads from raw_reads
    file hs2_indices from hs2_indices.collect()
    file alignment_splicesites from alignment_splicesites.collect()
    file wherearemyfiles from ch_where_hisat2.collect()

  output:
    file "${prefix}.bam" 
    /*TO CHECK*/
}

/*
 * STEP 5 Mark duplicates - nf-core
 */

process markDuplicates {

  when:
    !params.skip_dupradar

  input:
    file bam from bam_markduplicates

  output:
    file "${bam.baseName}.markDups.bam" into bam_md
    file "${bam.baseName}.markDups_metrics.txt" into picard_results
    file "${bam.baseName}.markDups.bam.bai"
}

/*
 * STEP 5 - preseq analysis - nf-core
 */

process preseq {
  when:
    !params.skip_qc && !params.skip_preseq

  input:
    file bam_preseq

  output:
    file "${bam_preseq.baseName}.ccurve.txt" into preseq_results
}


/*
 * STEP 7 - dupRadar - nf-core
 */

process dupradar {

  when:
    !params.skip_qc && !params.skip_dupradar

  input:
    file bam_md
    file gtf from gtf_dupradar.collect()

  output:
    file "*.{pdf,txt}" into dupradar_results
}

/*
 * STEP 8 Feature counts - nf-core
 */

process featureCounts {
  input:
    file bam_featurecounts
    file gtf from gtf_featureCounts.collect()
  
  output:
    file "${bam_featurecounts.baseName}_gene.featureCounts.txt" into geneCounts, featureCounts_to_merge
    file "${bam_featurecounts.baseName}_gene.featureCounts.txt.summary" into featureCounts_logs
}

/*
 * STEP 8 HTSeq counts - NEW !
 */

process HTseqCounts {
  input:
    file bam_featurecounts
    file gtf from gtf_featureCounts.collect()
  
  output:
    file "${bam_featurecounts.baseName}_gene.featureCounts.txt" into geneCounts, featureCounts_to_merge
    file "${bam_featurecounts.baseName}_gene.featureCounts.txt.summary" into featureCounts_logs
}


/*
 * STEP 8 STAR counts - NEW !

IMPORTANT
=========
Just add the option --quantMode GeneCounts of STAR is used for the mapping
========

 */

